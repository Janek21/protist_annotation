#!/usr/bin/env python3
"""
fetch_protist_rnaseq.py
=======================
Query the EBI ENA Portal API for RNA-seq runs across the protist supergroups and
attach the numeric NCBI lineage to each run. Choose which reads to keep with
--reads (default: long, i.e. the original longread_protists.tsv):

  long         long-read runs only            (ONT / PacBio)
  short        short-read runs only           (everything else)
  both         every RNA-seq run
  prefer-long  long-read runs, plus short-read runs ONLY for species that have
               no long-read run at all (i.e. if a species has both, its short
               reads are dropped). "Species" granularity is set by --group-by.

Convention follows Cobos-Bioinfo/Euka-Survey (db_builder/build_db/get_reads.py):
  - POST to the ENA Portal API `search` endpoint, result = read_run
  - query  = tax_tree(<taxid>) AND (library_strategy="rna-seq" OR library_strategy="fl-cdna")
  - long-read  == instrument_platform in {OXFORD_NANOPORE, PACBIO_SMRT}
  - short-read == any other platform (ILLUMINA, BGISEQ, ION_TORRENT, ...)

One query per supergroup pulls all rna-seq; the platform split is done
CLIENT-SIDE (Euka-Survey style). The numeric lineage (column 4) and the
species-rank grouping for prefer-long are resolved from NCBI taxonomy via ete3.

Output columns (tab-separated, no header):
  1  <Supergroup>.tsv:<run_accession>
  2  experiment_title          6  library_source        9  read_count
  3  tax_id                    7  library_strategy      10  first_public
  4  numeric NCBI lineage      8  instrument_platform
  5  scientific_name

Usage:
  pip install requests tenacity ete3
  python fetch_protist_rnaseq.py                                    # long -> longread_protists.tsv
  python fetch_protist_rnaseq.py --reads short                     # short
  python fetch_protist_rnaseq.py --reads both                      # both
  python fetch_protist_rnaseq.py --reads prefer-long               # prefer-long
  python fetch_protist_rnaseq.py --reads prefer-long --group-by tax_id
  python fetch_protist_rnaseq.py --taxids 3041 2763 33634          # specific tax IDs
  python fetch_protist_rnaseq.py --taxids-file my_taxids.txt       # one taxid per line
(First run downloads the NCBI taxdump for ete3, ~1-2 min.)

When --taxids or --taxids-file is given the hardcoded SUPERGROUPS list is ignored
and only the provided tax IDs are queried. Each taxid is used as its own group
label in the output (column 1 becomes "<taxid>.tsv:<run_accession>").
"""

from __future__ import annotations

import argparse
import sys
import time

import requests
from tenacity import retry, stop_after_attempt, wait_exponential

ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"

# Protist groups -> defining NCBI tax IDs. Every value is a root/clade node so
# one tax_tree() query captures all lineages beneath it (no within-clade gaps).
# Spans the eukaryote tree (Adl et al. 2019; Burki et al. 2020) EXCLUDING
# animals, land plants, and Fungi. Chytridiomycota and Microsporidia are true
# Fungi and are therefore omitted.
SUPERGROUPS: dict[str, int] = {
    # SAR
    "Stramenopiles":    33634,
    "Alveolata":        33630,
    "Rhizaria":         543769,
    # Archaeplastida (algal; Chlorophyta excludes land-plant Streptophyta)
    "Chlorophyta":      3041,
    "Rhodophyta":       2763,
    "Glaucophyta":      38254,   # NCBI: Glaucocystophyceae
    # Haptista / Cryptista / Telonemia
    "Haptista":         2608109,
    "Cryptophyceae":    3027,    # no NCBI "Cryptista" node; class is the broadest
    "Telonemia":        660925,
    # Discoba / Metamonada (former Excavata)
    "Discoba":          2611352,
    "Metamonada":       2611341,
    # Amoebozoa
    "Amoebozoa":        554915,
    # Unicellular Opisthokonta (non-animal, non-fungal)
    "Choanoflagellata": 28009,
    "Ichthyosporea":    127916,
    "Filasterea":       2687318,
    # Apusomonadida / Breviatea / CRuMs / Hemimastigophora
    "Apusomonadida":    2925400,
    "Breviatea":        1401294,
    "CRuMs":            2608240,
    "Hemimastigophora": 2489521,
}

# A run is "long" if its platform is one of these; otherwise it is "short".
LONG_READ_PLATFORMS = {"OXFORD_NANOPORE", "PACBIO_SMRT"}

# Default output filename per mode (used when -o/--out is not given).
DEFAULT_OUT = {
    "long":        "longread_protists.tsv",
    "short":       "shortread_protists.tsv",
    "both":        "allread_protists.tsv",
    "prefer-long": "preferlong_protists.tsv",
}

# ENA read_run fields we export, in the file's column order (lineage excluded;
# it is filled from NCBI). run_accession is requested but re-emitted with the
# supergroup-label prefix.
ENA_FIELDS = [
    "run_accession",
    "experiment_title",
    "tax_id",
    "scientific_name",
    "library_source",
    "library_strategy",
    "instrument_platform",
    "read_count",
    "first_public",
]


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=2, max=60))
def ena_rnaseq_runs(taxid: int) -> list[dict]:
    """All RNA-seq / FL-cDNA read_runs under `taxid` (limit=0 -> no cap), unfiltered by platform."""
    payload = {
        "result": "read_run",
        "query": f'tax_tree({taxid}) AND (library_strategy="rna-seq" OR library_strategy="fl-cdna")',
        "fields": ",".join(ENA_FIELDS),
        "format": "json",
        "limit": 0,
    }
    r = requests.post(
        ENA_BASE,
        data=payload,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=180,
    )
    r.raise_for_status()
    try:
        return r.json() or []
    except requests.exceptions.JSONDecodeError:
        return []


def _safe_taxid(rec: dict):
    try:
        return int(rec.get("tax_id"))
    except (TypeError, ValueError):
        return None


def _is_long(rec: dict) -> bool:
    return rec.get("instrument_platform", "") in LONG_READ_PLATFORMS


def species_with_long(records: list[dict], species_key_fn) -> set:
    """Set of species-keys that have at least one long-read run."""
    keys = set()
    for rec in records:
        if _is_long(rec):
            tid = _safe_taxid(rec)
            if tid is not None:
                keys.add(species_key_fn(tid))
    return keys


def keep_record(rec: dict, mode: str, long_species: set, species_key_fn) -> bool:
    """Decide whether a run belongs in the requested mode.

    For 'prefer-long', `long_species` is the precomputed set of species-keys
    that have >=1 long-read run, and short reads of those species are dropped.
    """
    is_long = _is_long(rec)
    if mode == "long":
        return is_long
    if mode == "short":
        return not is_long
    if mode == "both":
        return True
    # prefer-long
    if is_long:
        return True
    tid = _safe_taxid(rec)
    if tid is None:
        return True  # unknown species; row builder will skip if truly invalid
    return species_key_fn(tid) not in long_species


def make_taxonomy(group_by: str):
    """Return (lineage_of, species_key_of) backed by the local NCBI taxdump.

    group_by='species' collapses each tax_id to its species-rank ancestor (or
    itself if no species ancestor exists); group_by='tax_id' uses the tax_id as-is.
    """
    from ete3 import NCBITaxa

    ncbi = NCBITaxa()
    lin_cache: dict[int, str] = {}
    sp_cache: dict[int, int] = {}

    def lineage_of(taxid: int) -> str:
        if taxid not in lin_cache:
            try:
                ids = ncbi.get_lineage(taxid) or []
            except Exception:
                ids = []
            lin_cache[taxid] = ";".join(str(t) for t in ids)
        return lin_cache[taxid]

    def species_key_of(taxid: int) -> int:
        if group_by == "tax_id":
            return taxid
        if taxid not in sp_cache:
            key = taxid
            try:
                lin = ncbi.get_lineage(taxid) or []
                ranks = ncbi.get_rank(lin)
                for t in reversed(lin):
                    if ranks.get(t) == "species":
                        key = t
                        break
            except Exception:
                key = taxid
            sp_cache[taxid] = key
        return sp_cache[taxid]

    return lineage_of, species_key_of


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--reads", choices=("long", "short", "both", "prefer-long"), default="long",
        help="which RNA-seq reads to keep (default: long). 'prefer-long' keeps "
             "long reads and drops a species' short reads when that species also "
             "has long reads.",
    )
    ap.add_argument(
        "--group-by", choices=("species", "tax_id"), default="species",
        help="granularity of 'species' for --reads prefer-long: collapse to "
             "species rank (default) or use the exact tax_id.",
    )
    ap.add_argument("-o", "--out", default=None,
                    help="combined output TSV (default depends on --reads)")
    ap.add_argument("--per-group", action="store_true",
                    help="also write one <group>.tsv per group / taxid")
    taxid_src = ap.add_mutually_exclusive_group()
    taxid_src.add_argument(
        "--taxids", nargs="+", type=int, metavar="TAXID",
        help="one or more NCBI tax IDs to query instead of the built-in supergroups",
    )
    taxid_src.add_argument(
        "--taxids-file", metavar="FILE",
        help="file with one NCBI tax ID per line (comments with # are ignored)",
    )
    args = ap.parse_args()

    if args.taxids:
        groups = {str(t): t for t in args.taxids}
    elif args.taxids_file:
        groups = {}
        with open(args.taxids_file) as fh:
            for line in fh:
                line = line.split("#")[0].strip()
                if line:
                    groups[line] = int(line)
    else:
        groups = SUPERGROUPS

    out_path = args.out or DEFAULT_OUT[args.reads]
    lineage_of, species_key_of = make_taxonomy(args.group_by)

    combined: list[str] = []
    for group, taxid in groups.items():
        print(f"[ENA] {group} (tax_tree {taxid}) ...", file=sys.stderr)
        runs = ena_rnaseq_runs(taxid)

        long_species = (
            species_with_long(runs, species_key_of)
            if args.reads == "prefer-long" else set()
        )

        group_rows: list[str] = []
        for rec in runs:
            if not keep_record(rec, args.reads, long_species, species_key_of):
                continue
            tid = _safe_taxid(rec)
            if tid is None:
                continue
            row = [
                f"{group}.tsv:{rec.get('run_accession', '')}",
                rec.get("experiment_title", ""),
                str(tid),
                lineage_of(tid),
                rec.get("scientific_name", ""),
                rec.get("library_source", ""),
                rec.get("library_strategy", ""),
                rec.get("instrument_platform", ""),
                rec.get("read_count", ""),
                rec.get("first_public", ""),
            ]
            group_rows.append("\t".join(row))

        print(f"      kept {len(group_rows)} / {len(runs)} runs ({args.reads})",
              file=sys.stderr)
        if args.per_group:
            with open(f"{group}.tsv", "w") as fh:
                fh.write("\n".join(group_rows) + ("\n" if group_rows else ""))
        combined.extend(group_rows)
        time.sleep(0.3)  # be gentle with the API

    with open(out_path, "w") as fh:
        fh.write("\n".join(combined) + ("\n" if combined else ""))
    print(f"Wrote {len(combined)} rows -> {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
