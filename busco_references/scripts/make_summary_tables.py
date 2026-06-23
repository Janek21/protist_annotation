#!/usr/bin/env python3
"""Build summary tables from the BUSCO/counts outputs under summary/.

Produces three TSVs in the summary directory:
  - counts_summary.tsv   : gene/transcript counts per species
  - busco_summary.tsv    : BUSCO completeness (custom lineage vs eukaryote) per species/mode
  - general_summary.tsv  : merged view built from the two tables above, one row per species
"""
import argparse
import json
import os
import re
import sys
from glob import glob

import pandas as pd

MODES = ("gen", "prot")
C_RE = re.compile(r"C:\s*([\d.]+)")

# derived metrics appended to the counts table (same definitions as the
# annotator family's _metrics.tsv; see isoquant_annotator/derived_metrics.md).
# For the reference annotation the ORF-calling input is the annotation itself,
# so transcriptome_transcripts == transcript_count and coding_fraction is the
# share of transcript models carrying a CDS protein.
DERIVED_COLUMNS = [
    "coding_transcripts",
    "transcriptome_transcripts",
    "gene_density_per_mb",
    "transcript_density_per_mb",
    "isoforms_per_gene",
    "coding_fraction",
]


def _num(x):
    """Float or None (None for missing/blank/NaN), so metric guards never crash."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def add_derived_metrics(df):
    """Append transcriptome_transcripts + the four derived metrics, NA-safe."""
    def derive(row):
        gc = _num(row.get("gene_count"))
        tc = _num(row.get("transcript_count"))
        gs = _num(row.get("genome_size_bp"))
        cod = _num(row.get("coding_transcripts"))
        ntx = tc  # reference transcript models are the ORF-calling input
        return pd.Series({
            "transcriptome_transcripts": "NA" if ntx is None else row.get("transcript_count"),
            "gene_density_per_mb": f"{gc / (gs / 1e6):.2f}" if gc is not None and gs and gs > 0 else "NA",
            "transcript_density_per_mb": f"{tc / (gs / 1e6):.2f}" if tc is not None and gs and gs > 0 else "NA",
            "isoforms_per_gene": f"{tc / gc:.3f}" if tc is not None and gc and gc > 0 else "NA",
            "coding_fraction": f"{cod / ntx:.4f}" if cod is not None and ntx and ntx > 0 else "NA",
        })
    if df.empty:
        for col in DERIVED_COLUMNS:
            df[col] = pd.Series(dtype=str)
        return df
    return pd.concat([df, df.apply(derive, axis=1)], axis=1)


def parse_key(stem, suffix):
    """Split '<species>_<taxonID>_<mode><suffix>' from the right.

    Species names contain underscores, so we peel known fields off the end.
    Returns (species, taxonID, mode) or None if it does not match.
    """
    if not stem.endswith(suffix):
        return None
    parts = stem[: -len(suffix)].rsplit("_", 2)
    if len(parts) != 3:
        return None
    species, taxon, mode = parts
    if mode not in MODES:
        return None
    return species, taxon, mode


def read_score(path):
    """Return (lineage_name, completeness) from a BUSCO short_summary json.

    completeness is the Complete (C) number pulled from the one-line summary,
    e.g. 95.0 from 'C:95.0%[S:..,D:..],F:..,M:..,n:..'.
    """
    with open(path) as fh:
        data = json.load(fh)
    lineage = data.get("lineage_dataset", {}).get("name", "NA")
    m = C_RE.search(data.get("results", {}).get("one_line_summary", ""))
    completeness = float(m.group(1)) if m else float("nan")
    return lineage, completeness


def collect_counts(summary_dir):
    """DataFrame: species, taxonID, counts, and the derived metrics."""
    rows = {}
    for kind, col in (("gc", "gene_count"), ("tc", "transcript_count"),
                      ("gs", "genome_size_bp"), ("ct", "coding_transcripts")):
        for path in glob(os.path.join(summary_dir, "counts", f"*_{kind}.txt")):
            species, _, taxon = os.path.basename(path)[: -len(f"_{kind}.txt")].rpartition("_")
            if not species:
                continue
            with open(path) as fh:
                rows.setdefault((species, taxon), {})[col] = fh.read().strip()
    df = pd.DataFrame(
        [{"species": sp, "taxonID": tax, **v} for (sp, tax), v in rows.items()],
        columns=["species", "taxonID", "gene_count", "transcript_count",
                 "genome_size_bp", "coding_transcripts"],
    )
    df = add_derived_metrics(df)
    return df.sort_values(["species", "taxonID"], ignore_index=True)


def collect_busco(summary_dir):
    """DataFrame: species, taxonID, mode, lineage_used, busco_custom_lineage, busco_eukaryote."""
    rows = {}
    for sub, suffix, score_col, take_lineage in (
        ("busco_lineage", "_Lbusco.json", "busco_custom_lineage", True),
        ("busco_eukaryote", "_Ebusco.json", "busco_eukaryote", False),
    ):
        for path in glob(os.path.join(summary_dir, sub, f"*{suffix}")):
            parsed = parse_key(os.path.basename(path), suffix)
            if not parsed:
                continue
            species, taxon, mode = parsed
            lineage, completeness = read_score(path)
            row = rows.setdefault((species, taxon, mode), {})
            row[score_col] = completeness
            if take_lineage:
                row["lineage_used"] = lineage
    df = pd.DataFrame(
        [{"species": sp, "taxonID": tax, "mode": mode, **v}
         for (sp, tax, mode), v in rows.items()],
        columns=["species", "taxonID", "mode", "lineage_used",
                 "busco_custom_lineage", "busco_eukaryote"],
    )
    return df.sort_values(["species", "taxonID", "mode"], ignore_index=True)


def write_tsv(df, path):
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    print(f"wrote {path} ({len(df)} rows)")


_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)
from buscoPlot import run as _busco_plot


def run_busco_plot(glob_pattern, output_path):
    _busco_plot(glob_pattern, output_path)


def build_general(summary_dir):
    """Merge counts_summary.tsv and busco_summary.tsv into general_summary.tsv.

    Reads the two derived tables (not the source data) and pivots BUSCO rows
    by mode so each species has one row.
    """
    counts = pd.read_csv(os.path.join(summary_dir, "counts_summary.tsv"), sep="\t")
    busco = pd.read_csv(os.path.join(summary_dir, "busco_summary.tsv"), sep="\t")

    lineage = busco.groupby(["species", "taxonID"])["lineage_used"].first()
    wide = busco.pivot_table(
        index=["species", "taxonID"], columns="mode",
        values=["busco_custom_lineage", "busco_eukaryote"], aggfunc="first",
    )
    wide.columns = [f"{mode}_{val}" for val, mode in wide.columns]
    score_cols = [f"{m}_{v}" for v in ("busco_custom_lineage", "busco_eukaryote")
                  for m in MODES]
    wide = wide.reindex(columns=score_cols)

    general = (
        counts.merge(lineage, on=["species", "taxonID"], how="outer")
        .merge(wide, on=["species", "taxonID"], how="outer")
        .sort_values(["species", "taxonID"], ignore_index=True)
    )
    general = general[["species", "taxonID", "gene_count", "transcript_count",
                       "genome_size_bp"] + DERIVED_COLUMNS
                      + ["lineage_used"] + score_cols]
    write_tsv(general, os.path.join(summary_dir, "general_summary.tsv"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-s", "--summary-dir", default="summary",
                    help="summary directory (default: summary)")
    args = ap.parse_args()
    sd = args.summary_dir

    write_tsv(collect_counts(sd), os.path.join(sd, "counts_summary.tsv"))
    write_tsv(collect_busco(sd), os.path.join(sd, "busco_summary.tsv"))
    build_general(sd)

    lin_dir = os.path.join(sd, "busco_lineage")
    euk_dir = os.path.join(sd, "busco_eukaryote")
    plots = [
        (os.path.join(lin_dir, "*_prot_Lbusco.json"),
         os.path.join(lin_dir, "busco_prot_lineage_summary.png")),
        (os.path.join(euk_dir, "*_prot_Ebusco.json"),
         os.path.join(euk_dir, "busco_prot_eukaryote_summary.png")),
        (os.path.join(lin_dir, "*_gen_Lbusco.json"),
         os.path.join(lin_dir, "busco_gen_lineage_summary.png")),
        (os.path.join(euk_dir, "*_gen_Ebusco.json"),
         os.path.join(euk_dir, "busco_gen_eukaryote_summary.png")),
    ]
    for glob_pat, out_path in plots:
        run_busco_plot(glob_pat, out_path)


if __name__ == "__main__":
    main()
