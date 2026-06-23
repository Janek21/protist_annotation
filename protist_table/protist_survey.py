#!/usr/bin/env python3
"""
protist_survey.py  -  NCBI Protist Genomic Data Survey
=======================================================
Builds a five-column table of protist species from NCBI:

  Col 1 : Protist Group            (taxonomic clade)
  Col 2 : Species Name             (NCBI organism name)
  Col 3 : TaxID                    (NCBI Taxonomy identifier)
  Col 4 : Assembly Accession       (accession for any genome assembly;
                                    - if the species has no assembly)
  Col 5 : Annotation Accession     (GCF_… RefSeq / GCA_… GenBank /
                                    - assembly unannotated or no assembly)
  Col 6 : RNAseq Data              (short-read only / long-read only /
                                    short-read & long-read / None)
  Col 7 : GC%                       (genome GC content, from assembly_stats;
                                    present for annotated AND unannotated
                                    assemblies, blank only when no assembly)
  Col 8 : Long-read Experiments    (# distinct SRA experiments, long-read)
  Col 9 : Long-read Reads          (total reads across long-read runs)
  Col 10: Short-read Experiments   (# distinct SRA experiments, short-read)
  Col 11: Short-read Reads         (total reads across short-read runs)

  Tip: "assembly but no annotation" rows have a value in Assembly Accession
  and GC% but "-" in Annotation Accession. SRA-only rows (no genome) have "-"
  in both accession columns and a blank GC%. By DEFAULT only genome-backed
  species are listed (--assemblies-only); pass --no-assemblies-only to also
  include SRA-only species.

Experiment / read counting
---------------------------
  • "Experiments" = number of distinct SRA Experiment (SRX) accessions.
  • "Reads"       = sum of (spots + spots_with_mates) over the runs, i.e.
                    paired-end reads count as two. Reads are de-duplicated by
                    Run accession so nothing is counted twice.
  • Long-read AND short-read presence/counts are kept per exact taxid (strain):
    a strain assembly is credited only with runs deposited on its own taxid,
    never with a sibling strain's or the parent species' reads. (SRA-only rows,
    which have no assembly, roll their counts up to species rank.)
  • Long-read runinfo is fetched in full, so long-read counts are exact.
    Short-read runinfo is capped per group (--sra-run-cap, default 30000);
    when a group exceeds the cap, its short-read counts reflect that sample
    rather than the complete total.

Design features:

  • SRA-FIRST PASS - In addition to assembly-derived species, the program
    asks SRA directly which organisms have long-read and short-read RNAseq,
    and UNIONS them in. This captures transcriptome-only organisms that have
    no genome assembly (e.g. many cryptophytes / MMETSP taxa).

  • STRAIN-LEVEL RNASEQ - Both long-read and short-read RNAseq are attributed by
    EXACT taxid, so a strain assembly is flagged (and counted) only for the runs
    deposited on that same strain (e.g. 'Plasmodium falciparum 7G8' is NOT
    credited with reads deposited under 'P. falciparum' or strain '3D7'). Only
    SRA-only organisms, which have no assembly to pin a strain to, roll their
    presence/counts up to species rank.

Data sources
------------
  • NCBI Datasets API v2  – genome assemblies & annotation metadata
  • NCBI Entrez (SRA)     – RNAseq run inventory (runinfo)
  • NCBI Entrez (Taxonomy)– species-rank resolution
      short-read = Illumina or Ion Torrent  (RNA-Seq strategy)
      long-read  = Oxford Nanopore or PacBio (RNA-Seq or FL-cDNA strategy)

Usage
-----
  pip install requests pandas

  python protist_survey.py                       # full run
  python protist_survey.py --limit 10 --no-sra   # quick assembly-only test
  python protist_survey.py --api-key YOUR_KEY     # 10 req/s instead of 3
  python protist_survey.py --output protists.tsv --html
  python protist_survey.py --sra-run-cap 50000    # cap short-read runinfo / group
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from textwrap import dedent

import pandas as pd
import requests

# ──────────────────────────────────────────────────────────────────────────────
#  API endpoints
# ──────────────────────────────────────────────────────────────────────────────
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENTREZ_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ──────────────────────────────────────────────────────────────────────────────
#  Protist taxonomic groups  (display name → NCBI Taxonomy ID)
# ──────────────────────────────────────────────────────────────────────────────
# Every value is a root/clade node so one txid[Organism:exp] query captures all
# lineages beneath it (avoids the within-clade gaps of enumerating sub-classes,
# e.g. Blastocystis under Stramenopiles or Perkinsozoa under Alveolata). Spans
# the eukaryote tree (Adl et al. 2019; Burki et al. 2020) EXCLUDING animals,
# land plants, and Fungi (Chytridiomycota and Microsporidia are true Fungi).
PROTIST_GROUPS: dict[str, int] = {
    # ── SAR ─────────────────────────────────────────────────────────────────
    "Stramenopiles"      : 33634,
    "Alveolata"          : 33630,
    "Rhizaria"           : 543769,
    # ── Archaeplastida (algal; Chlorophyta excludes land-plant Streptophyta) ─
    "Chlorophyta"        : 3041,
    "Rhodophyta"         : 2763,
    "Glaucophyta"        : 38254,    # NCBI: Glaucocystophyceae
    # ── Haptista / Cryptista / Telonemia ─────────────────────────────────────
    "Haptista"           : 2608109,
    "Cryptophyceae"      : 3027,     # no NCBI "Cryptista" node; class is broadest
    "Telonemia"          : 660925,
    # ── Discoba / Metamonada (former Excavata) ───────────────────────────────
    "Discoba"            : 2611352,
    "Metamonada"         : 2611341,
    # ── Amoebozoa ─────────────────────────────────────────────────────────────
    "Amoebozoa"          : 554915,
    # ── Unicellular Opisthokonta (non-animal, non-fungal) ────────────────────
    "Choanoflagellata"   : 28009,
    "Ichthyosporea"      : 127916,
    "Filasterea"         : 2687318,
    # ── Apusomonadida / Breviatea / CRuMs / Hemimastigophora ─────────────────
    "Apusomonadida"      : 2925400,
    "Breviatea"          : 1401294,
    "CRuMs"              : 2608240,
    "Hemimastigophora"   : 2489521,
}

# ──────────────────────────────────────────────────────────────────────────────
#  Protist group → eukaryotic supergroup  (mirrors the sections above)
# ──────────────────────────────────────────────────────────────────────────────
SUPERGROUPS: dict[str, str] = {
    "Stramenopiles"    : "SAR",
    "Alveolata"        : "SAR",
    "Rhizaria"         : "SAR",
    "Chlorophyta"      : "Archaeplastida",
    "Rhodophyta"       : "Archaeplastida",
    "Glaucophyta"      : "Archaeplastida",
    "Haptista"         : "Haptista",
    "Cryptophyceae"    : "Cryptista",
    "Telonemia"        : "Telonemia",
    "Discoba"          : "Discoba",
    "Metamonada"       : "Metamonada",
    "Amoebozoa"        : "Amoebozoa",
    "Choanoflagellata" : "Opisthokonta",
    "Ichthyosporea"    : "Opisthokonta",
    "Filasterea"       : "Opisthokonta",
    "Apusomonadida"    : "Obazoa",
    "Breviatea"        : "Obazoa",
    "CRuMs"            : "CRuMs",
    "Hemimastigophora" : "Hemimastigophora",
}

# SRA Entrez platform sub-queries
_LONG_PLAT  = "(OXFORD_NANOPORE[Platform] OR PACBIO_SMRT[Platform])"
_SHORT_PLAT = "(ILLUMINA[Platform] OR ION_TORRENT[Platform])"


# ──────────────────────────────────────────────────────────────────────────────
#  HTTP helpers
# ──────────────────────────────────────────────────────────────────────────────

def make_session(api_key: str | None, email: str) -> tuple[requests.Session, float]:
    s = requests.Session()
    s.headers.update({
        "Accept"    : "application/json",
        "User-Agent": f"protist_survey/2.0 ({email})",
    })
    if api_key:
        s.headers.update({"api-key": api_key})
    delay = 0.12 if api_key else 0.40          # 10 vs 3 req/s
    return s, delay


def _request(session, url, params, delay, want_json, retries=4):
    """Shared rate-limited GET with backoff. Returns dict (json) or str (text)."""
    for attempt in range(retries):
        time.sleep(delay)
        try:
            r = session.get(url, params=params, timeout=180)
            if r.status_code == 429:
                wait = 2 ** (attempt + 2)
                print(f"\n  ⚡ Rate-limited – sleeping {wait}s …", file=sys.stderr)
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json() if want_json else r.text
        except requests.RequestException as exc:
            if attempt == retries - 1:
                print(f"\n  ⚠  Request failed: {exc}", file=sys.stderr)
                return {} if want_json else ""
            time.sleep(2 ** attempt)
    return {} if want_json else ""


def ncbi_get(session, url, params, delay):
    return _request(session, url, params, delay, want_json=True)


def ncbi_get_text(session, url, params, delay):
    return _request(session, url, params, delay, want_json=False)


# ──────────────────────────────────────────────────────────────────────────────
#  NCBI Datasets API v2 - genome assemblies
# ──────────────────────────────────────────────────────────────────────────────

def fetch_assemblies(session, taxid, limit, delay):
    url = f"{DATASETS_BASE}/genome/taxon/{taxid}/dataset_report"
    reports, page_token = [], None
    while True:
        remaining = (limit - len(reports)) if limit else 1000
        params = {
            "page_size"               : min(1000, remaining),
            "filters.exclude_atypical": "true",
            "returned_content"        : "COMPLETE",
        }
        if page_token:
            params["page_token"] = page_token
        data  = ncbi_get(session, url, params, delay)
        batch = data.get("reports", [])
        reports.extend(batch)
        if not batch:
            break
        if limit and len(reports) >= limit:
            reports = reports[:limit]
            break
        page_token = data.get("next_page_token")
        if not page_token:
            break
    return reports


def _annotation_accession(accession: str, ann_info: dict) -> str:
    if accession.startswith("GCF"):
        return accession
    if ann_info.get("name") or ann_info.get("release_date") or ann_info.get("stats"):
        return accession
    return "-"


def _gc_percent(report: dict):
    """Pull GC content (%) from assembly_stats; '' if unavailable."""
    stats = report.get("assembly_stats") or {}
    gc = stats.get("gc_percent")
    try:
        return round(float(gc), 2) if gc not in (None, "") else ""
    except (TypeError, ValueError):
        return ""


def parse_report(report: dict, group: str) -> dict | None:
    organism = report.get("organism", {})
    species  = organism.get("organism_name", "").strip()
    taxid    = organism.get("tax_id")
    if not species or taxid is None:
        return None
    accession = report.get("accession", "")
    ann_info  = report.get("annotation_info") or {}
    return {
        "_taxid"              : taxid,
        "_accession"          : accession,
        "Group"               : group,
        "Species"             : species,
        "TaxID"               : taxid,
        "Assembly Accession"  : accession or "-",
        "GC%"                 : _gc_percent(report),
        "Annotation Accession": _annotation_accession(accession, ann_info),
    }


# ──────────────────────────────────────────────────────────────────────────────
#  NCBI Entrez - SRA-first organism discovery (runinfo)
# ──────────────────────────────────────────────────────────────────────────────

def _to_int(s: str | None) -> int:
    """Parse an integer from a runinfo cell; 0 on blank/garbage."""
    s = (s or "").strip()
    if not s:
        return 0
    try:
        return int(s)
    except ValueError:
        return 0


def parse_runinfo(text: str) -> list[dict]:
    """Extract per-run records from an SRA runinfo CSV.

    Each record is a dict with keys:
      name       – ScientificName
      taxid      – TaxID (int) or None
      run        – Run accession (SRR/ERR/DRR…) or ""  (used to de-dup reads)
      experiment – Experiment accession (SRX…) or ""   (counted as 1 experiment)
      reads      – spots + spots_with_mates (paired reads counted as two)
    """
    out: list[dict] = []
    if not text or not text.strip():
        return out
    try:
        reader = csv.DictReader(io.StringIO(text))
        for row in reader:
            name = (row.get("ScientificName") or "").strip()
            tx   = (row.get("TaxID") or "").strip()
            if not name and not tx:
                continue
            reads = _to_int(row.get("spots")) + _to_int(row.get("spots_with_mates"))
            out.append({
                "name"      : name,
                "taxid"     : int(tx) if tx.isdigit() else None,
                "run"       : (row.get("Run") or "").strip(),
                "experiment": (row.get("Experiment") or "").strip(),
                "reads"     : reads,
            })
    except csv.Error:
        pass
    return out


def fetch_sra_organisms(session, term, delay, api_key, email, cap):
    """
    Run an SRA esearch (history) + paginated runinfo efetch.
    Returns (list[record], total_count), where each record is the per-run
    dict produced by parse_runinfo (name, taxid, run, experiment, reads).
    """
    esp = {
        "db": "sra", "term": term, "usehistory": "y",
        "retmode": "json", "retmax": "0",
        "tool": "protist_survey", "email": email,
    }
    if api_key:
        esp["api_key"] = api_key
    data = ncbi_get(session, f"{ENTREZ_BASE}/esearch.fcgi", esp, delay)
    res  = data.get("esearchresult", {}) if data else {}
    try:
        count = int(res.get("count", "0"))
    except (TypeError, ValueError):
        count = 0
    webenv, qk = res.get("webenv"), res.get("querykey")
    if not count or not webenv or not qk:
        return [], count

    n, batch, start, orgs = min(count, cap), 5000, 0, []
    while start < n:
        fp = {
            "db": "sra", "query_key": qk, "WebEnv": webenv,
            "rettype": "runinfo", "retmode": "text",
            "retstart": str(start), "retmax": str(min(batch, n - start)),
            "tool": "protist_survey", "email": email,
        }
        if api_key:
            fp["api_key"] = api_key
        chunk = parse_runinfo(ncbi_get_text(session, f"{ENTREZ_BASE}/efetch.fcgi", fp, delay))
        if not chunk:
            break
        orgs.extend(chunk)
        start += batch
    return orgs, count


def _absorb_runs(records, group, org_set, name_map, group_map,
                 exp_map, reads_map, seen_runs):
    """Fold SRA runinfo records into per-taxid experiment/read aggregates.

      org_set   : set[int]            – taxids seen for this read type
      name_map  : dict[int, str]      – taxid → first ScientificName
      group_map : dict[int, str]      – taxid → protist group
      exp_map   : dict[int, set[str]] – taxid → {Experiment accessions}
      reads_map : dict[int, int]      – taxid → total reads
      seen_runs : set[str]            – Run accessions already counted (read de-dup)
    """
    for rec in records:
        tx = rec["taxid"]
        if not tx:
            continue
        org_set.add(tx)
        name_map.setdefault(tx, rec["name"])
        group_map.setdefault(tx, group)

        run = rec["run"]
        if not run or run not in seen_runs:      # count reads once per run
            if run:
                seen_runs.add(run)
            reads_map[tx] = reads_map.get(tx, 0) + rec["reads"]

        exp = rec["experiment"] or run            # fall back to run if no SRX
        if exp:
            exp_map.setdefault(tx, set()).add(exp)


# ──────────────────────────────────────────────────────────────────────────────
#  NCBI Entrez - Taxonomy: resolve each taxid up to species rank
# ──────────────────────────────────────────────────────────────────────────────

def parse_taxonomy_xml(xml: str) -> dict[int, tuple[int, str]]:
    """
    Parse an Entrez taxonomy efetch XML payload.
    Returns {queried_taxid: (species_taxid, species_name)}.
    A taxon below species rank is mapped to its 'species' ancestor in LineageEx;
    a taxon at/above species rank maps to itself.
    """
    result: dict[int, tuple[int, str]] = {}
    if not xml:
        return result
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return result
    for tx in root.findall("Taxon"):
        tid_txt = tx.findtext("TaxId")
        if not tid_txt or not tid_txt.isdigit():
            continue
        tid  = int(tid_txt)
        rank = (tx.findtext("Rank") or "").lower()
        name = tx.findtext("ScientificName") or str(tid)
        sp_tid, sp_name = tid, name
        if rank != "species":
            lex = tx.find("LineageEx")
            if lex is not None:
                for anc in lex.findall("Taxon"):
                    if (anc.findtext("Rank") or "").lower() == "species":
                        a_tid = anc.findtext("TaxId")
                        if a_tid and a_tid.isdigit():
                            sp_tid  = int(a_tid)
                            sp_name = anc.findtext("ScientificName") or sp_name
                        break
        result[tid] = (sp_tid, sp_name)
    return result


def resolve_species(session, taxids, delay, api_key, email) -> dict[int, tuple[int, str]]:
    ids = sorted({int(t) for t in taxids if t})
    out: dict[int, tuple[int, str]] = {}
    BATCH = 180
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        params = {
            "db": "taxonomy", "id": ",".join(map(str, chunk)),
            "retmode": "xml", "tool": "protist_survey", "email": email,
        }
        if api_key:
            params["api_key"] = api_key
        parsed = parse_taxonomy_xml(
            ncbi_get_text(session, f"{ENTREZ_BASE}/efetch.fcgi", params, delay))
        out.update(parsed)
        for t in chunk:                      # fallback for anything unresolved
            out.setdefault(t, (t, str(t)))
        print(f"\r  resolved {len(out)}/{len(ids)} taxids …", end="", flush=True)
    print()
    return out


# ──────────────────────────────────────────────────────────────────────────────
#  Labels & progress
# ──────────────────────────────────────────────────────────────────────────────

def rnaseq_label(has_short: bool, has_long: bool) -> str:
    if has_short and has_long:
        return "short-read, long-read"
    if has_short:
        return "short-read only"
    if has_long:
        return "long-read only"
    return "None"


def print_progress(label, idx, total, n, what="assemblies"):
    width  = 26
    filled = int(width * idx / total)
    bar    = "█" * filled + "░" * (width - filled)
    print(f"\r[{bar}] {100*idx/total:5.1f}%  {label:<22} {n:>5} {what}",
          end="", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
#  HTML report
# ──────────────────────────────────────────────────────────────────────────────

_RNASEQ_BADGE = {
    "short-read, long-read": ("#cce5ff", "#004085"),
    "short-read only"       : ("#d6eaf8", "#1a5276"),
    "long-read only"        : ("#e8daef", "#512e5f"),
    "None"                  : ("#f2f3f4", "#616a6b"),
    "N/A (--no-sra)"        : ("#f2f3f4", "#616a6b"),
}


def _rnaseq_badge(text):
    bg, fg = _RNASEQ_BADGE.get(text, ("#fff3cd", "#856404"))
    return (f'<span style="background:{bg};color:{fg};padding:2px 9px;'
            f'border-radius:4px;font-size:.82em;white-space:nowrap">{text}</span>')


def _ann_badge(acc):
    if acc.startswith("GCF"):
        bg, fg = "#d4edda", "#155724"
    elif acc.startswith("GCA"):
        bg, fg = "#fff3cd", "#856404"
    else:
        bg, fg = "#f8d7da", "#721c24"
    return (f'<span style="background:{bg};color:{fg};padding:2px 9px;'
            f'border-radius:4px;font-size:.82em;font-family:monospace;'
            f'white-space:nowrap">{acc}</span>')


def _num_cell(v):
    """Right-aligned numeric cell with a raw `data-sort` value for sorting.
    Shows a muted dash when zero/blank/NA (sorts as 0)."""
    try:
        n = int(v)
    except (TypeError, ValueError):
        n = 0
    if n == 0:
        return ('<td data-sort="0" style="text-align:right;color:#aaa">–</td>')
    return (f'<td data-sort="{n}" style="text-align:right;'
            f'font-variant-numeric:tabular-nums">{n:,}</td>')


def _gc_cell(v):
    """Right-aligned GC% cell with raw `data-sort`; dash when unavailable
    (missing GC sorts below all real values via data-sort=-1)."""
    try:
        g = float(v)
    except (TypeError, ValueError):
        return '<td data-sort="-1" style="text-align:right;color:#aaa">–</td>'
    return (f'<td data-sort="{g:.4f}" style="text-align:right;'
            f'font-variant-numeric:tabular-nums">{g:.1f}%</td>')


def save_html(df, path):
    rows_html = "\n".join(
        f"<tr><td>{r['SuperGroup']}</td><td>{r['Group']}</td>"
        f"<td><em>{r['Species']}</em></td>"
        f"<td style='font-family:monospace'>{r['TaxID']}</td>"
        f"<td>{_ann_badge(r['Assembly Accession'])}</td>"
        f"{_gc_cell(r['GC%'])}"
        f"<td>{_ann_badge(r['Annotation Accession'])}</td>"
        f"<td>{_rnaseq_badge(r['RNAseq Data'])}</td>"
        f"{_num_cell(r['Long-read Experiments'])}"
        f"{_num_cell(r['Long-read Reads'])}"
        f"{_num_cell(r['Short-read Experiments'])}"
        f"{_num_cell(r['Short-read Reads'])}</tr>"
        for _, r in df.iterrows()
    )
    html = dedent(f"""\
    <!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
    <title>Protist Genomic Overview</title>
    <style>
      *{{box-sizing:border-box;margin:0;padding:0}}
      body{{font-family:"Segoe UI",Arial,sans-serif;background:#f4f6f9;color:#333;padding:2em}}
      h1{{font-size:1.5em;color:#2c3e50;margin-bottom:.3em}}
      .sub{{color:#777;font-size:.88em;margin-bottom:1em}}
      .legend{{display:flex;gap:.6em;margin-bottom:1em;flex-wrap:wrap;font-size:.8em}}
      .legend span{{padding:2px 9px;border-radius:4px}}
      #search{{padding:7px 12px;border:1px solid #ccc;border-radius:5px;font-size:.93em;
              width:340px;margin-bottom:.6em;outline:none}}
      #search:focus{{border-color:#3498db}}
      table{{border-collapse:collapse;width:100%;background:#fff;border-radius:8px;
            overflow:hidden;box-shadow:0 2px 8px rgba(0,0,0,.08);font-size:.86em}}
      thead th{{background:#2c3e50;color:#fff;padding:11px 14px;text-align:left;
               position:sticky;top:0;z-index:1}}
      thead th.sortable{{cursor:pointer;user-select:none;white-space:nowrap}}
      thead th.sortable:hover{{background:#3a5169}}
      thead th.num{{text-align:right}}
      th .arrow{{display:inline-block;width:1em;color:#9fb3c8;font-size:.85em}}
      th.sorted .arrow{{color:#fff}}
      tbody tr:nth-child(even){{background:#f7f9fc}}
      tbody tr:hover{{background:#eaf2ff}}
      td{{padding:7px 14px;border-bottom:1px solid #e8e8e8}}
      td:nth-child(3){{font-style:italic}}
      #count{{color:#555;font-size:.85em;margin-bottom:.6em}}
    </style></head><body>
      <h1>Protist Genomic Overview</h1>
      <p class="sub">{len(df):,} species across {df['Group'].nunique()} groups
        ({df['SuperGroup'].nunique()} supergroups) -
        NCBI Datasets v2, SRA &amp; Taxonomy. A row with an Assembly Accession but
        “-” in Annotation Accession is a genome assembly without annotation; “-”
        in both columns (with blank GC%) marks an SRA-only species (no genome).</p>
      <div class="legend"><b>Accession (assembly &amp; annotation):</b>
        <span style="background:#d4edda;color:#155724">GCF_ RefSeq</span>
        <span style="background:#fff3cd;color:#856404">GCA_ GenBank</span>
        <span style="background:#f8d7da;color:#721c24">- none</span></div>
      <p class="sub" style="margin-bottom:.6em">Tip: click a column header to sort.
        SuperGroup, Protist Group and Species sort alphabetically (A→Z, then Z→A);
        GC% and the read columns sort numerically (max→min, then min→max).</p>
      <input id="search" type="text" placeholder="Filter by group, species, accession …"
             oninput="filterTable()">
      <div id="count">{len(df):,} entries shown</div>
      <table id="tbl"><thead><tr>
        <th class="sortable" onclick="sortTxt(this.cellIndex)">SuperGroup<span class="arrow"></span></th>
        <th class="sortable" onclick="sortTxt(this.cellIndex)">Protist Group<span class="arrow"></span></th>
        <th class="sortable" onclick="sortTxt(this.cellIndex)">Species<span class="arrow"></span></th>
        <th>TaxID</th>
        <th>Assembly Accession</th>
        <th class="sortable num" onclick="sortNum(this.cellIndex)">GC%<span class="arrow"></span></th>
        <th>Annotation Accession</th>
        <th>RNAseq Data</th>
        <th class="sortable num" onclick="sortNum(this.cellIndex)">Long-read Exp.<span class="arrow"></span></th>
        <th class="sortable num" onclick="sortNum(this.cellIndex)">Long-read Reads<span class="arrow"></span></th>
        <th class="sortable num" onclick="sortNum(this.cellIndex)">Short-read Exp.<span class="arrow"></span></th>
        <th class="sortable num" onclick="sortNum(this.cellIndex)">Short-read Reads<span class="arrow"></span></th></tr></thead>
      <tbody>
    {rows_html}
      </tbody></table>
      <script>
        function filterTable(){{
          const q=document.getElementById('search').value.toLowerCase();
          const rows=document.querySelectorAll('#tbl tbody tr');let v=0;
          rows.forEach(r=>{{const s=r.textContent.toLowerCase().includes(q);
            r.style.display=s?'':'none';if(s)v++;}});
          document.getElementById('count').textContent=v+' entries shown';}}

        let sortCol=null, sortDir=-1;          // num: -1 max→min / +1 min→max; txt: +1 A→Z / -1 Z→A
        function sortBy(col, type){{
          if(col===sortCol){{sortDir*=-1;}}                       // same column → flip direction
          else{{sortCol=col;sortDir=(type==='txt')?1:-1;}}        // txt starts A→Z, num starts max→min
          const tbody=document.querySelector('#tbl tbody');
          const rows=Array.from(tbody.querySelectorAll('tr'));
          rows.sort((a,b)=>{{
            if(type==='txt'){{
              const av=a.children[col].textContent.trim().toLowerCase();
              const bv=b.children[col].textContent.trim().toLowerCase();
              return av<bv?-sortDir:av>bv?sortDir:0;
            }}
            const av=parseFloat(a.children[col].dataset.sort)||0;
            const bv=parseFloat(b.children[col].dataset.sort)||0;
            return (av-bv)*sortDir;}});
          rows.forEach(r=>tbody.appendChild(r));
          document.querySelectorAll('#tbl thead th').forEach((th,i)=>{{
            const ar=th.querySelector('.arrow');
            if(!ar)return;
            th.classList.toggle('sorted', i===col);
            ar.textContent = (i===col) ? (sortDir<0?'▼':'▲') : '';}});
        }}
        function sortNum(col){{sortBy(col,'num');}}
        function sortTxt(col){{sortBy(col,'txt');}}
      </script>
    </body></html>
    """)
    path.write_text(html, encoding="utf-8")


# ──────────────────────────────────────────────────────────────────────────────
#  Core pipeline
# ──────────────────────────────────────────────────────────────────────────────

_JUNK = ("uncultured", "environmental", "metagenome", "synthetic ", "mixed ")


def _finish(rows):
    if not rows:
        print("\n[ERROR] No data retrieved - check connectivity to "
              "api.ncbi.nlm.nih.gov / eutils.ncbi.nlm.nih.gov.", file=sys.stderr)
        sys.exit(1)
    df = pd.DataFrame(rows, columns=[
        "Group", "Species", "TaxID",
        "Assembly Accession", "GC%", "Annotation Accession", "RNAseq Data",
        "Long-read Experiments", "Long-read Reads",
        "Short-read Experiments", "Short-read Reads"])
    df.insert(0, "SuperGroup", df["Group"].map(SUPERGROUPS).fillna("Other"))
    return df.sort_values(["SuperGroup", "Group", "Species"]).reset_index(drop=True)


def build_table(api_key, limit, skip_sra, email, sra_run_cap, assemblies_only=True):
    session, delay = make_session(api_key, email)
    total = len(PROTIST_GROUPS)

    # ── Pass 1 - assemblies (dedup by name, prefer GCF over GCA) ──────────────
    print(f"\nPass 1/4 - genome assemblies ({total} groups) …\n")
    species_best: dict[str, dict] = {}
    for idx, (group, taxid) in enumerate(PROTIST_GROUPS.items(), 1):
        reports = fetch_assemblies(session, taxid, limit, delay)
        print_progress(group, idx, total, len(reports))
        for rpt in reports:
            row = parse_report(rpt, group)
            if row is None:
                continue
            nm, acc = row["Species"], row["_accession"]
            if nm not in species_best:
                species_best[nm] = row
            elif acc.startswith("GCF") and not species_best[nm]["_accession"].startswith("GCF"):
                species_best[nm] = row
    print()
    assembly_rows = list(species_best.values())

    # ── --no-sra short-circuit ────────────────────────────────────────────────
    if skip_sra:
        rows = []
        for r in assembly_rows:
            r["RNAseq Data"] = "N/A (--no-sra)"
            r["Long-read Experiments"]  = ""
            r["Long-read Reads"]        = ""
            r["Short-read Experiments"] = ""
            r["Short-read Reads"]       = ""
            r.pop("_taxid", None); r.pop("_accession", None)
            rows.append(r)
        return _finish(rows), {"sra_only": 0}

    # ── Pass 2 - SRA-first discovery (long-read full, short-read capped) ──────
    print(f"\nPass 2/4 - SRA RNAseq discovery (long-read priority) …\n")
    sra_long: set[int]  = set()
    sra_short: set[int] = set()
    taxid_name: dict[int, str]  = {}
    taxid_group: dict[int, str] = {}
    # per-taxid experiment / read aggregates
    long_exp:    dict[int, set[str]] = {}
    long_reads:  dict[int, int]      = {}
    short_exp:   dict[int, set[str]] = {}
    short_reads: dict[int, int]      = {}
    long_seen_runs:  set[str] = set()      # de-dup reads by run accession
    short_seen_runs: set[str] = set()
    for idx, (group, taxid) in enumerate(PROTIST_GROUPS.items(), 1):
        long_term  = (f"txid{taxid}[Organism:exp] AND "
                      f"(RNA-Seq[Strategy] OR FL-cDNA[Strategy]) AND {_LONG_PLAT}")
        short_term = (f"txid{taxid}[Organism:exp] AND "
                      f"RNA-Seq[Strategy] AND {_SHORT_PLAT}")
        long_recs, _  = fetch_sra_organisms(session, long_term,  delay, api_key, email, cap=200_000)
        short_recs, _ = fetch_sra_organisms(session, short_term, delay, api_key, email, cap=sra_run_cap)
        _absorb_runs(long_recs,  group, sra_long,  taxid_name, taxid_group,
                     long_exp,  long_reads,  long_seen_runs)
        _absorb_runs(short_recs, group, sra_short, taxid_name, taxid_group,
                     short_exp, short_reads, short_seen_runs)
        print_progress(group, idx, total,
                       len({r["taxid"] for r in long_recs if r["taxid"]}),
                       what="long-read taxa")
    print()

    # ── Pass 3 - resolve every taxid to species rank ──────────────────────────
    print("\nPass 3/4 - resolving taxids to species level …")
    all_taxids = {r["_taxid"] for r in assembly_rows if r["_taxid"]}
    all_taxids |= sra_long | sra_short
    resolved = resolve_species(session, all_taxids, delay, api_key, email)

    # Both long-read and short-read presence stay at the EXACT taxid (sra_long /
    # sra_short, per strain) for assembly rows. species_long / species_short are
    # kept only for SRA-only rows (no assembly => no strain to attribute to).
    species_short = {resolved.get(t, (t, ""))[0] for t in sra_short}
    species_long  = {resolved.get(t, (t, ""))[0] for t in sra_long}

    # experiment/read aggregates rolled up to species rank (SRA-only rows only)
    sp_short_exp:   dict[int, set[str]] = {}
    sp_short_reads: dict[int, int]      = {}
    sp_long_exp:    dict[int, set[str]] = {}
    sp_long_reads:  dict[int, int]      = {}

    def _rollup(exp_map, reads_map, sp_exp, sp_reads):
        for tx, exps in exp_map.items():
            sp = resolved.get(tx, (tx, ""))[0]
            sp_exp.setdefault(sp, set()).update(exps)
        for tx, n in reads_map.items():
            sp = resolved.get(tx, (tx, ""))[0]
            sp_reads[sp] = sp_reads.get(sp, 0) + n

    _rollup(short_exp, short_reads, sp_short_exp, sp_short_reads)
    _rollup(long_exp,  long_reads,  sp_long_exp,  sp_long_reads)

    def _assembly_counts(exact_tid):
        """(long_exp, long_reads, short_exp, short_reads) for an assembly row,
        all keyed on the EXACT taxid (this strain only)."""
        return (len(long_exp.get(exact_tid, ())),   long_reads.get(exact_tid, 0),
                len(short_exp.get(exact_tid, ())),   short_reads.get(exact_tid, 0))

    def _species_counts(sp_tid):
        """(long_exp, long_reads, short_exp, short_reads) rolled up to species,
        used for SRA-only rows (no assembly => no strain to attribute)."""
        return (len(sp_long_exp.get(sp_tid, ())),   sp_long_reads.get(sp_tid, 0),
                len(sp_short_exp.get(sp_tid, ())),   sp_short_reads.get(sp_tid, 0))

    # ── Pass 4 - assemble final rows ──────────────────────────────────────────
    print("\nPass 4/4 - building table …")
    rows, covered = [], set()
    for r in assembly_rows:
        exact_tid = r["_taxid"]
        sp_tid, _ = resolved.get(exact_tid, (exact_tid, r["Species"]))
        covered.add(sp_tid)
        has_long  = exact_tid in sra_long          # strain-level: this taxid only
        has_short = exact_tid in sra_short         # strain-level: this taxid only
        r["RNAseq Data"] = rnaseq_label(has_short, has_long)
        le, lr, se, sr = _assembly_counts(exact_tid)
        r["Long-read Experiments"]  = le
        r["Long-read Reads"]        = lr
        r["Short-read Experiments"] = se
        r["Short-read Reads"]       = sr
        r.pop("_taxid", None); r.pop("_accession", None)
        rows.append(r)

    # SRA-only organisms (RNAseq present, no genome assembly)
    sra_only: dict[int, tuple[str, str]] = {}
    for tx in (sra_long | sra_short):
        sp_tid, sp_name = resolved.get(tx, (tx, taxid_name.get(tx, str(tx))))
        if sp_tid in covered:
            continue
        sra_only.setdefault(sp_tid, (sp_name, taxid_group.get(tx, "(SRA only)")))

    n_sra_only = 0
    for sp_tid, (name, group) in sra_only.items():
        if any(j in name.lower() for j in _JUNK):
            continue
        hs, hl = sp_tid in species_short, sp_tid in species_long
        if not (hs or hl):
            continue
        n_sra_only += 1
        if assemblies_only:                 # keep genome-backed species only
            continue
        le, lr, se, sr = _species_counts(sp_tid)
        rows.append({
            "Group": group, "Species": name, "TaxID": sp_tid,
            "Assembly Accession": "-", "GC%": "", "Annotation Accession": "-",
            "RNAseq Data": rnaseq_label(hs, hl),
            "Long-read Experiments": le,  "Long-read Reads": lr,
            "Short-read Experiments": se, "Short-read Reads": sr,
        })

    if assemblies_only:
        return _finish(rows), {"sra_only": 0, "sra_only_suppressed": n_sra_only}
    return _finish(rows), {"sra_only": n_sra_only}


# ──────────────────────────────────────────────────────────────────────────────
#  CLI
# ──────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        prog="protist_survey.py",
        description="Build a protist genomic-data table from NCBI.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("--output", "-o", default="protists.tsv", metavar="FILE",
                   help="Output .tsv or .csv (default: protists.tsv)")
    p.add_argument("--limit", "-l", type=int, default=None, metavar="N",
                   help="Cap assemblies per group (default: all)")
    p.add_argument("--api-key", "-k", default=None, metavar="KEY",
                   help="NCBI API key → 10 req/s (free at ncbi.nlm.nih.gov/account/)")
    p.add_argument("--email", "-e", default="anonymous@example.com", metavar="EMAIL",
                   help="E-mail for NCBI Entrez (recommended)")
    p.add_argument("--sra-run-cap", type=int, default=900_000, metavar="N",
                   help="Max short-read runinfo rows fetched per group (default 30000). "
                        "long-read is never capped.")
    p.add_argument("--no-sra", action="store_true",
                   help="Skip SRA + taxonomy passes (assembly-only, fast)")
    p.add_argument("--assemblies-only", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Only list species that have a genome assembly; drop "
                        "SRA-only species (no genome). Enabled by DEFAULT. Pass "
                        "--no-assemblies-only to also include SRA-only species. "
                        "RNAseq columns are filled for genome-backed species either way.")
    p.add_argument("--html", action="store_true",
                   help="Also write a searchable HTML report")
    args = p.parse_args()

    rate = "10 req/s" if args.api_key else "3 req/s  (full run can take ~1-2 h)"
    print("=" * 66)
    print("  NCBI Protist Genomic Data Survey  (SRA-first + species resolution)")
    print("=" * 66)
    print(f"  Groups           : {len(PROTIST_GROUPS)}")
    print(f"  Assembly limit   : {args.limit or 'all'} per group")
    print(f"  SRA RNAseq       : {'disabled (--no-sra)' if args.no_sra else 'enabled (long-read full, short-read capped)'}")
    print(f"  short-read cap   : {args.sra_run_cap:,} runs/group")
    print(f"  Scope            : {'genome-backed species only (default)' if args.assemblies_only else 'genomes + SRA-only species (--no-assemblies-only)'}")
    print(f"  Rate             : {rate}")
    print("=" * 66)

    df, info = build_table(args.api_key, args.limit, args.no_sra,
                           args.email, args.sra_run_cap,
                           assemblies_only=args.assemblies_only)

    out = Path(args.output)
    sep = "\t" if out.suffix in (".tsv", ".txt") else ","

    # Output cleanup: blank out "-"/"None", and turn "&" into a comma.
    df_out = df.replace({"-": "", "None": ""})
    df_out["RNAseq Data"] = df_out["RNAseq Data"].str.replace(" & ", ", ", regex=False)

    # File headers carry no spaces (machine-friendly); HTML keeps readable names.
    df_out.columns = [c.lower().replace(" ", "_") for c in df_out.columns]

    df_out.to_csv(out, sep=sep, index=False)
    print(f"\n✓  {len(df):,} rows saved → {out.resolve()}")

    if args.html:
        hp = out.with_suffix(".html")
        save_html(df, hp)
        print(f"✓  HTML report   → {hp.resolve()}")

    # ── Summary ──────────────────────────────────────────────────────────────
    has_asm = int((df["Assembly Accession"] != "-").sum())
    gcf  = int(df["Annotation Accession"].str.startswith("GCF").sum())
    gca  = int(df["Annotation Accession"].str.startswith("GCA").sum())
    annotated  = gcf + gca
    asm_no_ann = has_asm - annotated
    sra_only_n = info.get("sra_only", 0)
    suppressed = info.get("sra_only_suppressed", 0)
    gc_vals = pd.to_numeric(df["GC%"], errors="coerce").dropna()
    longr  = df["RNAseq Data"].str.contains("Long",  na=False).sum()
    shortr = df["RNAseq Data"].str.contains("Short", na=False).sum()
    norna  = (df["RNAseq Data"] == "None").sum()

    print(f"\n── Summary {'─' * 53}")
    print(f"  Total species                 : {len(df):>6,}")
    print(f"  With genome assembly          : {has_asm:>6,}")
    print(f"    • annotated (GCF+GCA)       : {annotated:>6,}  (GCF {gcf:,} / GCA {gca:,})")
    print(f"    • assembly, no annotation   : {asm_no_ann:>6,}")
    if suppressed:
        print(f"  SRA-only species (suppressed) : {suppressed:>6,}  (hidden by --assemblies-only)")
    else:
        print(f"  SRA-only (no genome)          : {sra_only_n:>6,}")
    if len(gc_vals):
        print(f"  GC% reported                  : {len(gc_vals):>6,}  "
              f"(mean {gc_vals.mean():.1f}%, range {gc_vals.min():.1f}–{gc_vals.max():.1f}%)")
    if not args.no_sra:
        print(f"  With long-read  RNAseq        : {longr:>6,}")
        print(f"  With short-read RNAseq        : {shortr:>6,}")
        print(f"  No RNAseq data                : {norna:>6,}")
        le = int(pd.to_numeric(df["Long-read Experiments"],  errors="coerce").sum())
        lr = int(pd.to_numeric(df["Long-read Reads"],        errors="coerce").sum())
        se = int(pd.to_numeric(df["Short-read Experiments"], errors="coerce").sum())
        sr = int(pd.to_numeric(df["Short-read Reads"],       errors="coerce").sum())
        print(f"  Long-read  experiments (total): {le:>6,}")
        print(f"  Long-read  reads (total)      : {lr:>14,}")
        print(f"  Short-read experiments (total): {se:>6,}")
        print(f"  Short-read reads (total)      : {sr:>14,}")

    print(f"\n── Species per group {'─' * 43}")
    print(df.groupby("Group").size().rename("# species")
            .sort_values(ascending=False).to_string())

    print(f"\n── First 25 rows {'─' * 47}")
    pd.set_option("display.max_colwidth", 38)
    pd.set_option("display.width", 140)
    print(df.head(25).to_string(index=False))
    print("─" * 66)
    print(f"Full table → {out}")


if __name__ == "__main__":
    main()
