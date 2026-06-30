#!/usr/bin/env python3
"""Augment the protist table with per-species results from isoquant, LyRic and geneid.

For every species in the protist table that can be matched against a tool's
output, the matching result columns are appended (prefixed by the tool name).
Species with no result for a given tool keep NA in that tool's columns.

The three pipelines write incompatible summary layouts, so each source is read
with its own path and species-key format; everything is then joined back onto
the protist table.

Each module's summary table is always regenerated from its per-species outputs
before being read, so the final table is never stale (use --no-refresh to build
only missing tables, or --no-generate to read existing tables only):
  reference     busco_references/scripts/make_summary_tables.py -s <summary_dir>
  isoquant      isoquant_annotator/scripts/make_summary_tables.py <summary_dir>
  lyric         LyRic_annotator/scripts/make_summary_tables.py [--merged|--joint]
  geneid        geneid-training/scripts/counting_machine.py -r <results_dir>
These generators read each module's per-species evaluation outputs; if those
upstream files are absent too, nothing is written and the source is skipped.

With --html, a searchable HTML report (protist_table style) is written next to
the output TSV.

Sources (relative to the projects root, all optional - missing ones are skipped)
  reference      busco_references/summary/general_summary.tsv
                 base genome + annotation metrics (gene/transcript counts, BUSCO);
                 key: <Genus_species[_strain]> + a separate taxonID column
  isoquant       isoquant_annotator/summary/general_summary.tsv
                 key: <Genus_species[_strain]>_<taxid>
  lyric          LyRic_annotator/summary/general_summary.tsv
                 key: <Genus_species[_strain]>  (taxid already stripped)
  lyric_merged   LyRic_annotator/summary/merge_general_summary.tsv
                 key: <Genus_species[_strain]>
  lyric (joint)  LyRic_annotator/summary/joint_summary.tsv
                 unified table (regular+merged, columns suffixed _regular/_merged);
                 when this file exists it REPLACES the two separate lyric tables
  geneid_<flav>  geneid-training/results/summary/<flav>/general_summary.tsv
                 flav in {regular, merged}
                 key: <Genus_species[_strain]>_<taxid>
                 counts + pivoted BUSCO (own/git lineage and eukaryote completeness)

Matching strategy (per source row, first hit wins)
  1. exact TaxID  (only when the source key carries a trailing _<digits>)
  2. exact normalized species name
  3. normalized name as a strain-prefix of the source key
     (handles 'Cyanidioschyzon_merolae_strain_10D' vs 'Cyanidioschyzon merolae')

Sections can be deactivated when a tool's data is unavailable; a deactivated
section adds no columns and never runs its generator:
  --skip geneid            drop the geneid section
  --skip geneid lyric      drop several at once
  --only isoquant          keep only one section (overrides --skip)
sections: reference, isoquant, lyric, geneid.

Usage (from the projects root or anywhere; defaults resolve off this file)
  python3 result_analytics/build_results_table.py
  python3 result_analytics/build_results_table.py --out my_table.tsv
  python3 result_analytics/build_results_table.py --isoquant /path/general_summary.tsv
  python3 result_analytics/build_results_table.py --skip geneid
"""
import argparse
import os
import re
import subprocess
import sys
from html import escape
from textwrap import dedent

import pandas as pd

# projects root = parent of this script's directory (result_analytics/)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ISO_DIR = os.path.join(ROOT, "isoquant_annotator")
LYRIC_DIR = os.path.join(ROOT, "lyric_annotator")
GENEID_DIR = os.path.join(ROOT, "geneid-training")
BUSCO_REF_DIR = os.path.join(ROOT, "busco_references")

# Deactivatable sections (--skip / --only). Each maps to one source group below.
SECTIONS = ["reference", "isoquant", "lyric", "geneid"]


def normalize_name(text):
    """Lowercase, drop dots, spaces/underscores -> single underscore.

    'Apicomplexa sp.' -> 'apicomplexa_sp' ; 'Genus__species' -> 'genus_species'.
    """
    text = str(text).strip().lower().replace(".", "")
    text = re.sub(r"[\s_]+", "_", text)
    return text.strip("_")


def split_stem(stem):
    """'<name>_<taxid>' -> (normalized_name, taxid_str); no trailing digits -> (norm, None)."""
    m = re.match(r"^(.+)_(\d+)$", str(stem))
    if m:
        return normalize_name(m.group(1)), m.group(2)
    return normalize_name(stem), None


def run_generator(gen, ran):
    """Run a module's summary-table generator once (deduped via the `ran` set).

    gen = {"cwd": ..., "cmd": [...]}. Generators read each module's per-species
    outputs and write the summary tables this script later consumes; if those
    upstream inputs are absent the generator writes nothing (or an empty table),
    which is handled by the caller re-checking for the file.
    """
    key = (gen["cwd"], tuple(gen["cmd"]))
    if key in ran:
        return
    ran.add(key)
    print(f"[gen]  running {' '.join(gen['cmd'])} (cwd={gen['cwd']})")
    try:
        res = subprocess.run(gen["cmd"], cwd=gen["cwd"], capture_output=True, text=True)
    except Exception as exc:  # noqa: BLE001 - report any launch failure and continue
        print(f"[gen]  error launching generator: {exc}")
        return
    if res.returncode != 0:
        tail = (res.stderr.strip().splitlines() or [""])[-1]
        print(f"[gen]  generator exited {res.returncode}: {tail}")
    else:
        tail = (res.stdout.strip().splitlines() or ["done"])[-1]
        print(f"[gen]  {tail}")


# Badge palettes copied verbatim from protist_table/protist_survey.py so the
# shared columns (Annotation Accession, RNAseq Data) render identically here.
_RNASEQ_BADGE = {
    "short-read, long-read": ("#cce5ff", "#004085"),
    "short-read only": ("#d6eaf8", "#1a5276"),
    "long-read only": ("#e8daef", "#512e5f"),
    "None": ("#f2f3f4", "#616a6b"),
    "N/A (--no-sra)": ("#f2f3f4", "#616a6b"),
}


def _rnaseq_badge(text):
    bg, fg = _RNASEQ_BADGE.get(text, ("#fff3cd", "#856404"))
    return (f'<span style="background:{bg};color:{fg};padding:2px 9px;'
            f'border-radius:4px;font-size:.82em;white-space:nowrap">{escape(text)}</span>')


def _ann_badge(acc):
    if acc.startswith("GCF"):
        bg, fg = "#d4edda", "#155724"
    elif acc.startswith("GCA"):
        bg, fg = "#fff3cd", "#856404"
    else:
        bg, fg = "#f8d7da", "#721c24"
    return (f'<span style="background:{bg};color:{fg};padding:2px 9px;'
            f'border-radius:4px;font-size:.82em;font-family:monospace;'
            f'white-space:nowrap">{escape(acc)}</span>')


def _gc_cell(v):
    """Right-aligned GC% cell with raw `data-sort`; dash (sorts below) when missing."""
    try:
        g = float(v)
    except (TypeError, ValueError):
        return '<td data-sort="-1" style="text-align:right;color:#aaa">–</td>'
    return (f'<td data-sort="{g:.4f}" style="text-align:right;'
            f'font-variant-numeric:tabular-nums">{g:.1f}%</td>')


def _num_cell(v):
    """Right-aligned numeric cell with raw `data-sort`. Integers get thousands
    separators; non-integers keep up to 2 decimals; blank/NA sorts as 0."""
    try:
        g = float(v)
    except (TypeError, ValueError):
        return '<td data-sort="-1" style="text-align:right;color:#aaa">–</td>'
    if g == 0:
        return '<td data-sort="0" style="text-align:right;color:#aaa">–</td>'
    # integral values emit an int data-sort (matches protists.html exactly);
    # non-integers (e.g. BUSCO %) keep 2 decimals.
    sort_val = int(g) if g == int(g) else g
    body = f"{int(g):,}" if g == int(g) else f"{g:,.2f}"
    return (f'<td data-sort="{sort_val}" style="text-align:right;'
            f'font-variant-numeric:tabular-nums">{body}</td>')


def _is_numeric_column(series):
    """True when every non-blank value in the column parses as a float."""
    seen = False
    for v in series:
        if pd.isna(v) or str(v).strip() in ("", "NA"):
            continue
        try:
            float(v)
        except (TypeError, ValueError):
            return False
        seen = True
    return seen


def save_html(df, path):
    """Write a self-contained HTML report matching protist_table's protists.html.

    Same color scheme (header #2c3e50, sortable-hover #3a5169, zebra #f7f9fc,
    row-hover #eaf2ff, arrow #9fb3c8) and the same clickable column sorting
    (numeric headers carry `sortNum`; first click max->min, click again min->max).

    The shared base columns reuse protist_survey.py's badges (Annotation
    Accession, RNAseq Data) and GC% cell. Every other purely-numeric column
    (read counts, gene/transcript counts, BUSCO scores, ...) is auto-detected and
    made sortable with right-aligned tabular numbers.
    """
    columns = list(df.columns)

    # base protist columns get the exact same treatment as protist_table/protists.html.
    # keyed by lowercased name so it works whether the table is Title- or lower-cased.
    SPECIAL = {"species", "taxid", "assembly_accession", "annotation_accession",
               "rnaseq_data", "gc%"}
    numeric = {c: _is_numeric_column(df[c]) for c in columns if str(c).lower() not in SPECIAL}

    def render(col, val):
        key = str(col).lower()
        blank = pd.isna(val) or str(val).strip() in ("", "NA")
        if key == "gc%":
            return _gc_cell(None if blank else val)
        # accession columns always badge; missing -> red "none" badge with "-",
        # matching protist_table/protists.html (which renders missing as "-").
        if key in ("assembly_accession", "annotation_accession"):
            return f"<td>{_ann_badge('-' if blank else str(val))}</td>"
        # rnaseq always badges; blank -> "None" badge (protists.html stores missing as "None").
        if key == "rnaseq_data":
            return f"<td>{_rnaseq_badge('None' if blank else str(val))}</td>"
        if key in SPECIAL:
            if blank:
                return '<td style="color:#bbb;text-align:center">–</td>'
            text = escape(str(val))
            if key == "species":
                return f"<td><em>{text}</em></td>"
            if key == "taxid":
                return f"<td style='font-family:monospace'>{text}</td>"
        if numeric.get(col):
            return _num_cell(None if blank else val)
        if blank:
            return '<td style="color:#bbb;text-align:center">–</td>'
        return f"<td>{escape(str(val))}</td>"

    def sortable(col):
        return str(col).lower() == "gc%" or numeric.get(col, False)

    header_cells = []
    for i, col in enumerate(columns):
        label = escape(str(col))
        if sortable(col):
            header_cells.append(
                f'<th class="sortable num" onclick="sortNum({i})">'
                f'{label}<span class="arrow"></span></th>')
        else:
            header_cells.append(f"<th>{label}</th>")
    header = "".join(header_cells)

    rows_html = "\n".join(
        "<tr>" + "".join(render(c, r[c]) for c in columns) + "</tr>"
        for _, r in df.iterrows()
    )
    html = dedent(f"""\
    <!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
    <title>Protist Results Overview</title>
    <style>
      *{{box-sizing:border-box;margin:0;padding:0}}
      body{{font-family:"Segoe UI",Arial,sans-serif;background:#f4f6f9;color:#333;padding:2em}}
      h1{{font-size:1.5em;color:#2c3e50;margin-bottom:.3em}}
      .sub{{color:#777;font-size:.88em;margin-bottom:1em}}
      #search{{padding:7px 12px;border:1px solid #ccc;border-radius:5px;font-size:.93em;
              width:340px;margin-bottom:.6em;outline:none}}
      #search:focus{{border-color:#3498db}}
      .wrap{{overflow-x:auto;border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,.08)}}
      table{{border-collapse:collapse;width:100%;background:#fff;font-size:.84em}}
      thead th{{background:#2c3e50;color:#fff;padding:10px 13px;text-align:left;
               position:sticky;top:0;z-index:1;white-space:nowrap}}
      thead th.sortable{{cursor:pointer;user-select:none}}
      thead th.sortable:hover{{background:#3a5169}}
      thead th.num{{text-align:right}}
      th .arrow{{display:inline-block;width:1em;color:#9fb3c8;font-size:.85em}}
      th.sorted .arrow{{color:#fff}}
      tbody tr:nth-child(even){{background:#f7f9fc}}
      tbody tr:hover{{background:#eaf2ff}}
      td{{padding:6px 13px;border-bottom:1px solid #e8e8e8;white-space:nowrap}}
      #count{{color:#555;font-size:.85em;margin-bottom:.6em}}
    </style></head><body>
      <h1>Protist Results Overview</h1>
      <p class="sub">{len(df):,} species; base genome/annotation metrics plus
        isoquant, LyRic and geneid outputs merged per species.</p>
      <div class="legend" style="display:flex;gap:.6em;margin-bottom:1em;flex-wrap:wrap;font-size:.8em">
        <b>Accession (assembly &amp; annotation):</b>
        <span style="background:#d4edda;color:#155724;padding:2px 9px;border-radius:4px">GCF_ RefSeq</span>
        <span style="background:#fff3cd;color:#856404;padding:2px 9px;border-radius:4px">GCA_ GenBank</span>
        <span style="background:#f8d7da;color:#721c24;padding:2px 9px;border-radius:4px">- none</span></div>
      <p class="sub" style="margin-bottom:.6em">Tip: click a numeric column header
        to sort (first click max&rarr;min, click again for min&rarr;max).</p>
      <input id="search" type="text" placeholder="Filter by any column …" oninput="filterTable()">
      <div id="count">{len(df):,} entries shown</div>
      <div class="wrap"><table id="tbl"><thead><tr>{header}</tr></thead>
      <tbody>
    __ROWS__
      </tbody></table></div>
      <script>
        function filterTable(){{
          const q=document.getElementById('search').value.toLowerCase();
          const rows=document.querySelectorAll('#tbl tbody tr');let v=0;
          rows.forEach(r=>{{const s=r.textContent.toLowerCase().includes(q);
            r.style.display=s?'':'none';if(s)v++;}});
          document.getElementById('count').textContent=v+' entries shown';}}

        let sortCol=null, sortDir=-1;          // dir -1 = max->min, +1 = min->max
        function sortNum(col){{
          if(col===sortCol){{sortDir*=-1;}}     // same column -> flip direction
          else{{sortCol=col;sortDir=-1;}}        // new column -> start max->min
          const tbody=document.querySelector('#tbl tbody');
          const rows=Array.from(tbody.querySelectorAll('tr'));
          rows.sort((a,b)=>{{
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
      </script>
    </body></html>
    """)
    # insert rows after dedent so their unindented <tr> lines don't defeat it
    html = html.replace("__ROWS__", rows_html)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)


def find_col(columns, *candidates):
    """Resolve a column name case-insensitively (protists.tsv ships lowercased)."""
    lower = {str(c).lower(): c for c in columns}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    return None


def build_indexes(protists):
    """Return (taxid -> row-index, normalized-name -> row-index) for the protist table.

    First occurrence wins on collisions so the augmented row count never grows.
    The species/taxid columns are resolved case-insensitively so the table stays
    consistent whether the source uses 'Species'/'TaxID' or 'species'/'taxid'.
    """
    species_col = find_col(protists.columns, "species")
    taxid_col = find_col(protists.columns, "taxid")
    if species_col is None:
        raise KeyError("protist table has no 'species' column")
    taxid_index, name_index, name_by_idx = {}, {}, {}
    for idx, row in protists.iterrows():
        norm = normalize_name(row[species_col])
        name_by_idx[idx] = norm
        if taxid_col is not None:
            taxid = str(row[taxid_col]).strip()
            if taxid and taxid.lower() != "nan":
                taxid_index.setdefault(taxid, idx)
        name_index.setdefault(norm, idx)
    return taxid_index, name_index, name_by_idx


def names_compatible(src_name, prot_name):
    """True when a source species name is plausibly the same taxon as a protist row.

    Equal, or one is a strain/variant prefix of the other (handles
    'cyanidioschyzon_merolae_strain_10d' vs 'cyanidioschyzon_merolae'). Used to
    reject taxid hits where the source carries a wrong taxid for a different
    species (e.g. a mislabeled 'Euglena_gracilis' tagged with Chaetoceros' taxid).
    An empty name on either side is non-blocking (returns True).
    """
    if not src_name or not prot_name:
        return True
    if src_name == prot_name:
        return True
    return src_name.startswith(prot_name + "_") or prot_name.startswith(src_name + "_")


def match_row(source_key, taxid_index, name_index, name_keys, name_by_idx):
    """Resolve a source species key to a protist row index, or None.

    A taxid hit is only accepted when the source species name is compatible with
    the protist row's name; otherwise we fall through to name-based matching so a
    source row carrying the wrong taxid cannot land on an unrelated species.
    """
    name, taxid = split_stem(source_key)
    if taxid is not None and taxid in taxid_index:
        idx = taxid_index[taxid]
        if names_compatible(name, name_by_idx.get(idx, "")):
            return idx
        # taxid matched a different species -> ignore it, try the name instead
    if name in name_index:
        return name_index[name]
    # strain-prefix fallback: protist name is a prefix of the (longer) source key
    candidates = [pn for pn in name_keys if name == pn or name.startswith(pn + "_")]
    if candidates:
        return name_index[max(candidates, key=len)]
    return None


def attach_source(protists, taxid_index, name_index, name_keys, name_by_idx, label, path,
                  gen=None, generate=False, refresh=True, ran=None, taxid_col=None):
    """Read one summary TSV and append its columns (prefixed) onto the protist table.

    When `generate` is set the module's generator (`gen`) is run to (re)build the
    table: always when `refresh` is True, otherwise only when the table is
    missing. Returns the number of source rows that matched a protist species.

    `taxid_col` names a dedicated taxid column in the source (e.g. the
    busco_references summary, whose species key carries no trailing taxid). When
    given it is matched first against the protist TaxID and is excluded from the
    appended value columns.
    """
    if generate and gen is not None and (refresh or not os.path.isfile(path)):
        reason = "refreshing" if os.path.isfile(path) else "table missing, generating"
        print(f"[gen]  {label}: {reason} ...")
        run_generator(gen, ran if ran is not None else set())
    if not os.path.isfile(path):
        print(f"[skip] {label}: not found at {path}")
        return 0

    df = pd.read_csv(path, sep="\t", dtype=str)
    if "species" not in df.columns:
        print(f"[skip] {label}: no 'species' column in {path}")
        return 0

    drop_cols = {"species", taxid_col} if taxid_col else {"species"}
    value_cols = [c for c in df.columns if c not in drop_cols]
    for col in value_cols:
        protists[f"{label}_{col}"] = pd.NA

    matched, unmatched, conflicts = 0, [], 0
    for _, srow in df.iterrows():
        idx = None
        # dedicated taxid column (e.g. busco_references): accept only when the
        # source name is compatible with the taxid-matched protist row, so a row
        # carrying a wrong taxid does not clobber an unrelated species.
        if taxid_col and taxid_col in df.columns:
            tx = str(srow[taxid_col]).strip()
            if tx and tx.lower() != "nan":
                cand = taxid_index.get(tx)
                if cand is not None:
                    if names_compatible(normalize_name(srow["species"]), name_by_idx.get(cand, "")):
                        idx = cand
                    else:
                        conflicts += 1
                        print(f"[warn] {label}: '{srow['species']}' has taxid {tx} which "
                              f"belongs to '{name_by_idx.get(cand)}'; ignoring taxid, matching by name")
        if idx is None:
            idx = match_row(srow["species"], taxid_index, name_index, name_keys, name_by_idx)
        if idx is None:
            unmatched.append(srow["species"])
            continue
        # do not let an empty/NA source value overwrite a value already filled by
        # an earlier (e.g. correctly-matched) row for the same protist species.
        for col in value_cols:
            val = srow[col]
            if pd.isna(val) or str(val).strip() in ("", "NA"):
                continue
            protists.at[idx, f"{label}_{col}"] = val
        matched += 1

    note = f" ({len(unmatched)} unmatched: {', '.join(unmatched[:5])}{'...' if len(unmatched) > 5 else ''})" if unmatched else ""
    warn = f" [{conflicts} taxid conflict(s) resolved by name]" if conflicts else ""
    print(f"[ok]   {label}: {matched}/{len(df)} rows matched -> +{len(value_cols)} cols{note}{warn}")
    return matched


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protists", default=os.path.join(ROOT, "protist_table", "protists.tsv"))
    parser.add_argument("--reference", default=os.path.join(ROOT, "busco_references", "summary", "general_summary.tsv"),
                        help="base genome/annotation metrics (gene/transcript counts + BUSCO)")
    parser.add_argument("--isoquant", default=os.path.join(ROOT, "isoquant_annotator", "summary", "general_summary.tsv"))
    parser.add_argument("--lyric", default=os.path.join(ROOT, "LyRic_annotator", "summary", "general_summary.tsv"))
    parser.add_argument("--lyric-merged", default=os.path.join(ROOT, "LyRic_annotator", "summary", "merge_general_summary.tsv"))
    parser.add_argument("--lyric-joint", default=os.path.join(ROOT, "LyRic_annotator", "summary", "joint_summary.tsv"),
                        help="unified LyRic table; when present it replaces --lyric and --lyric-merged")
    parser.add_argument("--geneid-results", default=os.path.join(ROOT, "geneid-training", "results"),
                        help="geneid results dir; reads summary/<flavour>/counts_summary.tsv")
    parser.add_argument("--geneid-flavours", nargs="+", default=["regular", "merged"])
    parser.add_argument("--lyric-mode", choices=["separate", "joint"], default="separate",
                        help="which LyRic table(s) to use/generate when no joint table exists yet")
    parser.add_argument("--skip", nargs="+", default=[], choices=SECTIONS, metavar="SECTION",
                        help=f"deactivate whole sections (no columns added, no generators run); "
                             f"choices: {', '.join(SECTIONS)}. e.g. --skip geneid")
    parser.add_argument("--only", nargs="+", default=None, choices=SECTIONS, metavar="SECTION",
                        help="keep ONLY these sections (overrides --skip)")
    parser.add_argument("--no-generate", dest="generate", action="store_false",
                        help="do not run module generators at all (only read existing tables)")
    parser.add_argument("--no-refresh", dest="refresh", action="store_false",
                        help="generate only when a table is missing instead of always refreshing it")
    parser.add_argument("--html", action="store_true",
                        help="also write a searchable HTML report next to --out")
    parser.add_argument("--out", default=os.path.join(ROOT, "result_analytics", "protists_annotated.tsv"))
    args = parser.parse_args()

    if not os.path.isfile(args.protists):
        parser.error(f"protist table not found: {args.protists}")

    # resolve which sections are active. --only wins over --skip when both given.
    if args.only is not None:
        active = [s for s in SECTIONS if s in set(args.only)]
    else:
        active = [s for s in SECTIONS if s not in set(args.skip)]
    if not active:
        parser.error("no sections active (everything was skipped); nothing to do")
    skipped = [s for s in SECTIONS if s not in active]
    print(f"Sections active: {', '.join(active)}"
          + (f" | deactivated: {', '.join(skipped)}" if skipped else ""))

    protists = pd.read_csv(args.protists, sep="\t", dtype=str)
    print(f"Loaded {len(protists)} protist rows from {args.protists}")

    taxid_index, name_index, name_by_idx = build_indexes(protists)
    name_keys = list(name_index.keys())

    # generators: how to (re)build each module's summary tables when missing.
    # absolute script paths + cwd so behaviour is independent of where we run from.
    py = sys.executable
    iso_script = os.path.join(ISO_DIR, "scripts", "make_summary_tables.py")
    lyric_script = os.path.join(LYRIC_DIR, "scripts", "make_summary_tables.py")
    geneid_script = os.path.join(GENEID_DIR, "scripts", "counting_machine.py")
    ref_script = os.path.join(BUSCO_REF_DIR, "scripts", "make_summary_tables.py")
    # LyRic's generator writes to a hardcoded relative "summary" dir, so its cwd
    # must be the directory that contains summary/ (derived from the output path
    # so an overridden --lyric/--lyric-joint still generates to the right place).
    lyric_cwd = os.path.dirname(os.path.dirname(os.path.abspath(args.lyric)))
    lyric_joint_cwd = os.path.dirname(os.path.dirname(os.path.abspath(args.lyric_joint)))
    # busco_references generator reads/writes via -s/--summary-dir, so it is cwd-independent.
    ref_gen = {"cwd": BUSCO_REF_DIR, "cmd": [py, ref_script, "-s", os.path.dirname(os.path.abspath(args.reference))]}
    iso_gen = {"cwd": ISO_DIR, "cmd": [py, iso_script, os.path.dirname(os.path.abspath(args.isoquant))]}
    lyric_gen = {"cwd": lyric_cwd, "cmd": [py, lyric_script]}
    lyric_merged_gen = {"cwd": lyric_cwd, "cmd": [py, lyric_script, "--merged"]}
    lyric_joint_gen = {"cwd": lyric_joint_cwd, "cmd": [py, lyric_script, "--joint"]}
    geneid_gen = {"cwd": GENEID_DIR, "cmd": [py, geneid_script, "-r", os.path.abspath(args.geneid_results)]}

    # (label, path, generator, taxid_col) for every source; geneid flavours expand to one source each.
    # each group is gated by its section so deactivated sections add no columns and
    # never trigger their generator. reference goes first so base genome/annotation
    # metrics precede the tool columns; its species key carries no taxid, so it
    # matches via the separate "taxonID" column.
    sources = []
    if "reference" in active:
        sources.append(("reference", args.reference, ref_gen, "taxonID"))
    if "isoquant" in active:
        sources.append(("isoquant", args.isoquant, iso_gen, None))
    if "lyric" in active:
        # use the unified joint table if it exists or was requested, else the two separate ones
        if os.path.isfile(args.lyric_joint) or args.lyric_mode == "joint":
            if os.path.isfile(args.lyric_joint):
                print(f"LyRic: using unified table {args.lyric_joint}")
            sources.append(("lyric", args.lyric_joint, lyric_joint_gen, None))
        else:
            sources.append(("lyric", args.lyric, lyric_gen, None))
            sources.append(("lyric_merged", args.lyric_merged, lyric_merged_gen, None))
    if "geneid" in active:
        for flav in args.geneid_flavours:
            sources.append((
                f"geneid_{flav}",
                os.path.join(args.geneid_results, "summary", flav, "general_summary.tsv"),
                geneid_gen,  # one run produces every flavour; deduped by run_generator
                None,
            ))

    total = 0
    ran = set()  # dedupe generator runs across sources (e.g. geneid's single command)
    for label, path, gen, taxid_col in sources:
        total += attach_source(protists, taxid_index, name_index, name_keys, name_by_idx,
                               label, path, gen=gen, generate=args.generate,
                               refresh=args.refresh, ran=ran, taxid_col=taxid_col)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    protists.to_csv(args.out, sep="\t", index=False, na_rep="NA")
    print(f"\nWrote {args.out} ({len(protists)} rows, {len(protists.columns)} columns; {total} total result matches)")

    if args.html:
        html_path = os.path.splitext(args.out)[0] + ".html"
        save_html(protists, html_path)
        print(f"Wrote {html_path}")


if __name__ == "__main__":
    main()
