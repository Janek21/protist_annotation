#!/usr/bin/env python3
"""collect_busco_summaries.py

Find every BUSCO summary PNG produced by the annotation pipelines and copy each
one into a destination folder, organised by ORIGIN so nothing collides (the same
file name -- e.g. busco_figure.png -- appears in many pipelines).

Only the four result pipelines are searched, each a deactivatable SECTION
(mirrors result_analytics/build_results_table.py):
  reference   busco_references
  isoquant    isoquant_annotator
  lyric       lyric_annotator
  geneid      geneid-training
Deactivate sections when a tool's plots are not wanted:
  --skip geneid            drop the geneid section
  --skip geneid lyric      drop several at once
  --only reference         keep only one section (overrides --skip)

With --refresh, buscoPlot.py is re-run for each active section before collecting
(via that pipeline's summary generator, which calls buscoPlot for every folder)
so the gathered plots are freshly regenerated:
  reference   busco_references/scripts/make_summary_tables.py -s summary
  isoquant    isoquant_annotator/scripts/make_summary_tables.py summary
  lyric       lyric_annotator/scripts/make_summary_tables.py
  geneid      geneid-training/scripts/counting_machine.py -r results

"Origin" = where the PNG came from, relative to the projects root. Choose how it
is derived with --origin:
  dir     (default) full source directory path, sanitised  e.g.
                    lyric_annotator/summary/busco_lineage -> lyric_annotator_summary_busco_lineage
  top               first path component                 e.g. lyric_annotator
  parent            immediate parent directory name      e.g. busco_lineage

Layout produced (default --layout subdirs):
  <dest>/<origin>/<original_name>.png
With --layout flat:
  <dest>/<origin>__<original_name>.png

By default it matches pipeline summary plots and BUSCO's own figures:
  *summary*.png   busco_figure*.png
Override with one or more --glob.

Heavy pipeline-intermediate dirs that never hold BUSCO plots are pruned by name
(--prune-dir, default 'aux' and 'fastas_*') so the search stays fast on the
multi-terabyte tree; discovery uses the OS `find` with an os.scandir fallback.

Copy is additive and non-destructive. Use --dry-run to preview.
"""
import argparse
import os
import shutil
import subprocess
import sys
from fnmatch import fnmatch

# projects root = parent of this script's directory (result_analytics/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)

# Deactivatable sections (--skip / --only) -> the pipeline folder each searches.
SECTION_DIRS = {
    "reference": "busco_references",
    "isoquant":  "isoquant_annotator",
    "lyric":     "lyric_annotator",
    "geneid":    "geneid-training",
}
SECTIONS = ["reference", "isoquant", "lyric", "geneid"]


# How to re-run buscoPlot.py for each section. Each pipeline embeds the plotting
# in its summary generator (which calls buscoPlot.run() for every folder), so a
# refresh runs that generator -- the same commands result_analytics's
# build_results_table.py uses. cwd is the section dir; paths are relative to it.
def generator_for(section, sec_dir):
    py = sys.executable
    scripts = os.path.join(sec_dir, "scripts")
    if section == "reference":
        return [py, os.path.join(scripts, "make_summary_tables.py"), "-s", "summary"]
    if section == "isoquant":
        return [py, os.path.join(scripts, "make_summary_tables.py"), "summary"]
    if section == "lyric":
        return [py, os.path.join(scripts, "make_summary_tables.py")]
    if section == "geneid":
        return [py, os.path.join(scripts, "counting_machine.py"), "-r", "results"]
    return None


def run_generator(label, cmd, cwd):
    """Run one section's plot/summary generator (regenerates its BUSCO plots)."""
    print(f"[gen]  {label}: running {' '.join(cmd)} (cwd={cwd})")
    try:
        res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    except Exception as exc:  # report any launch failure and continue
        print(f"[gen]  {label}: error launching generator: {exc}")
        return
    if res.returncode != 0:
        tail = (res.stderr.strip().splitlines() or [""])[-1]
        print(f"[gen]  {label}: exited {res.returncode}: {tail}")
    else:
        tail = (res.stdout.strip().splitlines() or ["done"])[-1]
        print(f"[gen]  {label}: {tail}")


def sanitize(s: str) -> str:
    return s.strip(os.sep).replace(os.sep, "_") or "root"


def origin_of(rel_dir: str, mode: str) -> str:
    parts = [p for p in rel_dir.split(os.sep) if p not in ("", ".")]
    if not parts:
        return "root"
    if mode == "top":
        return parts[0]
    if mode == "parent":
        return parts[-1]
    return sanitize(rel_dir)            # mode == "dir"


def find_pngs(src, globs, dest, excludes, prune_dirs):
    """Locate matching PNGs under `src`.

    Prefers the OS `find` command: this tree can hold millions of files and a
    pure-Python walk is too slow (it gets watchdog-killed). Falls back to a lazy
    os.scandir walk when `find` is unavailable (e.g. Windows) or errors.

    `prune_dirs` are directory-name globs whose whole subtree is skipped. The
    defaults cut huge pipeline-intermediate stores (e.g. */OUT/aux, fastas_*)
    that never hold BUSCO summaries, turning a multi-minute walk into seconds.
    """
    matches = _find_via_os(src, globs, dest, excludes, prune_dirs)
    if matches is not None:
        return matches
    return _find_via_scandir(src, globs, dest, excludes, prune_dirs)


def _find_via_os(src, globs, dest, excludes, prune_dirs):
    # ( -type d ( -name .git -o -name <prune> ... -o -path <dest> -o -path *<exc>* ) -prune )
    # -o ( -type f ( -iname g1 -o -iname g2 ... ) -print )
    prune = ["-name", ".git", "-o", "-path", dest]
    for p in prune_dirs:
        prune += ["-o", "-name", p]
    for e in excludes:
        prune += ["-o", "-path", f"*{e}*"]
    name_tests = []
    for g in globs:
        if name_tests:
            name_tests.append("-o")
        name_tests += ["-iname", g]
    cmd = (["find", src, "(", "-type", "d", "("] + prune
           + [")", "-prune", ")", "-o", "(", "-type", "f", "("]
           + name_tests + [")", "-print", ")"])
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except (OSError, ValueError):
        return None
    if proc.returncode != 0:
        return None
    return [p for p in proc.stdout.split("\n") if p]


def _find_via_scandir(src, globs, dest, excludes, prune_dirs):
    matches = []
    stack = [src]
    while stack:
        cur = stack.pop()
        if os.path.abspath(cur) == dest or any(e in cur for e in excludes):
            continue
        try:
            scan = os.scandir(cur)
        except OSError:
            continue
        with scan:
            for entry in scan:
                try:
                    is_dir = entry.is_dir(follow_symlinks=False)
                except OSError:
                    is_dir = False
                if is_dir:
                    if (entry.name == ".git"
                            or os.path.abspath(entry.path) == dest
                            or any(fnmatch(entry.name, p) for p in prune_dirs)):
                        continue
                    stack.append(entry.path)
                elif entry.name.lower().endswith(".png") and any(
                        fnmatch(entry.name, g) for g in globs):
                    matches.append(entry.path)
    return matches


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=ROOT,
                    help="projects root holding the pipeline folders (default: "
                         "parent of this script)")
    ap.add_argument("--dest", default=os.path.join(SCRIPT_DIR, "busco_summaries"),
                    help="output folder (default: result_analytics/busco_summaries)")
    ap.add_argument("--glob", action="append", default=None,
                    help="filename glob to match (repeatable; "
                         "default: '*summary*.png' and 'busco_figure*.png')")
    ap.add_argument("--origin", choices=("dir", "top", "parent"), default="dir",
                    help="how to derive the origin name (default: dir)")
    ap.add_argument("--layout", choices=("subdirs", "flat"), default="subdirs",
                    help="subdirs: <dest>/<origin>/<name>; flat: <dest>/<origin>__<name>")
    ap.add_argument("--skip", nargs="+", default=[], choices=SECTIONS, metavar="SECTION",
                    help=f"deactivate whole sections; choices: {', '.join(SECTIONS)}. "
                         f"e.g. --skip geneid")
    ap.add_argument("--only", nargs="+", default=None, choices=SECTIONS, metavar="SECTION",
                    help="keep ONLY these sections (overrides --skip)")
    ap.add_argument("--exclude", action="append", default=[],
                    help="path substring to skip (repeatable), e.g. old_lyric, remove")
    ap.add_argument("--prune-dir", action="append", default=None,
                    help="directory-name glob whose whole subtree is skipped "
                         "(repeatable; default: 'aux' and 'fastas_*', the heavy "
                         "pipeline-intermediate dirs). Pass '' once to disable.")
    ap.add_argument("--refresh", action="store_true",
                    help="re-run buscoPlot.py for each active section first "
                         "(via that pipeline's summary generator) so the "
                         "collected plots are freshly regenerated")
    ap.add_argument("--dry-run", action="store_true", help="preview, copy nothing")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    dest = os.path.abspath(args.dest)
    globs = args.glob or ["*summary*.png", "busco_figure*.png"]
    prune_dirs = [p for p in (args.prune_dir if args.prune_dir is not None
                              else ["aux", "fastas_*"]) if p]

    # resolve which sections are active. --only wins over --skip when both given.
    if args.only is not None:
        active = [s for s in SECTIONS if s in set(args.only)]
    else:
        active = [s for s in SECTIONS if s not in set(args.skip)]
    if not active:
        ap.error("no sections active (everything was skipped); nothing to do")
    skipped = [s for s in SECTIONS if s not in active]
    print(f"Sections active: {', '.join(active)}"
          + (f" | deactivated: {', '.join(skipped)}" if skipped else ""))

    if args.refresh:
        for section in active:
            sec_dir = os.path.join(root, SECTION_DIRS[section])
            if not os.path.isdir(sec_dir):
                print(f"[skip] refresh {section}: not found at {sec_dir}")
                continue
            run_generator(section, generator_for(section, sec_dir), sec_dir)

    matches = []
    for section in active:
        sec_dir = os.path.join(root, SECTION_DIRS[section])
        if not os.path.isdir(sec_dir):
            print(f"[skip] {section}: not found at {sec_dir}")
            continue
        found = find_pngs(sec_dir, globs, dest, args.exclude, prune_dirs)
        print(f"[ok]   {section}: {len(found)} png under {SECTION_DIRS[section]}")
        matches.extend(found)

    if not matches:
        print("No BUSCO summary PNGs matched.", file=sys.stderr)
        return

    # build copy plan, resolving any within-origin name clashes. origin is taken
    # relative to the projects root so the section folder leads each origin.
    plan = []                      # (src_path, dst_path, origin)
    used = set()
    for path in sorted(matches):
        rel_dir = os.path.relpath(os.path.dirname(path), root)
        origin = origin_of(rel_dir, args.origin)
        name = os.path.basename(path)
        if args.layout == "flat":
            dst = os.path.join(dest, f"{origin}__{name}")
        else:
            dst = os.path.join(dest, origin, name)
        stem, ext = os.path.splitext(dst)
        i = 1
        while dst in used:                      # different sources, identical target
            dst = f"{stem}_{i}{ext}"
            i += 1
        used.add(dst)
        plan.append((path, dst, origin))

    by_origin = {}
    for _, _, origin in plan:
        by_origin[origin] = by_origin.get(origin, 0) + 1
    for origin in sorted(by_origin):
        print(f"  {origin:<45} {by_origin[origin]} png")
    print(f"{'[DRY-RUN] ' if args.dry_run else ''}"
          f"{len(plan)} PNGs from {len(by_origin)} origins -> {dest}")

    if args.dry_run:
        for s, d, _ in plan:
            print(f"  {os.path.relpath(s, root)}  ->  {os.path.relpath(d, dest)}")
        return

    for s, d, _ in plan:
        os.makedirs(os.path.dirname(d), exist_ok=True)
        shutil.copy2(s, d)
    print(f"Copied {len(plan)} files.")


if __name__ == "__main__":
    main()
