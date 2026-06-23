#!/usr/bin/env python3
"""
rename_to_longread_taxid.py

Rename the species folders in a data dir (default: species, i.e. run from data/) so that each
folder's trailing taxid suffix becomes the taxid that longread_protists.tsv
assigns to that organism.

Why a bridge other than the taxid is needed
--------------------------------------------
The whole point is that the folder taxid and the longread taxid disagree, so we
cannot join on the taxid itself. We join on the *organism name*:

  longread_protists.tsv  col5 = organism name   col3 = taxid

Matching policy (most specific first, never guess):
  1. STRAIN match  - the folder's full reconstructed name (binomial + strain,
                     minus the trailing taxid and the word "strain") matches a
                     longread organism name exactly (normalized). This keeps
                     strain-level taxids when longread actually has strain rows
                     (e.g. Toxoplasma_gondii_ME49 -> 508771).
  2. SPECIES match - fall back to the binomial (genus species), but ONLY if that
                     binomial maps to a single taxid in longread. This fixes
                     e.g. Plasmodium_berghei_ANKA_5823 -> _5821.
  3. SKIP          - no match, or the binomial is ambiguous and no strain match
                     was found. Reported, never renamed.

Safety:
  * dry-run by default; pass --apply to actually rename.
  * refuses any rename whose target already exists or collides with another
    planned target.
  * writes an undo script so every applied rename can be reverted.
"""
import argparse
import os
import re
import sys
from collections import defaultdict


def norm(s: str) -> str:
    """lowercase, non-alphanumerics -> single spaces, drop the token 'strain'."""
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", " ", s).strip()
    toks = [t for t in s.split() if t != "strain"]
    return " ".join(toks)


TAXID_SUFFIX = re.compile(r"_(\d+)$")


def split_folder(name: str):
    """return (base_without_taxid, taxid_or_None)."""
    m = TAXID_SUFFIX.search(name)
    if m:
        return name[: m.start()], m.group(1)
    return name, None


def load_longread(path: str):
    """build name->taxid (strain level) and binomial->set(taxids)."""
    name2tax = {}
    name_conflict = set()
    bino2tax = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            taxid, organism = cols[2].strip(), cols[4].strip()
            if not re.fullmatch(r"\d+", taxid):
                continue
            n = norm(organism)
            if not n:
                continue
            if n in name2tax and name2tax[n] != taxid:
                name_conflict.add(n)
            else:
                name2tax[n] = taxid
            toks = n.split()
            if len(toks) >= 2:
                bino2tax[" ".join(toks[:2])].add(taxid)
    for n in name_conflict:
        # ambiguous full name: refuse to use it
        name2tax.pop(n, None)
    return name2tax, bino2tax


def resolve(base: str, name2tax: dict, bino2tax: dict):
    """return (new_taxid, match_type) or (None, reason)."""
    fn = norm(base)
    if fn in name2tax:
        return name2tax[fn], "strain"
    toks = fn.split()
    if len(toks) >= 2:
        bino = " ".join(toks[:2])
        taxids = bino2tax.get(bino)
        if taxids and len(taxids) == 1:
            return next(iter(taxids)), "species"
        if taxids and len(taxids) > 1:
            return None, "ambiguous-binomial:" + "/".join(sorted(taxids))
    return None, "no-match"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-dir", default="species")
    ap.add_argument("--longread", default="longread_protists.tsv")
    ap.add_argument("--apply", action="store_true",
                    help="perform the renames (default: dry-run)")
    ap.add_argument("--undo-script", default="rename_undo.sh")
    args = ap.parse_args()

    if not os.path.isdir(args.species_dir):
        sys.exit(f"ERROR: not a directory: {args.species_dir}")
    if not os.path.isfile(args.longread):
        sys.exit(f"ERROR: no such file: {args.longread}")

    name2tax, bino2tax = load_longread(args.longread)

    folders = sorted(d for d in os.listdir(args.species_dir)
                     if os.path.isdir(os.path.join(args.species_dir, d)))

    plans = []          # (folder, new_name, old_taxid, new_taxid, match_type)
    skips = []          # (folder, reason)
    for folder in folders:
        base, old_tax = split_folder(folder)
        new_tax, info = resolve(base, name2tax, bino2tax)
        if new_tax is None:
            skips.append((folder, info))
            continue
        new_name = f"{base}_{new_tax}"
        if new_name == folder:
            plans.append((folder, new_name, old_tax, new_tax, info + " (nochange)"))
        else:
            plans.append((folder, new_name, old_tax, new_tax, info))

    # collision detection
    targets = defaultdict(list)
    for folder, new_name, *_ in plans:
        targets[new_name].append(folder)
    existing = set(folders)
    blocked = set()
    for folder, new_name, *_ in plans:
        if new_name == folder:
            continue
        if len(targets[new_name]) > 1:
            blocked.add(folder)
        elif new_name in existing:
            blocked.add(folder)

    # report
    print(f"{'FOLDER':<48} {'OLD':>9} -> {'NEW':<9} {'MATCH':<10} ACTION")
    print("-" * 92)
    to_apply = []
    for folder, new_name, old_tax, new_tax, info in plans:
        if new_name == folder:
            action = "keep"
        elif folder in blocked:
            action = "BLOCKED(collision)"
        else:
            action = "rename"
            to_apply.append((folder, new_name))
        print(f"{folder:<48} {str(old_tax):>9} -> {new_tax:<9} {info.split(' ')[0]:<10} {action}")
    for folder, reason in skips:
        print(f"{folder:<48} {'?':>9} -> {'?':<9} {'-':<10} SKIP ({reason})")

    print("-" * 92)
    print(f"folders: {len(folders)} | renames: {len(to_apply)} | "
          f"keep: {sum(1 for p in plans if p[0]==p[1])} | "
          f"blocked: {len(blocked)} | skipped: {len(skips)}")

    if not args.apply:
        print("\nDRY-RUN. Re-run with --apply to perform the renames above.")
        return

    if not to_apply:
        print("\nNothing to rename.")
        return

    with open(args.undo_script, "w", encoding="utf-8") as undo:
        undo.write("#!/bin/bash\n# undo log for rename_to_longread_taxid.py\nset -e\n")
        for folder, new_name in to_apply:
            src = os.path.join(args.species_dir, folder)
            dst = os.path.join(args.species_dir, new_name)
            os.rename(src, dst)
            undo.write(f'mv {dst!r} {src!r}\n')
            print(f"renamed: {folder} -> {new_name}")
    os.chmod(args.undo_script, 0o755)
    print(f"\nApplied {len(to_apply)} renames. Undo with: bash {args.undo_script}")


if __name__ == "__main__":
    main()
