#!/usr/bin/env bash

#Clone (or update) every engine listed in repos.json into the repo root.
#the engines are independent repositories, this script assembles them side by side for local development
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
MANIFEST="$ROOT_DIR/repos.json"

if [[ ! -f "$MANIFEST" ]]; then
  echo "Manifest not found: $MANIFEST" >&2
  exit 1
fi

#emit "name<TAB>url<TAB>branch" per repo, using jq when available and falling back to python3 otherwise.
read_manifest() {
  if command -v jq >/dev/null 2>&1; then
    jq -r '.repos[] | [.name, .url, .branch] | @tsv' "$MANIFEST"
  else
    python3 -c '
import json, sys
with open(sys.argv[1]) as fh:
    for r in json.load(fh)["repos"]:
        print("\t".join((r["name"], r["url"], r.get("branch", "main"))))
' "$MANIFEST"
  fi
}

while IFS=$'\t' read -r name url branch; do
  [[ -z "$name" ]] && continue
  dest="$ROOT_DIR/$name"
  if [[ -d "$dest/.git" ]]; then
    echo "Updating $name ($branch)"
    git -C "$dest" fetch --quiet origin "$branch"
    git -C "$dest" checkout --quiet "$branch"
    git -C "$dest" pull --quiet --ff-only origin "$branch"
  else
    echo "Cloning $name ($branch)"
    git clone --quiet --branch "$branch" "$url" "$dest"
  fi
done < <(read_manifest)

echo "Done."
