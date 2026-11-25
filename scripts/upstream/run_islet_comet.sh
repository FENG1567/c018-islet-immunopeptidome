#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 COMET_EXECUTABLE PARAMS_FILE OUTPUT_PEPXML_DIR INPUT.mgf [INPUT.mgf ...]" >&2
  exit 2
fi

comet=$(realpath "$1")
params=$(realpath "$2")
pep_dir=$3
shift 3
mkdir -p "$pep_dir"
pep_dir=$(realpath "$pep_dir")

for mgf in "$@"; do
  source=$(realpath "$mgf")
  stem=$(basename "$source" .mgf)
  work="$pep_dir/.${stem}.work"
  if [[ -e "$pep_dir/$stem.pep.xml" || -e "$work" ]]; then
    echo "Refusing to overwrite existing output for $stem" >&2
    exit 1
  fi
  mkdir -p "$work"
  printf 'SEARCH_START %s %s\n' "$stem" "$(date -Is)"
  (cd "$(dirname "$params")" && "$comet" "-P$params" "-N$work/$stem" "$source") \
    >"$pep_dir/$stem.comet.log" 2>&1
  test -s "$work/$stem.pep.xml"
  mv "$work/$stem.pep.xml" "$pep_dir/$stem.pep.xml"
  printf 'SEARCH_DONE %s %s %s\n' "$stem" "$(date -Is)" "$(stat -c %s "$pep_dir/$stem.pep.xml")"
done
