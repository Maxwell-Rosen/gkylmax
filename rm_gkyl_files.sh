#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 [-n] [-f] [-r] <directory>
 -n  dry-run (list files, don't delete)
 -f  force (no confirmation)
 -r  recursive (descend into subdirectories). Without -r only files in the top-level directory are considered.
Example: $0 -n stellar-lorentzian1x-exploration
EOF
  exit 1
}

dry_run=0
force=0
recursive=0
while getopts "nfr" opt; do
  case $opt in
    n) dry_run=1;;
    f) force=1;;
    r) recursive=1;;
    *) usage;;
  esac
done
shift $((OPTIND-1))

if [ $# -lt 1 ]; then
  usage
fi

dir="$1"
[ -d "$dir" ] || { echo "Directory not found: $dir"; exit 2; }

# Build find command depending on recursive flag
if [ "$recursive" -eq 1 ]; then
  find_cmd=(find "$dir" -type f \( -name '*.gkyl' -o -name '*.json' \) -print0)
else
  find_cmd=(find "$dir" -maxdepth 1 -type f \( -name '*.gkyl' -o -name '*.json' \) -print0)
fi

# Collect files safely (handle spaces/newlines)
files=()
while IFS= read -r -d '' f; do
  files+=("$f")
done < <("${find_cmd[@]}")

if [ ${#files[@]} -eq 0 ]; then
  echo "No *.gkyl or *.json files found under $dir"
  exit 0
fi

# Compute total size in bytes
total_bytes=0
for f in "${files[@]}"; do
  # stat -c%s works on Linux; ignore errors by treating size as 0
  sz=$(stat -c%s -- "$f" 2>/dev/null || echo 0)
  total_bytes=$((total_bytes + sz))
done

# Convert to GB with 3 decimal places
total_gb=$(awk -v b="$total_bytes" 'BEGIN{printf "%.3f", b/1024/1024/1024}')

echo "Found ${#files[@]} files (total ${total_gb} GB):"
for f in "${files[@]}"; do printf '  %s\n' "$f"; done

if [ $dry_run -eq 1 ]; then
  echo "Dry-run: no files will be removed."
  exit 0
fi

if [ $force -eq 0 ]; then
  printf "Delete these %d files (total %s GB)? [y/N] " "${#files[@]}" "$total_gb"
  read -r ans
  case "$ans" in
    [yY]|[yY][eE][sS]) ;;
    *) echo "Aborted."; exit 0;;
  esac
fi

# Show size summary then delete
printf "Removing %s GB of files...\n" "$total_gb"
for f in "${files[@]}"; do
  rm -f -- "$f"
done

echo "Done."