#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 [-n] [-f] [-r] <directory>
 -n  dry-run (list files, don't delete)
 -f  force (no confirmation)
 -r  recursive (descend into subdirectories). Without -r only files in the top-level directory are considered.
     Set DELETE_JOBS to control parallel delete workers (default: 4).
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

# Build find command depending on recursive flag.
if [ "$recursive" -eq 1 ]; then
  find_base=(find "$dir" -type f \( -name '*.gkyl' -o -name '*.json' \))
else
  find_base=(find "$dir" -maxdepth 1 -type f \( -name '*.gkyl' -o -name '*.json' \))
fi

# Use streaming approaches to avoid per-file stat loops and large arrays.
file_count=$("${find_base[@]}" -printf '.' | wc -c)
if [ "$file_count" -eq 0 ]; then
  echo "No *.gkyl or *.json files found under $dir"
  exit 0
fi

# Compute total size efficiently without materializing filenames.
total_bytes=$("${find_base[@]}" -printf '%s\n' | awk '{sum += $1} END {print sum + 0}')
total_gb=$(awk -v b="$total_bytes" 'BEGIN{printf "%.3f", b/1024/1024/1024}')

echo "Found ${file_count} files totaling about ${total_gb} GB."

if [ $dry_run -eq 1 ]; then
  echo "Dry-run: no files will be removed."
  exit 0
fi

if [ $force -eq 0 ]; then
  printf "Delete these %d files (total %s GB)? [y/N] " "$file_count" "$total_gb"
  read -r ans
  case "$ans" in
    [yY]|[yY][eE][sS]) ;;
    *) echo "Aborted."; exit 0;;
  esac
fi

# Delete efficiently using parallel batched rm processes.
# DELETE_JOBS can be tuned if the filesystem prefers fewer or more workers.
delete_jobs=${DELETE_JOBS:-16}
printf "Removing %s GB of files...\n" "$total_gb"
if [ "$recursive" -eq 1 ]; then
  find "$dir" -type f \( -name '*.gkyl' -o -name '*.json' \) -print0 | xargs -0 -r -n 2000 -P "$delete_jobs" rm -f --
else
  find "$dir" -maxdepth 1 -type f \( -name '*.gkyl' -o -name '*.json' \) -print0 | xargs -0 -r -n 2000 -P "$delete_jobs" rm -f --
fi

echo "Done."