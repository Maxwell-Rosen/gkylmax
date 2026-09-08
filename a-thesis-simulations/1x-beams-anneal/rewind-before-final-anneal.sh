#!/usr/bin/env bash

# Rewind an output directory made with the pre-final-anneal schedule to the
# end of frame 70. The obsolete final FDP occupied frames 71 through 85.

set -euo pipefail

apply=false
if [[ ${1:-} == --apply ]]; then
  apply=true
  shift
fi

if (( $# != 1 )); then
  echo "Usage: $0 [--apply] RUN_DIR" >&2
  exit 2
fi

run_dir_input=$1
app_name=${APP_NAME:-zzim}
keep_frame=70
old_last_frame=85

if [[ ! -d $run_dir_input ]]; then
  echo "Run directory does not exist: $run_dir_input" >&2
  exit 2
fi
run_dir=$(cd -- "$run_dir_input" && pwd -P)
if [[ $run_dir == / ]]; then
  echo "Refusing to operate on the filesystem root." >&2
  exit 2
fi

for species in ion elc; do
  checkpoint="$run_dir/${app_name}-${species}_${keep_frame}.gkyl"
  if [[ ! -f $checkpoint ]]; then
    echo "Required restart checkpoint is missing: $checkpoint" >&2
    exit 1
  fi
done

declare -a obsolete_frames=()
declare -a diagnostic_files=()

while IFS= read -r -d '' path; do
  name=${path##*/}
  if [[ $name =~ _([0-9]+)\.gkyl$ ]]; then
    frame=$((10#${BASH_REMATCH[1]}))
    if (( frame > keep_frame && frame <= old_last_frame )); then
      obsolete_frames+=("$path")
    fi
  else
    # These time-series files contain data from the obsolete FDP. Move them
    # aside so the restarted process does not append after that stale tail.
    case $name in
      "$app_name"-dt.gkyl|\
      "$app_name"-field_energy*.gkyl|\
      "$app_name"-*integrated*.gkyl|\
      "$app_name"-*omegaH_dt.gkyl|\
      "$app_name"-*L2norm.gkyl)
        diagnostic_files+=("$path")
        ;;
    esac
  fi
done < <(find "$run_dir" -maxdepth 1 -type f -name "$app_name-*.gkyl" -print0)

echo "Run directory: $run_dir"
echo "Keeping frames 0 through $keep_frame."
echo "Frame-indexed files to delete: ${#obsolete_frames[@]}"
printf '  %s\n' "${obsolete_frames[@]}"
echo "Diagnostic files to archive: ${#diagnostic_files[@]}"
printf '  %s\n' "${diagnostic_files[@]}"

if ! $apply; then
  echo
  echo "Dry run only. Re-run with --apply to perform the cleanup."
  exit 0
fi

if (( ${#obsolete_frames[@]} == 0 )); then
  echo "No obsolete frame files were found; nothing was deleted."
  exit 0
fi

if (( ${#diagnostic_files[@]} > 0 )); then
  archive_dir="$run_dir/diagnostics-before-final-anneal-$(date +%Y%m%dT%H%M%S)-$$"
  mkdir -- "$archive_dir"
  mv -- "${diagnostic_files[@]}" "$archive_dir/"
  echo "Archived old diagnostic time series in: $archive_dir"
fi

rm -- "${obsolete_frames[@]}"
echo "Deleted obsolete frame files. Restart from frame $keep_frame and run through frame 100."
