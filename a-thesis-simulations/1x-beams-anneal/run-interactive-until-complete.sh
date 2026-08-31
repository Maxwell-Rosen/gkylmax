#!/bin/bash

# Repeatedly run a checkpointable Gkeyll simulation in four-hour Perlmutter
# interactive allocations. With no arguments, this script runs the simulation
# in its own directory and expects its final checkpoint to be frame 65.
#
# It can also orchestrate another run directory:
#   ./run-interactive-until-complete.sh RUN_DIR FINAL_FRAME [FILE_PREFIX]

set -uo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
script_path="$script_dir/$(basename -- "${BASH_SOURCE[0]}")"

mode=start
if [[ ${1:-} == "--worker" ]]; then
  mode=worker
  shift
elif [[ ${1:-} == "--foreground" ]]; then
  mode=worker
  shift
fi

run_dir_input=${1:-$script_dir}
final_frame=${2:-85}
file_prefix=${3:-gk_lorentzian_mirror-ion_}
job_script_name=${JOB_SCRIPT_NAME:-jobscript-gkyl-perlmutter}
max_no_progress=${MAX_NO_PROGRESS:-3}
interactive_session_limit=${INTERACTIVE_SESSION_LIMIT:-2}
session_poll_seconds=${SESSION_POLL_SECONDS:-600}

if [[ ! $final_frame =~ ^[0-9]+$ ]]; then
  echo "FINAL_FRAME must be a non-negative integer, not '$final_frame'." >&2
  exit 2
fi
if [[ ! $max_no_progress =~ ^[0-9]+$ ]]; then
  echo "MAX_NO_PROGRESS must be a non-negative integer, not '$max_no_progress'." >&2
  exit 2
fi
if [[ ! $interactive_session_limit =~ ^[1-9][0-9]*$ ]]; then
  echo "INTERACTIVE_SESSION_LIMIT must be a positive integer, not '$interactive_session_limit'." >&2
  exit 2
fi
if [[ ! $session_poll_seconds =~ ^[1-9][0-9]*$ ]]; then
  echo "SESSION_POLL_SECONDS must be a positive integer, not '$session_poll_seconds'." >&2
  exit 2
fi
if [[ ! -d $run_dir_input ]]; then
  echo "Run directory does not exist: $run_dir_input" >&2
  exit 2
fi

run_dir=$(cd -- "$run_dir_input" && pwd -P)
job_script="$run_dir/$job_script_name"
completion_file="$run_dir/${file_prefix}${final_frame}.gkyl"
master_log="$run_dir/interactive-run.log"

mkdir -p -- "$run_dir/.interactive-logs"
exec > >(tee -a "$master_log") 2>&1

latest_frame()
{
  local path name number
  local newest=-1
  local -a checkpoint_files=()

  shopt -s nullglob
  checkpoint_files=("$run_dir/${file_prefix}"[0-9]*.gkyl)
  shopt -u nullglob

  for path in "${checkpoint_files[@]}"; do
    name=${path##*/}
    number=${name#"$file_prefix"}
    number=${number%.gkyl}
    if [[ $number =~ ^[0-9]+$ ]] && (( 10#$number > newest )); then
      newest=$((10#$number))
    fi
  done
  printf '%d\n' "$newest"
}

is_complete()
{
  local newest
  [[ -f $completion_file ]] && return 0
  newest=$(latest_frame)
  (( newest >= final_frame ))
}

if [[ $mode == start ]]; then
  if ! command -v tmux >/dev/null 2>&1; then
    echo "tmux is not available. Run this in the foreground with:" >&2
    echo "  $script_path --foreground '$run_dir' '$final_frame' '$file_prefix'" >&2
    exit 1
  fi

  if is_complete; then
    echo "Simulation is already complete (frame $final_frame exists)."
    exit 0
  fi

  run_tag=$(basename -- "$run_dir" | tr -c '[:alnum:]_-' '-')
  path_tag=$(printf '%s' "$run_dir" | cksum | awk '{print $1}')
  session_name="gkyl-${run_tag:0:24}-${path_tag}"

  if tmux has-session -t "=$session_name" 2>/dev/null; then
    echo "Already running in tmux session: $session_name"
    echo "Attach with: tmux attach -t '$session_name'"
    exit 0
  fi

  printf -v worker_command '%q ' \
    env "MAX_NO_PROGRESS=$max_no_progress" "JOB_SCRIPT_NAME=$job_script_name" \
    "INTERACTIVE_SESSION_LIMIT=$interactive_session_limit" \
    "SESSION_POLL_SECONDS=$session_poll_seconds" \
    "$script_path" --worker "$run_dir" "$final_frame" "$file_prefix"
  if ! tmux new-session -d -s "$session_name" "$worker_command"; then
    echo "Could not create tmux session '$session_name'." >&2
    exit 1
  fi

  echo "Started: $session_name"
  echo "Attach:  tmux attach -t '$session_name'"
  echo "Log:     $master_log"
  exit 0
fi

if ! command -v salloc >/dev/null 2>&1; then
  echo "salloc is not available; run this script on a Perlmutter login node." >&2
  exit 1
fi
if ! command -v squeue >/dev/null 2>&1; then
  echo "squeue is not available; cannot check for an interactive-session slot." >&2
  exit 1
fi
if [[ ! -r $job_script ]]; then
  echo "Cannot read job script: $job_script" >&2
  exit 1
fi
if [[ ! -x $run_dir/sim ]]; then
  echo "Simulation executable is missing or not executable: $run_dir/sim" >&2
  exit 1
fi

echo "[$(date --iso-8601=seconds)] Interactive runner started"
echo "Run directory: $run_dir"
echo "Completion checkpoint: ${file_prefix}${final_frame}.gkyl"
echo "Interactive-session limit: $interactive_session_limit; poll interval: ${session_poll_seconds}s"

interactive_session_count()
{
  local count

  if ! count=$(squeue \
    --me \
    --qos=interactive \
    --states=PENDING,RUNNING \
    --noheader \
    --format='%i' | awk 'NF { count++ } END { print count + 0 }'); then
    echo "Could not query active interactive sessions with squeue." >&2
    return 1
  fi

  printf '%d\n' "$count"
}

wait_for_interactive_session()
{
  local active_sessions

  while true; do
    if ! active_sessions=$(interactive_session_count); then
      return 1
    fi

    if (( active_sessions < interactive_session_limit )); then
      echo "[$(date --iso-8601=seconds)] Interactive slot available ($active_sessions/$interactive_session_limit in use)."
      return 0
    fi

    echo "[$(date --iso-8601=seconds)] All $interactive_session_limit interactive slots are in use."
    echo "Polling again in $session_poll_seconds seconds ..."
    sleep "$session_poll_seconds"
  done
}

attempt=0
no_progress=0
run_id="$(date +%Y%m%dT%H%M%S)-$$"

while ! is_complete; do
  if ! wait_for_interactive_session; then
    exit 1
  fi

  attempt=$((attempt + 1))
  before=$(latest_frame)
  attempt_log="$run_dir/.interactive-logs/${run_id}-attempt-${attempt}.log"

  echo
  echo "[$(date --iso-8601=seconds)] Attempt $attempt; latest frame: $before"
  echo "Requesting a two-node, eight-GPU, four-hour interactive allocation ..."

  salloc \
    --nodes=2 \
    --ntasks=64 \
    --constraint=gpu \
    --gpus=8 \
    --account=m4487 \
    --time=04:00:00 \
    --qos=interactive \
    --job-name="gkyl-${run_id:0:24}" \
    --chdir="$run_dir" \
    bash -l "$job_script" 2>&1 | tee "$attempt_log"
  allocation_status=${PIPESTATUS[0]}

  after=$(latest_frame)
  echo "[$(date --iso-8601=seconds)] Allocation ended with status $allocation_status; latest frame: $after"

  if is_complete; then
    echo "Simulation complete. The allocation was released immediately."
    exit 0
  fi

  if grep -Eq \
    'Update method failed|Time-step was below .*Aborting simulation|Failed to read restart file' \
    "$attempt_log"; then
    echo "A numerical or restart failure was detected; not requesting another allocation." >&2
    echo "See: $attempt_log" >&2
    exit 1
  fi

  if (( allocation_status == 0 )); then
    echo "The simulation command exited normally without producing final frame $final_frame." >&2
    echo "See: $attempt_log" >&2
    exit 1
  fi

  if (( after > before )); then
    no_progress=0
  else
    no_progress=$((no_progress + 1))
    echo "No new checkpoint was produced in this allocation ($no_progress consecutive)."
    if (( max_no_progress > 0 && no_progress >= max_no_progress )); then
      echo "Stopping after $max_no_progress allocations without checkpoint progress." >&2
      echo "Set MAX_NO_PROGRESS=0 to disable this safety limit." >&2
      exit 1
    fi
  fi

  echo "The run is incomplete; requesting its next interactive allocation."
done

echo "Simulation complete (frame $final_frame exists)."
