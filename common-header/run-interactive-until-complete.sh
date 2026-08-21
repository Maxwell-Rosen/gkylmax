#!/bin/bash

# Repeatedly run a checkpointable Gkeyll simulation in Perlmutter interactive
# allocations until its final checkpoint exists.

set -uo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
script_path="$script_dir/$(basename -- "${BASH_SOURCE[0]}")"

show_help()
{
  cat <<'EOF'
Usage: run-interactive-until-complete.sh [OPTIONS] [RUN_DIR]

Run a checkpointable Gkeyll simulation in successive Perlmutter interactive
allocations. RUN_DIR defaults to the current directory.

Required:
  -p, --file-prefix PREFIX       Checkpoint prefix before the frame number
                                 (for example, gk_mirror-ion_)

Simulation options:
  -d, --run-dir DIR             Simulation directory (default: current directory)
  -f, --final-frame NUM         Completion frame (default: 65)
  -j, --job-script FILE         Script run inside each allocation, relative to
                                 RUN_DIR unless absolute
                                 (default: jobscript-gkyl-perlmutter)
      --foreground              Run in this shell instead of a detached tmux session
      --max-no-progress NUM     Stop after this many allocations without a new
                                 checkpoint; 0 disables the limit (default: 3)

Allocation options:
  -A, --account ACCOUNT         Slurm account (default: m4487)
  -N, --nodes NUM               Nodes per allocation (default: 2)
  -n, --ntasks NUM              Tasks per allocation (default: 64)
  -G, --gpus NUM                GPUs per allocation (default: 8)
  -t, --time D-HH:MM:SS         Allocation wall time (default: 04:00:00)
  -C, --constraint NAME         Slurm constraint (default: gpu)
  -q, --qos NAME                Slurm QOS (default: interactive)
      --session-limit NUM       Maximum concurrent interactive sessions
                                 (default: 2)
      --poll-seconds NUM        Seconds between session-slot checks (default: 600)

Other:
  -h, --help                    Show this help message

Examples:
  run-interactive-until-complete.sh -p gk_lorentzian_mirror-ion_
  run-interactive-until-complete.sh -p gk_mirror-ion_ -f 80 /path/to/run
  run-interactive-until-complete.sh -p gk_mirror-ion_ -A myproject -N 1 -G 4
EOF
}

missing_value()
{
  echo "Option $1 requires a value." >&2
  echo "Use --help for usage information." >&2
  exit 2
}

mode=start
run_dir_input=$PWD
run_dir_was_set=0
final_frame=65
file_prefix=
job_script_input=jobscript-gkyl-perlmutter
max_no_progress=3
interactive_session_limit=2
session_poll_seconds=600
account=m4487
nodes=2
ntasks=64
gpus=8
allocation_time=04:00:00
constraint=gpu
qos=interactive

while (( $# > 0 )); do
  case $1 in
    --worker)
      mode=worker
      shift
      ;;
    --foreground)
      mode=worker
      shift
      ;;
    -d|--run-dir)
      (( $# >= 2 )) || missing_value "$1"
      run_dir_input=$2
      run_dir_was_set=1
      shift 2
      ;;
    --run-dir=*)
      run_dir_input=${1#*=}
      run_dir_was_set=1
      shift
      ;;
    -f|--final-frame)
      (( $# >= 2 )) || missing_value "$1"
      final_frame=$2
      shift 2
      ;;
    --final-frame=*)
      final_frame=${1#*=}
      shift
      ;;
    -p|--file-prefix)
      (( $# >= 2 )) || missing_value "$1"
      file_prefix=$2
      shift 2
      ;;
    --file-prefix=*)
      file_prefix=${1#*=}
      shift
      ;;
    -j|--job-script)
      (( $# >= 2 )) || missing_value "$1"
      job_script_input=$2
      shift 2
      ;;
    --job-script=*)
      job_script_input=${1#*=}
      shift
      ;;
    --max-no-progress)
      (( $# >= 2 )) || missing_value "$1"
      max_no_progress=$2
      shift 2
      ;;
    --max-no-progress=*)
      max_no_progress=${1#*=}
      shift
      ;;
    -A|--account)
      (( $# >= 2 )) || missing_value "$1"
      account=$2
      shift 2
      ;;
    --account=*)
      account=${1#*=}
      shift
      ;;
    -N|--nodes)
      (( $# >= 2 )) || missing_value "$1"
      nodes=$2
      shift 2
      ;;
    --nodes=*)
      nodes=${1#*=}
      shift
      ;;
    -n|--ntasks)
      (( $# >= 2 )) || missing_value "$1"
      ntasks=$2
      shift 2
      ;;
    --ntasks=*)
      ntasks=${1#*=}
      shift
      ;;
    -G|--gpus)
      (( $# >= 2 )) || missing_value "$1"
      gpus=$2
      shift 2
      ;;
    --gpus=*)
      gpus=${1#*=}
      shift
      ;;
    -t|--time)
      (( $# >= 2 )) || missing_value "$1"
      allocation_time=$2
      shift 2
      ;;
    --time=*)
      allocation_time=${1#*=}
      shift
      ;;
    -C|--constraint)
      (( $# >= 2 )) || missing_value "$1"
      constraint=$2
      shift 2
      ;;
    --constraint=*)
      constraint=${1#*=}
      shift
      ;;
    -q|--qos)
      (( $# >= 2 )) || missing_value "$1"
      qos=$2
      shift 2
      ;;
    --qos=*)
      qos=${1#*=}
      shift
      ;;
    --session-limit)
      (( $# >= 2 )) || missing_value "$1"
      interactive_session_limit=$2
      shift 2
      ;;
    --session-limit=*)
      interactive_session_limit=${1#*=}
      shift
      ;;
    --poll-seconds)
      (( $# >= 2 )) || missing_value "$1"
      session_poll_seconds=$2
      shift 2
      ;;
    --poll-seconds=*)
      session_poll_seconds=${1#*=}
      shift
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    --)
      shift
      if (( $# > 1 )); then
        echo "Only one positional RUN_DIR may be supplied." >&2
        exit 2
      fi
      if (( $# == 1 )); then
        if (( run_dir_was_set )); then
          echo "RUN_DIR was supplied both positionally and with --run-dir." >&2
          exit 2
        fi
        run_dir_input=$1
        run_dir_was_set=1
      fi
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      echo "Use --help for usage information." >&2
      exit 2
      ;;
    *)
      if (( run_dir_was_set )); then
        echo "Only one RUN_DIR may be supplied." >&2
        exit 2
      fi
      run_dir_input=$1
      run_dir_was_set=1
      shift
      ;;
  esac
done

if [[ -z $file_prefix ]]; then
  echo "--file-prefix is required." >&2
  echo "Use --help for usage information." >&2
  exit 2
fi
if [[ ! $final_frame =~ ^[0-9]+$ ]]; then
  echo "--final-frame must be a non-negative integer, not '$final_frame'." >&2
  exit 2
fi
if [[ ! $max_no_progress =~ ^[0-9]+$ ]]; then
  echo "--max-no-progress must be a non-negative integer, not '$max_no_progress'." >&2
  exit 2
fi
for option_value in \
  "nodes:$nodes" \
  "ntasks:$ntasks" \
  "gpus:$gpus" \
  "session-limit:$interactive_session_limit" \
  "poll-seconds:$session_poll_seconds"; do
  option_name=${option_value%%:*}
  value=${option_value#*:}
  if [[ ! $value =~ ^[1-9][0-9]*$ ]]; then
    echo "--$option_name must be a positive integer, not '$value'." >&2
    exit 2
  fi
done
final_frame=$((10#$final_frame))
max_no_progress=$((10#$max_no_progress))
nodes=$((10#$nodes))
ntasks=$((10#$ntasks))
gpus=$((10#$gpus))
interactive_session_limit=$((10#$interactive_session_limit))
session_poll_seconds=$((10#$session_poll_seconds))
for option_value in \
  "account:$account" \
  "time:$allocation_time" \
  "constraint:$constraint" \
  "qos:$qos" \
  "job-script:$job_script_input"; do
  option_name=${option_value%%:*}
  value=${option_value#*:}
  if [[ -z $value ]]; then
    echo "--$option_name must not be empty." >&2
    exit 2
  fi
done
if [[ ! -d $run_dir_input ]]; then
  echo "Run directory does not exist: $run_dir_input" >&2
  exit 2
fi

run_dir=$(cd -- "$run_dir_input" && pwd -P)
if [[ $job_script_input == /* ]]; then
  job_script=$job_script_input
else
  job_script="$run_dir/$job_script_input"
fi
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
    echo "tmux is not available. Add --foreground to run in this shell." >&2
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

  worker_args=(
    --worker
    --run-dir "$run_dir"
    --final-frame "$final_frame"
    --file-prefix "$file_prefix"
    --job-script "$job_script"
    --max-no-progress "$max_no_progress"
    --account "$account"
    --nodes "$nodes"
    --ntasks "$ntasks"
    --gpus "$gpus"
    --time "$allocation_time"
    --constraint "$constraint"
    --qos "$qos"
    --session-limit "$interactive_session_limit"
    --poll-seconds "$session_poll_seconds"
  )
  printf -v worker_command '%q ' "$script_path" "${worker_args[@]}"
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
  echo "salloc is not available; run this command on a Perlmutter login node." >&2
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

echo "[$(date --iso-8601=seconds)] Interactive runner started"
echo "Run directory: $run_dir"
echo "Completion checkpoint: ${file_prefix}${final_frame}.gkyl"
echo "Job script: $job_script"
echo "Interactive-session limit: $interactive_session_limit; poll interval: ${session_poll_seconds}s"

interactive_session_count()
{
  local count

  if ! count=$(squeue \
    --me \
    --qos="$qos" \
    --states=PENDING,RUNNING \
    --noheader \
    --format='%i' | awk 'NF { count++ } END { print count + 0 }'); then
    echo "Could not query active $qos sessions with squeue." >&2
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
  echo "Requesting $nodes node(s), $gpus GPU(s), and $allocation_time of wall time ..."

  salloc \
    --nodes="$nodes" \
    --ntasks="$ntasks" \
    --constraint="$constraint" \
    --gpus="$gpus" \
    --account="$account" \
    --time="$allocation_time" \
    --qos="$qos" \
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

  allocation_timed_out=0
  if grep -Eqi \
    'DUE TO TIME LIMIT|exceeded its time limit|CANCELLED.*TIME LIMIT' \
    "$attempt_log"; then
    allocation_timed_out=1
    echo "The allocation reached its wall-time limit; the run will resume from its latest checkpoint."
  fi

  # salloc can return zero after Slurm revokes an allocation at its wall-time
  # limit, even though the srun job step was killed. Only treat a zero status as
  # a clean early exit when the allocation log does not show a timeout.
  if (( allocation_status == 0 && allocation_timed_out == 0 )); then
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
      echo "Set --max-no-progress 0 to disable this safety limit." >&2
      exit 1
    fi
  fi

  echo "The run is incomplete; requesting its next interactive allocation."
done

echo "Simulation complete (frame $final_frame exists)."
