#!/bin/bash

#==============================================================================
# submit-restarts-stellar-ai.sh
# 
# Submit a chain of dependent jobs to the SLURM scheduler. Each job depends
# on the previous job completing, allowing for automatic restart chains.
#==============================================================================

set -o pipefail

# Default values
num_jobs=1
restart_from=""
jobscript="jobscript-gkyl-stellar-ai"

#==============================================================================
# Help function
#==============================================================================
show_help() {
  cat <<EOF
Usage: submit-restarts-stellar-ai.sh [OPTIONS]

Options:
  --jobs NUM, -j NUM              Number of jobs to submit in the chain
                                  (default: 1)
  
  --restart-from JOB_ID, -r JOB_ID
                                  Start the chain from an existing job ID.
                                  The first new job will depend on JOB_ID.
  
  --jobscript FILE                Path to jobscript to submit
                                  (default: jobscript-gkyl-stellar-ai)
  
  --help, -h                      Show this help message and exit
EOF
}

#==============================================================================
# Argument parsing
#==============================================================================
while [[ $# -gt 0 ]]; do
  case $1 in
    --jobs|-j)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 1; }
      num_jobs="$2"
      shift 2
      ;;
    --restart-from|-r)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 1; }
      restart_from="$2"
      shift 2
      ;;
    --jobscript)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 1; }
      jobscript="$2"
      shift 2
      ;;
    --help|-h)
      show_help
      exit 0
      ;;
    *)
      echo "Error: Unknown option: $1" >&2
      echo "Use --help for usage information." >&2
      exit 1
      ;;
  esac
done

#==============================================================================
# Validation
#==============================================================================
if ! [[ "$num_jobs" =~ ^[0-9]+$ ]] || [ "$num_jobs" -lt 1 ]; then
  echo "Error: Number of jobs must be a positive integer." >&2
  exit 1
fi

if [ -n "$restart_from" ] && ! [[ "$restart_from" =~ ^[0-9]+$ ]]; then
  echo "Error: Job ID must be a positive integer." >&2
  exit 1
fi

if [ ! -f "$jobscript" ]; then
  echo "Error: Jobscript not found: $jobscript" >&2
  exit 1
fi

#==============================================================================
# Submit jobs
#==============================================================================
job_ids=()
previous_job_id="$restart_from"

if [[ -n "$previous_job_id" ]]; then
  echo "Starting chain from job: $previous_job_id"
fi

for ((i = 1; i <= num_jobs; i++)); do
  submit_args=(--parsable)
  if [[ -n "$previous_job_id" ]]; then
    submit_args+=("--dependency=afterany:$previous_job_id")
  fi
  if ! submission=$(sbatch "${submit_args[@]}" "$jobscript"); then
    echo "Error: Failed to submit job $i. Previously submitted jobs: ${job_ids[*]:-none}" >&2
    exit 1
  fi
  # Slurm may append a cluster name to parsable output.
  next_job_id=${submission%%;*}
  if ! [[ "$next_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Unexpected sbatch output: $submission" >&2
    exit 1
  fi
  job_ids+=("$next_job_id")
  echo "Submitted job $i: $next_job_id${previous_job_id:+ (depends on $previous_job_id)}"
  previous_job_id=$next_job_id
done

echo ""
echo "Successfully submitted job chain:"
echo "${job_ids[*]}"
