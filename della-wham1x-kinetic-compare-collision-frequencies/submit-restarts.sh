#!/bin/bash

# Array to store job IDs
job_ids=()

# Submit the first job and capture its job ID
first_job_id=$(sbatch jobscript-gkyl-della | awk '{print $4}')
job_ids+=($first_job_id)

# Submit subsequent jobs with dependency on the previous job
for i in {2..6} # Adjust the range as needed for the number of jobs
do
  previous_job_id=${job_ids[-1]}
  next_job_id=$(sbatch --dependency=afterany:$previous_job_id jobscript-gkyl-della | awk '{print $4}')
  job_ids+=($next_job_id)
done

# Print all job IDs
echo "Submitted jobs with IDs: ${job_ids[@]}"