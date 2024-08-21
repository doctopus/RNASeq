#!/bin/bash
#SBATCH --job-name=testTime
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/testTime_%A_%a.out
#SBATCH -e ./o/testTime_%A_%a.err

# Start time
start_time=$(date +%s)
# Capture the default output filename specified by #SBATCH -o
output_file=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

# Your job commands here
# For example:
sleep 65

# End time
end_time=$(date +%s)

# Calculate time taken in seconds
time_taken=$((end_time - start_time))

# Convert time taken to human-readable format
time_taken_formatted=$(date -u -d @${time_taken} +"%H:%M:%S")

# Append the time taken to the default output file
echo "Time taken by job name:$SLURM_JOB_NAME \
    with job ID: $SLURM_JOB_ID \
    having array task ID: $SLURM_ARRAY_TASK_ID \
    : $time_taken_formatted" | tee -a "$output_file"