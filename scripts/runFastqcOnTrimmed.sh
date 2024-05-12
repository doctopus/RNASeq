#!/bin/bash

#SBATCH --job-name=FastQC_Trimmed
#SBATCH --ntasks=1
#SBATCH --partition=xtreme
#SBATCH --cpus-per-task=1  
#SBATCH --mem-per-cpu=15000
#SBATCH -o ./o/TrimmedFastQC_%A_%a.out
#SBATCH -e ./o/TrimmedFastQC_%A_%a.err

# start_time=$(date +%s)
# # Capture the default output filename specified by #SBATCH -o
# output_file=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')


echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load FastQC/0.11.9
module load python/2.7.18

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data
trimmedfqcdir=${base}/FastQC_Trimmed
mkdir -p ${trimmedfqcdir}

#TRIMMED_FASTQ_DIR="${base}/FASTQ_TRIMMED"

filename=${wdir}/trimmedFastq.list.txt

index=${SLURM_ARRAY_TASK_ID}

count=1
while IFS=" " read -r f1 f2 f3 f4
do
	if [ "$index" == "$count" ]; then
		fwdfq=$f3
		revfq=$f4

		fastqc ${fwdfq} --outdir ${trimmedfqcdir}
		fastqc ${revfq} --outdir ${trimmedfqcdir}
	fi

	echo ${count}

	count=$(expr $count + 1)

done  < ${filename}

# # End time
# end_time=$(date +%s)
# time_taken=$((end_time - start_time))
# # Convert time taken to human-readable format
# time_taken_formatted=$(date -u -d @${time_taken} +"%H:%M:%S")
# # Append the time taken to the default output file
# echo "Time taken by job name:$SLURM_JOB_NAME \
#     with job ID: $SLURM_JOB_ID \
#     having array task ID: $SLURM_ARRAY_TASK_ID :\
#     $time_taken_formatted" | tee -a "$output_file"

# #### When finished, there is a multiqc package to combine FastQC results
# #### module load python/3.9.5
# #### cd FastQC_Trimmed
# #### multiqc .