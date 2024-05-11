#!/bin/bash

#SBATCH --job-name=FastQC
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2  
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/FastQC_%A_%a.out
#SBATCH -e ./o/FastQC_%A_%a.err

start_time=$(date +%s)

echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load FastQC/0.11.9
module load python/2.7.18

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/output
input=${pdir}/input
fqcdir=${base}/FastQC
mkdir -p ${fqcdir}

filename=${input}/fastqFiles.list.txt
index=${SLURM_ARRAY_TASK_ID}

count=1
while IFS="	" read -r f1 f2 f3 f4
do
	if [ "$index" == "$count" ]; then
		fwdfq=$f3
		revfq=$f4

		fastqc ${fwdfq} --outdir ${fqcdir}
		fastqc ${revfq} --outdir ${fqcdir}
	fi

	echo ${count}

	count=$(expr $count + 1)

done  < ${filename}

end_time=$(date +%s)
time_taken=$((end_time - start_time))
time_taken_formatted=$(date -u -d @${time_taken} +"%H:%M:%S")
# Capture the default output filename specified by #SBATCH -o
output_file=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
echo "Time taken by $SLURM_JOB_NAME: $time_taken_formatted" | tee -a "$output_file"

# #### When finished, there is a multiqc package to combine FastQC results
# #### module load python/3.9.5
# #### cd FastQC
# #### multiqc .