#!/bin/bash

#SBATCH --job-name=Trimmomatic
#SBATCH --ntasks=25
#SBATCH --partition=xtreme
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=10000
#SBATCH -o ./o/Trimmomatic_%A_%a.out
#SBATCH -e ./o/Trimmomatic_%A_%a.err

start_time=$(date +%s)
# Capture the default output filename specified by #SBATCH -o
output_file=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load python/3.9.5
module load java/11.0.9
module load trimmomatic/0.39

wdir=$(pwd)
pdir=$(dirname "$wdir")
base=${pdir}/data
TRIMMED_FASTQ_DIR=${base}/FASTQ_TRIMMED

mkdir -p ${TRIMMED_FASTQ_DIR}

filename=${wdir}/fastqFiles.list.txt


export trimmomatic_dir=/cm/shared/apps/trimmomatic/0.39

count=0
while IFS= read -r line; do
    ((count++))
	
    # Extract fields from the line
    serial=$(echo "$line" | awk '{print $1}')
    SAMPLE_NAME=$(echo "$line" | awk '{print $2}')
    FASTQ_R1=$(echo "$line" | awk '{print $3}')
    FASTQ_R2=$(echo "$line" | awk '{print $4}')

    java -jar ${trimmomatic_dir}/trimmomatic-0.39.jar \
        PE \
        -threads 25 \
        -phred33 \
        ${FASTQ_R1} ${FASTQ_R2} \
        ${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_1.trimmed_PE.fq.gz \
        ${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_1.trimmed_SE.fq.gz \
        ${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_2.trimmed_PE.fq.gz \
        ${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_2.trimmed_SE.fq.gz \
        -trimlog ${base}/trimmomaticLog.txt \
        ILLUMINACLIP:${trimmomatic_dir}/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:50
done < ${filename}

#Calculate Time Taken by the Script
end_time=$(date +%s)
time_taken=$((end_time - start_time))
time_taken_formatted=$(date -u -d @${time_taken} +"%H:%M:%S")
echo "Time Taken: $time_taken_formatted" | tee -a "$output_file"
