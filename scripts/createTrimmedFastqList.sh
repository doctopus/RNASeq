#!/bin/bash

#SBATCH --job-name=createTrimmedFastqList
#SBATCH --time=6-23:00:00
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/FastqTrimmedList_%A_%a.out
#SBATCH -e ./o/FastqTrimmedList_%A_%a.err

### cleaning the rawbam folder and getting Bamfiles file

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data

TRIMMED_FASTQ_DIR="${base}/FASTQ_TRIMMED"

#mkdir -p ${TRIMMED_FASTQ_DIR}

input_file="${wdir}/Fastq.list.txt"
output_file="${wdir}/trimmedFastq.list.txt"

count=0
while IFS= read -r line; do
    ((count++))

    # Extract fields from the line
    serial=$(echo "$line" | awk '{print $1}')
    SAMPLE_NAME=$(echo "$line" | awk '{print $2}')
    FASTQ_R1=$(echo "$line" | awk '{print $3}')
    FASTQ_R2=$(echo "$line" | awk '{print $4}')

    if [[ -z "$serial" ]]; then
        continue
    fi
    # Construct paths for trimmed files
    trimmed_fastq_r1="${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_1.trimmed_PE.fq.gz"
    trimmed_fastq_r2="${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_2.trimmed_PE.fq.gz"

    # Output the new list with trimmed file paths
    echo -e "${serial} ${SAMPLE_NAME} ${trimmed_fastq_r1} ${trimmed_fastq_r2}" >> "${output_file}"

done < "${input_file}"

# Append an empty line
echo >> "${output_file}"