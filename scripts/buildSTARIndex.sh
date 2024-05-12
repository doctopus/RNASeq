#!/bin/bash
# Building a STAR index file

#SBATCH --job-name=STARIndex
#SBATCH --ntasks=25
#SBATCH --partition=xtreme
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH -o ./o/STARIndex_%A_%a.out
#SBATCH -e ./o/STARIndex_%A_%a.err

######
start_time=$(date +%s)
# Capture the default output filename specified by #SBATCH -o
output_file=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load python/2.7.18
module load STAR/2.7.6a

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data
resources=${pdir}/resources
#resoucesGTF=${resouurces}/GTF
#resourcesFASTA=${resources}/FASTA

mkdir -p ${resources}/STARIndex

reference=${resources}/STARIndex
#gtf=${resources}/GTF/gencode.v36.chr_patch_hapl_scaff.annotation.gtf
#filename=${wdir}/Fastq.list.txt  

# index=${SLURM_ARRAY_TASK_ID}  
# threads=16
# count=1


# set the runThreadN to be one less than your NCPU request!
# sjdbOverhang should be readlength-1. The read length will be in fastqc file here it is 150
# sjdbOverhang value of 100 will work similarly to the ideal value
STAR --runThreadN 25 \
--runMode genomeGenerate \
--genomeDir $reference \
--genomeFastaFiles /mnt/beegfs/beherat2/kinase_RNASeq/resources/FASTA/GRCh38.p14.genome.fa \
--sjdbGTFfile /mnt/beegfs/beherat2/kinase_RNASeq/resources/GTF/gencode.v45.chr_patch_hapl_scaff.annotation.gtf \
--outFileNamePrefix $resources/StarIndex_Log.out \
--sjdbOverhang 149 

# End time
end_time=$(date +%s)
time_taken=$((end_time - start_time))
# Convert time taken to human-readable format
time_taken_formatted=$(date -u -d @${time_taken} +"%H:%M:%S")
# Append the time taken to the default output file
echo "Time taken by job name:$SLURM_JOB_NAME \
    with job ID: $SLURM_JOB_ID \
    having array task ID: $SLURM_ARRAY_TASK_ID \
    : $time_taken_formatted" | tee -a "$output_file"
    
#TODO Fix Log.out file location. It throws warning in the err file that it could not transfer the log file