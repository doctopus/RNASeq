#!/bin/bash
# Building a STAR index file if not indexed before

#SBATCH --job-name=STARreference
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/STARIndex_%A_%a.out
#SBATCH -e ./o/STARIndex_%A_%a.err

######
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load python/2.7.18
module load STAR/2.7.6a

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data
##resources=/mnt/beegfs/training/CITIWorkshops/RNASeq/resources
resources=${pdir}/resources
#resoucesGTF=${resouurces}/GTF
#resourcesFASTA=${resources}/FASTA

mkdir -p ${resources}/STARv45p14

reference=${resources}/STARv45p14
#gtf=${resources}/GTF/gencode.v36.chr_patch_hapl_scaff.annotation.gtf
#filename=${wdir}/Fastq.list.txt

index=${SLURM_ARRAY_TASK_ID}
threads=16
count=1


# set the runThreadN to be one less than your NCPU request!
# sjdbOverhang should be readlength-1. The read length will be in fastqc file here it is 150
# sjdbOverhang value of 100 will work similarly to the ideal value
STAR --runThreadN 15 \
--runMode genomeGenerate \
--genomeDir $reference \
--genomeFastaFiles /mnt/beegfs/beherat2/kinase_RNASeq/resources/FASTA/GRCh38.p14.genome.fa \
--sjdbGTFfile /mnt/beegfs/beherat2/kinase_RNASeq/resources/GTF/gencode.v45.chr_patch_hapl_scaff.annotation.gtf \
--outFileNamePrefix $resources/Log.out \
#--sjdbOverhang 149