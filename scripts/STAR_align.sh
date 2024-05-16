#!/bin/bash

#SBATCH --job-name=STAR
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/STAR_%A_%a.out
#SBATCH -e ./o/STAR_%A_%a.err

######
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load python/2.7.18
module load STAR/2.7.6a
module load samtools/1.8

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data
##resources=/mnt/beegfs/training/CITIWorkshops/RNASeq/resources
resources=${pdir}/resources

mkdir -p ${base}/rawbam

reference=${resources}/STARv45p14
#gtf=${resources}/GTF/gencode.v36.chr_patch_hapl_scaff.annotation.gtf
gtf=${resources}/GTF/gencode.v45.chr_patch_hapl_scaff.annotation.gtf

filename=${wdir}/fastqFiles.list.txt
uniquesample=(318-Control_A1 318-Control_A2 318-Control_A3 \
				318-128-10_B1 318-128-10_B2 318-128-10_B3 \
				318-130_C1 318-130_C2 318-130_C3 \
				318-128-13_D1 318-128-13_D2 318-128-13_D3 \
				358-Control_E1 358-Control_E2 358-Control_E3 \
				358-128-10_F1 358-128-10_F2 358-128-10_F3 \
				358-130_G1 358-130_G2 358-130_G3 \
				358-128-13_H1 358-128-13_H2 358-128-13_H3)
index=${SLURM_ARRAY_TASK_ID}  
wrkgdir=${base}/rawbam
threads=16
count=1

while IFS="	" read -r f1 f2 f3 f4
do
	if [ "$index" == "$count" ]; then
		fastqfiles="$f3 $f4"
		wrkgdir_f="$wrkgdir/$f2."
                bamfile="${wrkgdir_f}Aligned.sortedByCoord.out.bam"
	fi

	STAR --runThreadN 16 \
	--outSAMattributes NM \
	--genomeDir ${reference} \
	--sjdbGTFfile ${gtf} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--twopassMode Basic \
	--readFilesCommand gunzip -c \
	--readFilesIn ${fastqfiles} \
	--outFileNamePrefix ${wrkgdir_f}

        samtools index ${bamfile}

	count=$(expr $count + 1)

done  < ${filename}


