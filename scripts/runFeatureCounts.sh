#!/bin/bash

#SBATCH --job-name=FeatureCounts
#SBATCH --ntasks=1
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=2 
#SBATCH --mem-per-cpu=15000
#SBATCH -o ./o/FeatureCounts_%A_%a.out
#SBATCH -e ./o/FeatureCounts_%A_%a.err

echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load python/2.7.18
module load R/4.3.2
module load subread/2.0.6

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data
##resources=/mnt/beegfs/training/CITIWorkshops/RNASeq/resources
resources=${pdir}/resources

index=${SLURM_ARRAY_TASK_ID}
##bamfiles=${wdir}/Bamfiles.txt
bamfiles=${base}/bamFiles.list.txt
# output=${base}/subread/featureCounts
output=${base}/FeatureCounts
mkdir -p ${output}
#gtf=${resources}/GTF/gencode.v36.chr_patch_hapl_scaff.annotation.gtf
gtf=${resources}/GTF/gencode.v45.chr_patch_hapl_scaff.annotation.gtf
format=GTF


### start array index at 2 cause we have column name
count=1
while read line
do
	if [ "$index" == "$count" ]; then
		/cm/shared/apps/subread/2.0.1/bin/featureCounts -p -B -s 0 -a ${gtf} -F ${format} -o ${output}/featureCounts_0.${index} ${line}
     	/cm/shared/apps/subread/2.0.1/bin/featureCounts -p -B -s 1 -a ${gtf} -F ${format} -o ${output}/featureCounts_s1.${index} ${line}
     	/cm/shared/apps/subread/2.0.1/bin/featureCounts -p -B -s 2 -a ${gtf} -F ${format} -o ${output}/featureCounts_s2.${index} ${line}
    fi

	count=$(expr $count + 1)

done < ${bamfiles}


