#!/bin/bash

#SBATCH --job-name=createBamList
#SBATCH --time=6-23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/bamList_%A_%a.out
#SBATCH -e ./o/bamList_%A_%a.err

### cleaning the rawbam folder and getting Bamfiles file

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data
input_dir=/mnt/beegfs/training/CITIWorkshops/RNASeq/data

ls ${input_dir}/rawbam/*Aligned.sortedByCoord.out.bam > ${base}/Bamfiles.txt
