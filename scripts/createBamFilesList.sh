#!/bin/bash

#SBATCH --job-name=createBamList
#SBATCH --time=6-23:00:00
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/bamList_%A_%a.out
#SBATCH -e ./o/bamList_%A_%a.err

### cleaning the rawbam folder and getting Bamfiles file

wdir=$(realpath "$(pwd)")
pdir=$(dirname "$wdir")
base=${pdir}/data
scripts=${wdir}
input_dir=${base}/BAMIndex

##ls ${input_dir}/rawbam/*Aligned.sortedByCoord.out.bam > ${base}/Bamfiles.txt
ls "${input_dir}"/*Aligned.sortedByCoord.out.bam > "${scripts}/bamFiles.list.txt"
# Add an extra line in the end
#echo >> "${scripts}/bamFiles.list.txt"