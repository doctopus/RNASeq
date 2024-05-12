#!/bin/bash

#SBATCH --job-name=Downstream
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2 
#SBATCH --mem-per-cpu=3900
#SBATCH -o ./o/downstream_%A_%a.out
#SBATCH -e ./o/downstream_%A_%a.err

echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST

module load python/2.7.18
module load R/4.3.1
module load geos/3.4.3

wdir=$('pwd')
pdir=$(dirname "$wdir")
base=${pdir}/data
#resources=/mnt/beegfs/training/CITIWorkshops/RNASeq/resources
resources=${pdir}/resources

mkdir -p ${base}/outs/counts

hall_file=${resources}/GSEA_gene_sets/h.all.v7.2.symbols.gmt

Rscript ./ScanandMerge.R ${base}

Rscript ./createDESeqObj.R ${base}

Rscript ./downstreamAnalysis.R ${base} ${hall_file} 
