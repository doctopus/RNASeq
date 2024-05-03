###login on beegfs
#Open terminal
ssh usr@10.88.8.31
pwd:

###Setting up your space on Beegfs and copy over data and scripts 
cd /mnt/beegfs/bangalp2/
mkdir RNAseq
cd RNAseq

###Copy scripts to this directory:
rsync -r --progress /mnt/beegfs/training/CITIWorkshops/RNASeq/scripts .
mkdir data

cd scripts
mkdir o
chmod 775 o
##Aray size should be the total number of rows in fastq.list.txt file ?or n+1
sbatch --array=1-4 FastQC.sh
cd ../data/FastQC
module load python/3.9.5
multiqc .

cd ../scripts
sbatch --array=1-4 STAR_align.sh

sbatch --array=1-4 featureCounts.sh

## run downstream analysis
sbatch runDownstream.sh