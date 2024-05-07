###login on beegfs
#Open terminal
ssh usr@10.88.8.31
pwd:

###Setting up your space on Beegfs and copy over data and scripts 
cd /mnt/beegfs/{user}/
mkdir RNAseq
cd RNAseq

###Copy scripts to this directory:
rsync -r --progress /mnt/beegfs/training/CITIWorkshops/RNASeq/scripts .
mkdir output

cd scripts
mkdir o
chmod 775 o
#Check Octal file permission by $ stat -c "%n %a" {file/folder}
##Aray size should be the total number of rows in fastqFiles.list.txt file ?or n+1
sbatch --array=1-25 runFastQC.sh
cd ../output/FastQC
module load python/3.9.5
multiqc .

cd ../scripts
sbatch --array=1-25 runSTARalign.sh

#Run create bamFiles.list; check if they have an extra line in the end
sbatch --array=1-4 createBamFilesList.sh


sbatch --array=1-4 featureCounts.sh

## run downstream analysis
sbatch runDownstream.sh