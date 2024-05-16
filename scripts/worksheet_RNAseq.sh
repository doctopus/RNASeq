
ssh <user>@10.88.8.31
pwd:

###Setting up your space on Beegfs and copy over data and scripts 
cd /mnt/beegfs/beherat2/
mkdir RNAseq
cd RNAseq

###Copy scripts to this directory: Scripts are run from the scripts directory
rsync -r --progress /mnt/beegfs/beherat2/kinaseRNASeq/scripts .
cd scripts
mkdir o
chmod 775 o

mkdir data
mkdir resources

# Build Star Genome Reference Index
sbatch buildSTARIndex.sh
#sbatch buildSTARIndex.sh #chmod 775 on resources where the Log.out file resides
#TODO looks like it still not able to move the log file so revert the resources chmod to 755

#Run Initial QC
#TODO In current code fastqFiles.list.txt needs to have an extra empty line in the end else last sample is not processed
#Array size needs to be +1 the lines of fastqFiles.list.txt
sbatch --array=1-25 FastQC.sh
cd ../data/FastQC
module load python/3.9.5
multiqc .

#Only if Trimming not needed then run alignment (Index is taken from training, based on old gtf nd fa)
cd ../scripts
sbatch --array=1-25 STAR_align.sh

#Check MultiQC File: If Trimming needed then
#Run Trimmomatic
sbatch runTrimmomatic.sh

#Create a trimmedFastq.list.txt for the next QC run
sbatch createTrimmedFastqList.sh
#This creates a trimmedFastq.list.txt in scripts
#TODO For now added a trailing empty line, but fix this code so no empty lines are needed
#TODO creates paths as /home/beherat2.. not wrt /mnt/
# Looks like the default line break is alright. No need to add a extra line

#Run QC Again
sbatch --array=1-25 runFastqcOnTrimmed.sh
#Output is in data/FastQC_Trimmed
#While in FASTQC_Trimmed directory
cd ../data/FastQC_Trimmed
module load python/3.9.5
multiqc .

# Check If Trimmed MultiQC is OK 


#Run STAR alignment (STARIndex is made from latest gtf v45 and fa v14 files)
cd ../scripts
sbatch --array=1-25 STARIndex_align.sh


## Tapas - Need to run createBamFilesList.sh and source the bamfiles from /data in the next step
sbatch --array=1 createBamFilesList.sh
##Tapas Check if bamfiles.list.txt has a extra line at the end else add one.


sbatch --array=1-25 runFeatureCounts.sh

## run downstream analysis
sbatch runDownstream.sh