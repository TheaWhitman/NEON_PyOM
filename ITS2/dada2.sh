#!/bin/bash

# Get the fastq sequences tar from Gluster
# This should be a folder called ITS2, which has all PEAR merged reads in it
# The whole thing is zipped as a .tar.gz
cp /mnt/gluster/twhitman/Cornell_16S/Cornell16S.tar.gz ./

# unzip the .fastq sequencing files
tar -xzf Cornell16S.tar.gz
# Should yield a folder called Cornell16S with all the fwd and rev reads in it
# in the form of .fastq.gz

# Get rid of the .tar.gz file so it's not transferred back automatically 
#rm *.tar.gz

# Get the R installation tar from Gluster
cp /mnt/gluster/twhitman/R/R.tar.gz ./
# May not work, so brought R file over as well.

# untar R installation and remove .tar.gz file
tar -xzf R.tar.gz
rm R.tar.gz

# Get rid of the 16S .tar.gz file so it's not transferred back automatically 
rm *.tar.gz

# Make sure script will use R installation
export PATH=$(pwd)/R/bin:$PATH

# run R script for dada2 analysis
R CMD BATCH dada2.R

# Zip up the relevant files
mkdir Dada2_Results
mv OTUtab.rds Dada2_Results/
mv OTUtab.nochim.rds Dada2_Results/
mv DADA2_seqs.fasta Dada2_Results/
mv DADA2_seqs_nochim.fasta Dada2_Results/
mv track.rds Dada2_Results/
mv pErrR.rds Dada2_Results/
mv pErrF.rds Dada2_Results/

tar -czvf Dada2_Results.tar.gz Dada2_Results/

mv Dada2_Results.tar.gz Dada2_Results_Cornell16S_Pooled_Full.tar.gz

# Move it out to Gluster
cp Dada2_Results_Cornell16S_Pooled_Full.tar.gz /mnt/gluster/twhitman

# Remove remaining files
rm -r Cornell16S
rm -r Dada2_Results
rm *.gz
rm *.fastq
