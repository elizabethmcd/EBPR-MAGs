#! /bin/bash

# bin genomes with metabat from hybrid assembly

# Read in assembly variable, make a binning directory, and move over corresponding 
ref=$1
mkdir bins

# Copy over assembly and BAM files from Gluster

cp $1 bins
cp /mnt/gluster/emcdaniel/EBPR-Mapping-Results/2013-05-23-illumina-pacbio-hybrid-assembly-contigs-vs-2013-05-23-EBPR.qced.sorted.bam bins/
cd bins

# Get depth 
jgi_summarize_bam_contig_depths --outputDepth hybrid-assembly-depth.txt *.bam

# Run metabat
metabat2 -i $ref -a hybrid-assembly-depth.txt -o bin

# Zip up
cd ..
tar -cvf hybrid-assembly-bins.tar.gz bins/