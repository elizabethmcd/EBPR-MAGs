#! /bin/bash

# Binning with MetaBat

# Read in assembly variable, make a binning directory, and move over corresponding 
ref=$1
refbase=$(basename $1 .a.fna)
mkdir $refbase-bins

# Copy over assembly and BAM files from Gluster

cp $1 $refbase-bins/
cp /mnt/gluster/emcdaniel/EBPR-Mapping-Results/$refbase*.sorted.bam $refbase-bins/
cd $refbase-bins/
mkdir $refbase-EBPR-bins

# Get depth 
jgi_summarize_bam_contig_depths --outputDepth $refbase-depth.txt *.bam

# Run metabat
metabat2 -i $refbase.a.fna -a $refbase-depth.txt -o $refbase-EBPR-bins/bin

# Zip up
tar -cvf $refbase-EBPR-bins.tar.gz $refbase-EBPR-bins/
