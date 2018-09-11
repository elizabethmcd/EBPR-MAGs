#! /bin/bash

# Run QC of metagenomic reads on Gluster

# Read in variables for running mapping
meta=$1
metabase=$(basename $1)

mkdir metagenomes

# Programs
tar -xvzf BBMap_38.07.tar.gz

# Copy over metagenomic read files from Gluster and decompress
cp $1 metagenomes/
cd metagenomes
tar -xzf $metabase
cd ..
metarun=$(basename $metabase .fastq.tar.gz)

# Filter the reads
bbmap/bbduk.sh in=metagenomes/$metarun.fastq out=metagenomes/$metarun.qced.fastq qtrim=r trimq=10 maq=10

# Move back only the sorted BAM files to Gluster for binning purposes
cp metagenomes/$metarun.qced.fastq /mnt/gluster/emcdaniel/EBPR-Metagenomes/