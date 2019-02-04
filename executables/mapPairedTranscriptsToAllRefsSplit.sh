#! /bin/bash

# Run mapping of metatranscriptomes to reference genomes split up using bbsplit, will map against best hit
# Queues by metatranscriptomic sample
# Perform statistics and count reads mapped with HTSeq 

# Setup 
mkdir transcriptomes
mkdir mappingResults 

# Programs 
tar -xvzf BBMap_38.07.tar.gz
tar -xvzf samtools.tar.gz

# Copy over metagenomic timepoint and bin
r1=$1
r2=$2
r1base=$(basename $r1)
r2base=$(basename $r2)
tranname=$3
r1file=$(basename $r1 .tar.gz)
r2file=$(basename $r2 .tar.gz)
r1name=$(basename $r1 .fixed.qced.fastq.tar.gz)
r2name=$(basename $r2 .fixed.qced.fastq.tar.gz)
cp $1 transcriptomes/
cp $2 transcriptomes/
cd transcriptomes
tar -xzvf *.gz
cd ..
cp /mnt/gluster/emcdaniel/EBPR-ORFs/*.ffn .

# Perform mapping  
bbmap/bbsplit.sh ref=3300009517-bin.12.ffn,3300009517-bin.13.ffn,3300009517-bin.1.ffn,3300009517-bin.29.ffn,3300009517-bin.30.ffn,3300009517-bin.31.ffn,3300009517-bin.3.ffn,3300009517-bin.42.ffn,3300009517-bin.47.ffn,3300009517-bin.52.ffn,3300009517-bin.6.ffn,3300009517-bin.7.ffn,3300026282-bin.4.ffn,3300026282-bin.5.ffn,3300026283-bin.18.ffn,3300026283-bin.19.ffn,3300026283-bin.21.ffn,3300026283-bin.28.ffn,3300026284-bin.6.ffn,3300026284-bin.9.ffn,3300026287-bin.17.ffn,3300026287-bin.29.ffn,3300026287-bin.38.ffn,3300026287-bin.4.ffn,3300026288-bin.15.ffn,3300026288-bin.19.ffn,3300026288-bin.23.ffn,3300026288-bin.30.ffn,3300026288-bin.32.ffn,3300026288-bin.34.ffn,3300026288-bin.43.ffn,3300026288-bin.6.ffn,3300026289-bin.23.ffn,3300026289-bin.24.ffn,3300026289-bin.28.ffn,3300026289-bin.38.ffn,3300026289-bin.41.ffn,3300026299-bin.12.ffn,3300026299-bin.22.ffn,3300026299-bin.26.ffn,3300026299-bin.49.ffn,3300026302-bin.10.ffn,3300026302-bin.20.ffn,3300026302-bin.24.ffn,3300026302-bin.25.ffn,3300026302-bin.27.ffn,3300026302-bin.31.ffn,3300026302-bin.32.ffn,3300026302-bin.46.ffn,3300026302-bin.47.ffn,3300026302-bin.59.ffn,3300026302-bin.61.ffn,3300026302-bin.62.ffn,3300026303-bin.32.ffn,3300026303-bin.38.ffn,3300026303-bin.42.ffn,3300026303-bin.46.ffn in1=transcriptomes/$r1file in2=transcriptomes/$r2file basename=mappingResults/_out%.sam

# Move back to gluster
mkdir $tranname-mapped
mv mappingResults/*.sam $tranname-mapped
tar -czvf $tranname-mapped.tar.gz $tranname-mapped
cp $tranname-mapped.tar.gz /mnt/gluster/emcdaniel

# cleanup 
rm *.ffn
rm *.tar.gz
rm */*.tar.gz
rm */*.fastq


