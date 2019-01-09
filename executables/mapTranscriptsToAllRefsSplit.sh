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
transcriptome=$1
tranbase=$(basename $transcriptome .tar.gz)
tranname=$(basename $transcriptome .fastq.bbduk.qced-nonrRNA.fastq.tar.gz)
cp $1 transcriptomes/
cd transcriptomes
tar -xzvf *.gz
cd ..
cp /mnt/gluster/emcdaniel/EBPR-Bins-NUCS/*.fna .

# Perform mapping  
bbmap/bbsplit.sh ref=3300009517-bin.12.fna, 3300009517-bin.13.fna, 3300009517-bin.1.fna,3300009517-bin.29.fna, 3300009517-bin.30.fna, 3300009517-bin.31.fna,3300009517-bin.3.fna, 3300009517-bin.42.fna, 3300009517-bin.47.fna,3300009517-bin.52.fna, 3300009517-bin.6.fna, 3300009517-bin.7.fna,3300026282-bin.4.fna, 3300026282-bin.5.fna, 3300026283-bin.18.fna,3300026283-bin.19.fna, 3300026283-bin.21.fna, 3300026283-bin.28.fna,3300026284-bin.6.fna, 3300026284-bin.9.fna, 3300026287-bin.17.fna,3300026287-bin.29.fna, 3300026287-bin.38.fna, 3300026287-bin.4.fna,3300026288-bin.15.fna, 3300026288-bin.19.fna, 3300026288-bin.23.fna,3300026288-bin.30.fna, 3300026288-bin.32.fna, 3300026288-bin.34.fna,3300026288-bin.43.fna, 3300026288-bin.6.fna, 3300026289-bin.23.fna,3300026289-bin.24.fna,3300026289-bin.28.fna, 3300026289-bin.38.fna,3300026289-bin.41.fna, 3300026299-bin.12.fna, 3300026299-bin.22.fna,3300026299-bin.26.fna, 3300026299-bin.49.fna, 3300026302-bin.10.fna,3300026302-bin.20.fna, 3300026302-bin.24.fna, 3300026302-bin.25.fna,3300026302-bin.27.fna,3300026302-bin.31.fna, 3300026302-bin.32.fna,3300026302-bin.46.fna, 3300026302-bin.47.fna,3300026302-bin.59.fna,3300026302-bin.61.fna, 3300026302-bin.62.fna, 3300026303-bin.32.fna,3300026303-bin.38.fna, 3300026303-bin.42.fna, 3300026303-bin.46.fna in=transcriptomes/$tranbase basename=mappingResults/$tranname-_out%.sam minratio=0.56 minhits=1 maxindel=16000

# Move back to gluster
mkdir $tranname-mapped
mv mappingResults/*.sam $tranname-mapped
tar -czvf $tranname-mapped.tar.gz $tranname-mapped
cp $tranname-mapped.tar.gz /mnt/gluster/emcdaniel

# cleanup 
rm *.fna
rm *.tar.gz
rm */*.tar.gz
rm */*.fastq


