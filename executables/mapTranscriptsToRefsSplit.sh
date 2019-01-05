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
bbmap/bbsplit.sh ref=3300026282-bin.4.fna,3300026284-bin.9.fna,3300026288-bin.43.fna in=transcriptomes/trans/$tranbase basename=mappingResults/out_%.sam minratio=0.56 minhits=1 maxindel=16000

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


