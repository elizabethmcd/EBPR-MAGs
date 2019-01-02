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
tranname=$(basename $transcriptome .fastq.bbduk.qced.fastq.tar.gz)
cp $1 transcriptomes/
cd transcriptomes
tar -xzvf *.gz
cd ..
cp /mnt/gluster/emcdaniel/EBPR-Bins-NUCS/*.fna .

# Perform mapping  
bbmap/bbsplit.sh ref=3300026282-bin.4,3300026284-bin.9,3300026288-bin.43 in=transcriptomes/$tranbase basename=mappingResults/out_%.sam minratio=0.56 minhits=1 maxindel=16000

# Move back to gluster
mkdir $tranname
mv mappingResults/*.sam $tranname
tar -czvf $tranname.tar.gz $tranname
cp $tranname.tar.gz /mnt/gluster/emcdaniel

# cleanup 
rm *.fna
rm *.tar.gz
rm */*.tar.gz
rm */*.fastq


