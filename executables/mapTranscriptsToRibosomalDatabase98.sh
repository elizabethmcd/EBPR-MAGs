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
tar -xzvf $r1base
tar -xzvf $r2base
cd ..
cp /mnt/gluster/emcdaniel/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta .

# Perform mapping  
bbmap/bbmap.sh ref=SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta in1=transcriptomes/$r1file in2=transcriptomes/$r2file outm=mappingResults/$tranname.bam -Xmx48g minid=0.98

# Move back to gluster
mkdir $tranname-ribosomal-mapped
mv mappingResults/*.bam $tranname-ribosomal-mapped
tar -czvf $tranname-ribosomal-mapped.tar.gz $tranname-ribosomal-mapped
cp $tranname-ribosomal-mapped.tar.gz /mnt/gluster/emcdaniel

# cleanup 
rm *.fasta
rm *.tar.gz
rm */*.tar.gz
rm */*.fastq


