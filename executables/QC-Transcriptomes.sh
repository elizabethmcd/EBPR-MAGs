#! /bin/bash 

# Metatranscriptome variable
meta=$1
mkdir trans

# Programs
tar -xvzf BBMap_38.07.tar.gz

# setup
cp $1 trans/
cd trans
tar -xzf $meta
cd ..
metabase=$(basename $1 .tar.gz)
metarun=$(basename $metabase .fastq.tar.gz)

# QC 
bbmap/bbduk.sh in=trans/$metabase out=trans/$metarun.bbduk.qced.fastq qtrim=r trimq=10 maq=10

# copy back to gluster
tar -czvf $metarun.bbduk.qced.fastq.tar.gz trans/$metarun.bbduk.qced.fastq
cp $metarun.bbduk.qced.fastq.tar.gz /mnt/gluster/emcdaniel/


