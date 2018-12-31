#! /bin/bash 

# Reformat fastq files and toss broken reads to not mess up downstream analyses

# File and software setup
meta=$1
mkdir trans
tar -xvzf BBMap_38.07.tar.gz
cp $1 trans/
cd trans
tar -xzf $meta
cd ..
metabase=$(basename $meta .tar.gz)
metarun=$(basename $meta .fastq.tar.gz)

# Reformat reads 
bbmap/reformat.sh in=trans/$metabase out=trans/$metarun.fixed.fastq tossbrokenreads

# Copy back to gluster
tar -czvf $metarun.fixed.fastq.tar.gz trans/$metarun.fixed.fastq
cp $metarun.fixed.fastq.tar.gz /mnt/gluster/emcdaniel