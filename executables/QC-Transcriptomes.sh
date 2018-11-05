#! /bin/bash 

# setup files
cp $1 .
tar -xzvf *.gz 
file=$(basename $1 .tar.gz)
name=$(basename $1 .fastq.tar.gz)

# software setup
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# QC 
./fastp -i $file -o $name.qced.fastq --cut_by_quality3

# copy back to gluster
tar -czvf $name.qced.fastq.tar.gz $name.qced.fastq
cp $name.qced.fastq.tar.gz /mnt/gluster/emcdaniel/


