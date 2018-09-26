#! /bin/bash

# Transfer over qced fastq files from Gluster 
cp /mnt/gluster/emcdaniel/EBPR-Metagenomes/*.qced.fastq .
mkdir assembly

# SPAdes setup 
tar -xvf SPAdes-3.12.0-Linux.tar.gz
# export PATH=$(pwd)/SPAdes-3.12.0/bin:$PATH

# Run SPAdes assembly 
SPAdes-3.12.0-Linux/bin/./spades.py -t 16 -m 500 -k 21,33,55,77,99,127 --dataset EBPR-reads-to-assemble-test.yaml -o assembly/

# Move back to Gluster
cp -r assembly/ /mnt/gluster/emcdaniel/

rm *.fastq

