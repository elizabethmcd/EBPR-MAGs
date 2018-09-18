#! /bin/bash

# Transfer over qced fastq files from Gluster 
cp /mnt/gluster/emcdaniel/EBPR-Metagenomes/*.qced.fastq .

# SPAdes setup 
tar -xvf SPAdes-3.12.0-Linux.tar.gz
export PATH=$(pwd)/SPAdes-3.12.0/bin:$PATH

# Run SPAdes assembly 
spades.py -t 12 -m 250 -k 21,33,55,77,99,127 --dataset EBPR-reads-to-assemble.yaml -o assembly/

# Move back to Gluster
cp -r assembly/ /mnt/gluster/emcdaniel/

rm *.fastq

