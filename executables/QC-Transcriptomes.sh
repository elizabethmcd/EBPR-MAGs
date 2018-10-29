#! /bin/bash 

# setup files
mkdir transcriptomes
cp $1 transcriptomes/
cd transcriptomes/
tar -xzvf *.gz 
file=$(basename $1 .tar.gz)
name=$(basename $1 .fastq.tar.gz)
cd ..

# software setup
tar -xzvf sortmerna-2.1-linux-64-multithread.tar.gz
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# QC 
./fastp -i transcriptomes/$file -o transcriptomes/$name.qced.fastq

# Index databases 
cd sortmerna-2.1b

./indexdb_rna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db

#Run the sorting program
./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db --reads ../transcriptomes/$name.qced.fastq  --fastx --aligned ../$name-rRNA --other ../$name-nonrRNA --log -v -m 1 -a 1

# copy back to gluster
cd ..
tar -czvf $name-rRNA.tar.gz $name-rRNA/
tar -czvf $name-nonrRNA.tar.gz $name-nonrRNA/
cp $name-rRNA.tar.gz $name-nonrRNa.tar.gz /mnt/gluster/emcdaniel/


