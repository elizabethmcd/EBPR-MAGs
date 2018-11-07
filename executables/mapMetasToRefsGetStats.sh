#! /bin/bash

# Run mapping of short metagenomic reads to reference genomes, calculate statistics
# Queues by metagenomic timepoint sample

# Setup 
mkdir metagenomes
mkdir refs
mkdir mappingResults 

# Programs 
tar -xvzf BBMap_38.07.tar.gz
tar -xvzf samtools.tar.gz
tar -xvzf python.tar.gz

# Python
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
chmod u+x *.py

# Copy over metagenomic timepoint and bin
ref=$1
meta=$2
outname=$3
refbase=$(basename $ref)
metabase=$(basename $meta)
refname=$(basename $ref .fna)
metaname=$(basename $meta .qced.fastq)
cp $1 refs/
cp $2 metagenomes/

# Perform mapping  
bbmap/bbmap.sh ref=refs/$refbase in=metagenomes/$metabase outm=$outname idtag minid=0.95 nodisk -Xmx48g

# Sorted BAM files
for file in mappingResults/*.qced.bam; do
    outsort=$(basename $file .qced.bam).sorted.bam;
    ./samtools/bin/samtools sort $file -o mappingResults/$outsort;
done

# Get depth 
for file in mappingResults/*.sorted.bam; do
    outdepth=$(basename $file .sorted.bam).depth;
    ./samtools/bin/samtools depth $file > mappingResults/$outdepth;
done 

# Reference lengths file 
for file in refs/*.fna; do
    python countBases.py $file; 
done
cat refs/*.len > refGenomes.len

# Metagenomic reads file 
for file in metagenomes/*.fastq; do
    awk '{s++}END{print FILENAME,s/4}' $file >> metaReads.txt;
done

# Create stats file 
for file in mappingResults/*.depth; do
    python calc-mapping-stats.py $file; 
done

# Bring back statistics and the BAM files with sorted/indexed BAM file, put into one directory, and zip to Gluster
cp *.coverage.txt /mnt/gluster/emcdaniel/

rm *.fastq
rm *.bam 
