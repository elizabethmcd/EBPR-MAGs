#! /bin/bash

# Run mapping of short metagenomic reads to reference genomes, calculate statistics
# Queues by metagenomic timepoint sample

# Metareads
meta=$1
metabase=$(basename $1 .qced.fastq)

# Setup 
mkdir metagenomes
mkdir refs
resultsDir=$metabase-mappingResults 
mkdir $resultsDir

# Programs 
tar -xvzf BBMap_38.07.tar.gz
tar -xvzf samtools.tar.gz
tar -xvzf python.tar.gz

# Python
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
chmod u+x *.py

# Copy over metagenomic reads and refined bin set
cp $1 metagenomes/
cp /mnt/gluster/emcdaniel/all-man-refined-bins/*.fna refs/.

# Perform mapping for every reference 
for file in refs/*.fna; do
    refbase=$(basename $file .fna);
    bbmap/bbmap.sh ref=$file in=metagenomes/$meta outm=$resultsDir/$refbase-vs-$metabase.bam idtag minid=0.95 nodisk -Xmx48g;
done 

# Sorted BAM files
for file in mappingResults/*.bam; do
    ./samtools/bin/samtools sort $file -o $resultsDir/${file%.bam}.sorted.bam

# Get depth 
for file in mappingResults/*.sorted.bam; do
    outname="${filename%.*}".depth;
    samtools depth $file > $resultsDir/$outname;
done 

# Reference lengths file 
for file in refs/*.fna; do
    python countBases.py $file; 
done
cat refs/*.len > refGenomes.len

# Metagenomic reads file 
for file in metagenomes/*.fastq; do
    awk '{s++}END{print FILENAME,s/4}' $filename >> metaReads.fastq.txt;
done

# Create stats file 
for file in $resultsDir/*.depth; do
    python calc-mapping-stats.py $file; 
done

# Bring back statistics, don't need the mapped results when just want depth/relative abundance/ANI/PID metrics
# If need mapped reads in directories with bins, see the RefiningBins scripts/workflow
cp coverage.txt /mnt/gluster/emcdaniel/.

# Cleanup
rm *.py
rm *.tar.gz
rm -rf metagenomes/ refs/ $resultsDir