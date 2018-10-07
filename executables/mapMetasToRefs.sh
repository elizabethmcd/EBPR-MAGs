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

# Copy over metagenomic reads and refined bin set
ref=$1
cp $1 refs/
refbase=$(basename $1 .fna)
cp /mnt/gluster/emcdaniel/EBPR-Metagenomes/*.qced.fastq metagenomes/.

# Perform mapping for every reference 
for file in metagenomes/*.qced.fastq; do
    metabase=$(basename $file .qced.fastq);
    bbmap/bbmap.sh ref=refs/$ref in=metagenomes/$file outm=mappingResults/$refbase-vs-$metabase.bam idtag minid=0.95 nodisk -Xmx48g;
done 

# Sorted BAM files
for file in mappingResults/*.bam; do
    ./samtools/bin/samtools sort $file -o mappingResults/${file%.bam}.sorted.bam

# Get depth 
for file in $resultsDir/*.sorted.bam; do
    outname="${filename%.*}".depth;
    ./samtools/bin/samtools depth $file > mappingResults/$outname;
done 

# Sorted, indexed BAM file
for file in mappingResults/*.sorted.bam; do
    outname="${filename%.*}".sorted.indexed.bam;
    ./samtools/bin/samtools index $file > mappingResults/$outname;
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
for file in mappingResults/*.depth; do
    python calc-mapping-stats.py $file; 
done

# Bring back statistics, don't need the mapped results when just want depth/relative abundance/ANI/PID metrics
# If need mapped reads in directories with bins, see the RefiningBins scripts/workflow
mkdir $refbase-results
mv coverage.txt $refbase-results/
mv mappingResults/*.sorted.index.bam $refbase-results/
cp -r $refbase-results/ /mnt/gluster/emcdaniel/.

# Cleanup
rm *.py
rm *.tar.gz
rm -rf metagenomes/ refs/ $resultsDir/ $refbase-results/