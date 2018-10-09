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
    outsort="${filename%.*}".sorted.bam;
    ./samtools/bin/samtools sort $file -o mappingResults/$outsort;
done

# Get depth 
for file in mappingResults/*.sorted.bam; do
    outdepth="${filename%.*}".depth;
    ./samtools/bin/samtools depth $file > mappingResults/$outdepth;
done 

# Sorted, indexed BAM file
for file in mappingResults/*.sorted.bam; do
    outindex="${filename%.*}".sorted.indexed.bam;
    ./samtools/bin/samtools index $file > mappingResults/$outindex;
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

# Bring back statistics and the sorted/indexed BAM file, put into one directory, and zip to Gluster
mkdir $refname-vs-$metaname
mv *.coverage.txt $refname-vs-$metaname/
mv mappingResults/*.sorted.index.bam $refname-vs-$metaname/
tar -cvzf $refname-vs-$metaname.tar.gz $refname-vs-$metaname/
cp $refname-vs-$metaname.tar.gz /mnt/gluster/emcdaniel/

