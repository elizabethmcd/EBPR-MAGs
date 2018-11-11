#! /bin/bash

# Run mapping of metatranscriptomes to concatenated reference genomes
# Queues by metatranscriptomic sample
# Perform statistics and count reads mapped with HTSeq 

# Setup 
mkdir transcriptomes
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
transcriptome=$1
tranbase=$(basename $transcriptome .tar.gz)
tranname=$(basename $transcriptome _mRNA-nonrRNA.fastq.tar.gz)
cp $1 transcriptomes/
cd transcriptomes
tar -xzvf *.gz
cd ..
cp /mnt/gluster/emcdaniel/all-ebpr-nucs.fna .

# Perform mapping  
bbmap/bbmap.sh ref=all-ebpr-nucs.fna in=transcriptomes/$tranbase outm=$tranname.bam idtag minid=0.95 nodisk -Xmx48g

# Sorted BAM files
for file in mappingResults/*.bam; do
    outsort=$(basename $file .bam).sorted.bam;
    ./samtools/bin/samtools sort $file -o mappingResults/$outsort;
done

# Move back to gluster
mkdir $tranname
mv mappingResults/*.sorted.bam $tranname
tar -czvf $tranname.tar.gz $tranname
cp $tranname.tar.gz /mnt/gluster/emcdaniel

# cleanup 
rm *.fna
rm *.tar.gz
rm */*.tar.gz
rm */*.fastq


