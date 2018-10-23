#! /bin/bash 

# Perform anvio metagenomic workflow to get contig/profile DBs for sets of bins to refine

# copy over bin directory, unpack everything 
file=$1
mkdir wd
cp $file wd/
cd wd
tar -xzvf $file
bin=$(basename $file .tar.gz)
binName="${bin//./}"
cd $bin
for file in *.tar.gz; do
    tar -xzvf $file; 
done

# contigs database
anvi-gen-contigs-database -f $bin.fna -o $bin-contigs.db 
anvi-run-hmms -c $bin-contigs.db 

# profile (sorted) BAM files 
for file in */*.sorted.bam; do
    anvi-profile -i $file -c $bin-contigs.db -M 500; 
done 

# merge profile databases 
anvi-merge */*/PROFILE.db -o $binName-SAMPLES-MERGED -c $bin-contigs.db

# collect and transfer back results 
mkdir $bin-anvi-results 
mv $binName-SAMPLES-MERGED/ $bin-contigs.db $bin-anvi-results/
tar -czvf $bin-anvi-results.tar.gz $bin-anvi-results/
cp $bin-anvi-results.tar.gz /mnt/gluster/emcdaniel
