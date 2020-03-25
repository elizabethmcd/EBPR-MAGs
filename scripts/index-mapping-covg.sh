#! /bin/bash

# Append genome name to scaffolds file for each bin, combine into one file
for file in *.fa; do GENOME=`basename ${file%.fa}`; sed -i "s|^>|>${GENOME}~|" $file; done
cat *.fa > all-bins.fasta

# make list of scaffold-to-bins pairings
# for use with InStrain's genome_wide workflow to profile genome wide SNPs and covg
grep '>' bins.fasta | sed 's|[<>,]||g' | awk -F '~' '{print $1"~"$2"\t"$1}' > scaffolds-to-bins.tsv

# index the combined reference file
bowtie2-build bins.fasta bt2/combined-bins.fasta

# queue mapping 
# create file with list of complete paths of each metagenome files, second column the output destination
# If have just a few metagenomes, can write a script of a for loop such as 

for file in ~/EBPR/coassembly/metagenomes/cleaned_fastqs/*.qced.fastq; do
	name=$(basename $file .qced.fastq);
	bowtie2 -p 4 -x bt2/R1R2-EBPR-bins-index.fasta -q $file > mapping/$name-mapping.sam;
done

# stats of % reads mapping back to an assembly or reference within BAM
samtools idxstats your_file.bam
awk 'FNR==NR{sum+=$3;next}{print $1 "\t" $3/sum}' <(samtools idxstats your_file.bam) <(samtools idxstats your_file.bam)

## getting mapping %'s from reads mapping to the assembly from MetaBat depth files
# create counts file
for file in *-depth.txt; 
    do name=$(basename $file -depth.txt); 
    awk '{print $1"\t"$4"\t"$6"\t"$8}' $file > $name-counts.txt; 
done

# sum the relevant columns
awk '{sum +=$4} END {print sum}' 3300026302-counts-raw.txt

## mapping to indvidiual bins stats
# queue mapping with coverM to get relative abundance
coverm genome \
    --reference POS-bins-combined.fasta \
    -s "~" \
    -m relative_abundance \
    --interleaved ../cleaned_fastqs/*-interleaved.fastq \
    --min-read-aligned-percent 0.75 \
    --min-read-percent-identity 0.95 \
    --min-covered-fraction 0 \
    -x fasta -t 5 &> log.txt

# coverage values

coverm genome \
    --reference POS-bins-combined.fasta \
    -s "~" \
    -m covered_fraction \
    --interleaved ../cleaned_fastqs/*-interleaved.fastq \
    --min-read-aligned-percent 0.75 \
    --min-read-percent-identity 0.95 \
    --min-covered-fraction 0 \
    -x fasta -t 5 &> POS_coverage_calcs.txt


# queue inStrain profiling

# parse inStrain TSVs with python dictionary to create combined output table