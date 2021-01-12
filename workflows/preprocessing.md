## Preprocessing Metatranscriptomes and Functional Annotations

This tutorial for preprocessing your MAGs to use in TbasCO assumes that you have assembled your MAGs, quality-checked them, and have (optionally) confirmed the taxonomic assignment of each MAG. 

Once you have assembled and curated your representative, non-redundant MAG database from your community, you will need to: 

1. Perform functional annotations with KofamKOALA
2. Competitively map metatranscriptomic reads to your MAGs

### Dependent Tools: 

- Prokka
- KofamKOALA
- Kallisto
- R packages: tidyverse, tximport

By default, TbasCO uses the KEGG database for functional comparisons. Annotation with the KEGG database can either be done through the online [GhostKOALA](https://www.kegg.jp/ghostkoala/) portal or through the [KofamKOALA](https://www.genome.jp/tools/kofamkoala/) distribution that uses HMMER for searches against HMM profiles. KofamKOALA can either be run through the online portal or ran locally through the command line.

### Formatting Files

First, you will want to format your fasta files in a way where locus tags for each annotated coding region/protein will match between your functional predictions and the metatranscriptome mapping. One way to do this is to get all the necessary file formats (.faa, .ffn, .gbk, etc.) for each assembled genome using Prokka: 

```
for file in *.fna; do 
    N=$(basename $file .fna);
    prokka --outdir $N --prefix $N --cpus 15 $file; 
done
```

Then, concatenate the collection of .ffn and .faa files to use for competitively mapping your metatranscriptomes to and performing functional annotations, respectively: 

_Coding regions:_

```
for filename in */*.ffn; 
    do GENNAME=`basename ${filename%.ffn}`; 
    sed "s|^>|>${GENNAME}_|" $filename; 
done > all-coding-regions.ffn
```

_Proteins:_

```
for filename in */*.faa; 
    do GENNAME=`basename ${filename%.faa}`; 
    sed "s|^>|>${GENNAME}_|" $filename; 
done > all-proteins.faa
```

What these commands will do is create a concatenated file of all coding regions/proteins for all genomes, and appends the name of the genome to each locus tag. So make sure the names of your original fasta files/genome names are something informative or you know what they are (tied to the taxonomic identity for example). 

### Functional Annotation with KofamKOALA

Then annotate the concatenated protein FASTA file with KofamKOALA (download the executables and the required HMM profiles): 

```
./exec_annotation all-proteins.faa;
    -p profiles/; 
    -k ko_list; 
    -o all-proteins-kofamkoala-annots.txt; 
    --cpu 8;
```

What is nice about this pipeline and HMMs over the BLAST-based method that GhostKOALA annotates datasets with the KEGG database is that the output will give you confidence scores according to the threshold cutoff of each HMM and if the annotation for a certain protein is significant with a `*`, from which you can parse out only the significant hits and non-duplicates (i.e. taking the annotation with the highest confidence score for a particular protein). For our purposes, we want the highest confident annotation for a particular protein. You can parse this out with a series of `awk` commands: 

```
grep '*' all-proteins-kofamkoala-annots.txt | awk '{print $2"\t"$3"\t"$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14}' > kofamkoala-annots-sig-mod.txt
sort -u -k1,1 kofamkoala-annots-sig-mod.txt > kofamkoala-annots-sig-mod-nodups.txt
awk '{print $1"\t"$2}' kofamkoala-annots-sig-mod-nodups.txt > ebpr-kofamkoala-annots-sig-mod-nodups-ko-list.txt
```

This will give you both a tab-delimited file of significant annotations, their scores and functional annotations, and then a list of each locus tag and KEGG identifier to use for downstream purposes, and is exactly the output that GhostKOALA gives, except without really knowing quantitatively why things were annotated the way they were. 

### Competitvely Mapping Metatranscriptome Samples to MAGs

These preprocessing steps with kallisto assume that you have performed any necessary quality filtering of metatranscriptomic reads and/or removed rRNA reads. There are always ways post-mapping to remove rRNA reads/loci, which you can do by removing genes that have abnormally high amounts of reads mapping to them or removing genes that are predicted to be ribosomal RNA by Prokka/barrnap when you performed functional annotations. 

For metatranscriptomic mapping, we used the tool `kallisto` for performing pseudoalignment of reads and quantification to coding regions. This software was originally [developed for mapping single-cell RNA-seq data](https://www.nature.com/articles/nbt.3519), but has also been demonstrated to work well with microbial community metatranscriptomics data. Using the concatenated coding regions file, you will first want to index your reference database: 

```
kallisto index -i mag-database-index all-coding-regions.ffn
```

The command for mapping for each sample is formatted as such: 

```
kallisto quant -i index -o outdir fastqfiles
```

You can also create a tab-delimited metadata file, where the first two columns are the names of the R1 and the R2 fastq files, and the 3rd column is the sample name. Such as: 

```
# SRX4072504.qced.R1.fastq	SRX4072504.qced.R2.fastq	sample1
```

Then you can pass that metadata file to a bash script that queues the mapping for each sample, so you can complete the mapping in one step for many samples. For example the bash script looks like: 

```
#! /usr/bin/bash
# Queue metatranscriptomic mapping with kallisto to reference MAG index with a metadata file specifying sample name to filenames of R1 and R2 files

index="$1"
MetadataFile="$2"

echo "Running kallisto for $MetadataFile !"

while read -r a b c; do
	echo $c $a $b;
	kallisto quant -i "$index" -o "$c" "$a" "$b";
done < "$MetadataFile"

echo "Finished running kallisto for $MetadataFile !"
```

And the command to run the script is `bash queue-kallisto.sh $index $metadataFile`. Where you input the name of your index file you created and the metadata file formatted with the fastq filenames and the sample name. This will create a directory for each sample, and in each of those directories are the kallisto mapping results that you will then parse. 

### Processing Kallisto Output 

The following steps can be performed using R. You will need an R package called `tximport` for parsing the kallisto output files. Either place the script within the directory where the output folders are (in each directory is a .tsv and .h5 file for every sample) or if you ran the jobs on a server, `scp` back all the resulting directories and lead your R script to that path. You can find the original documentation for `tximport` for parsing `kallisto` output files and other mapping programs [here](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto).

```
library(tximport)

dir <- "/path/to/your/transcriptomes/"
samples <- read.table(file.path(dir, "samples.txt"), header=TRUE)
files <- file.path(dir, samples$experiment, "abundance.h5")
# number of samples you have
names(files) <- paste0("sample", 1:6)
txi.kallisto <- tximport(files, type="kallisto", txOut = TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var="ID")

```

To use in TbasCO, you can either use raw counts or normalized TPM counts. TbasCO performs its' own normalization calculation to control for the abundance/differing transcript counts mapping back to each genome. 

### Merging Count Table with Functional Annotations
