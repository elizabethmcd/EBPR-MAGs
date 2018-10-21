# Extracting Metagenome Assembled Genomes from EBPR Enrichment Reactor Time-Series 

This series of workflows demonstrate how to extract, refine, and utilize metagenome assembled genomes (MAGs) from a metagenomic time-series of EBPR reactors. Most of the workflow was constructed to be run on a high-throughput computing system, specifically [CHTC at UW-Madison](http://chtc.cs.wisc.edu/). 

## Dependencies 

- [BBTools Suite](https://sourceforge.net/projects/bbmap/) 
- [Samtools](http://www.htslib.org/download/) 
- [MetaBAT](https://bitbucket.org/berkeleylab/metabat) 
- CheckM
- [Prodigal](https://github.com/hyattpd/Prodigal)
- ANI Calculator 
- [SPAdes](http://cab.spbu.ru/files/release3.12.0/manual.html)
- Anvi'o
- GTDBK-tk
- Prokka
- GhostKOALA

### Filter Raw Metagenomic Sequences 

To begin with, the raw metagenomic sequences need to be quality filtered using `bbduk`. This workflow uses the transfer protocol system Gluster, therefore all of the input files and resulting large fastq/BAM/SAM files will be deposited into Gluster as well. Make a list of the metagenomic reads to be filtered from Gluster, and deposit the text file into your home folder of CHTC to run the submission script from. From the directory the metagenomic reads are in, run `ls $PWD/*.fastq.tar.gz > ~/metaToQC.txt` to create the list to queue from for filtering. Run the submission script `EBPR-QC.sub` which will execute `EBPR.QC.sh`. 

### Perform Mapping of Metagenomic Reads to Assemblies 

Now that the metagenomic reads are filtered, map the reads to the assembled contigs to get differential coverage for binning. For individual assemblies, each individual metagenomic timepoint is mapped to every assembled timepoint. To create the lists for the list of references, or the assembled contigs, metagenomes, and then mapping combinations of refs > metagenomes, youw will do the following: 

```
# From the directories in which your metagenomic reads/references are: 

ls $PWD/*.qced.fastq > ~/metagenomeList.txt

ls $PWD/*.a.fna > ~/refList.txt

```

Then run the script `makeMappingCombos.py` from the submit node, otherwise don't create a submission script. This simple python script just makes combination of every ref > metagenome pair for mapping. This will creat the file mappingCombos.txt, from which the submission job `EBPR-Mapping.sub` will queue from. 

Before running the job, you will need to make an installation of Samtools using an interactive job. You will need to download [Samtools](http://www.htslib.org/download/) and place it in your home directory. Use the submission script `Samtools-CHTC-install.sub` to run an interactive job using `condor_submit -i Samtools-CHTC-install.sub`. Once the interactive job has started, do the following: 

```
tar xvfj samtools-1.9.1.tar.bz2
cd samtools-1.9.1
make
make prefix=../samtools install
cd ..
ls samtools
tar czvf samtools.tar.gz samtools/
ls
exit

#Move samtools to the zipped/ folder
mv samtools.tar.gz zipped/samtools.tar.gz

```

You will now have a new zipped folder of samtools that has been installed on CHTC, and when unzipped during your job submission, will work appropriately. You will only need to make this once, and not everytime you need to use samtools. Unless you remove it or CHTC updates something. Then run the job with `condor_submit EBPR-Mapping.sub`. These scripts use a percent identity sequence cutoff of 95%, which is not stringent enough to split into strains, supposedly. This will create consensus MAGs after binning. This will queue from mappingCombos.txt, and put the sorted BAM files onto Gluster. From there you will bin by differential coverage with MetaBAT. 

### Binning with MetaBAT

We will bin using MetaBAT and will use the handy dandy docker container for MetaBAT. The MetaBAT pipeline is run as such: 

```
runMetaBat.sh <options> assembly.fasta sample1.bam sample2.bam [....]
```

Where for each assembly of contigs, we are calculating the depths of reads mapped from each sample. The arguments thus need to be changed for every assembly, and then all the samples mapped to that assembly. For CHTC to queue from a list of file paths, I believe you have to manually enter the arguments for transfer/running the binning process. Additionally MetaBAT might need python 2.7, therefore python will have to be manually installed on CHTC using an interactive submission. First download python 2.7 from [here](https://www.python.org/downloads/release/python-2715/) and move the tarball to your home directory of CHTC. After submitting the submission script with `condor_submit -i Python-CHTC-install.sub` do the following: 

```
mkdir python
tar -xzf Python-2.7.15.tgz
cd Python-2.7.15
./configure --prefix=$(pwd)/../python
make
make install
cd ..
ls python
ls python/bin
export PATH=$(pwd)/python/bin:$PATH
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py
pip install numpy
pip install pandas
tar -czvf python.tar.gz python/

```
We can also use the Docker version of MetaBAT, so that should work without this, but now we have the python installation in case pandas is needed in the future. 

The script `makeBinningCombos.py` will manually make a text file to queue from, but you can also use bash to call assembly ID's from the BAM files in order to bin the corresponding samples based on their differential coverage of mapped reads to the particular assembly. Either way will work, and having a nice python 2.7 installation will be handy at some point. To use the pipeline, pull it down from Github with `git clone https://github.com/sstevens2/checkm-chtc-pipeline.git` into a folder on your home directory of CHTC. 

### Quality Check 

There will be multiple quality check steps throughout this process, but it's good before any refinement or bin comparisons to have a high-level view of how good these bins are across assemblies/samples. We will use CheckM to assess the completeness and redundancy of the bins, and also get a sense of classifications based on marker genes.

To run CheckM, we will use Sarah Steven's Docker image and submission scripts found [here](https://github.com/sstevens2/checkm-chtc-pipeline). Move all the bins to a references directory, but rename them before with the `makeANIcombos.sh` script, which will give each bin the name of the assembly timepoint it came from. Additionally, can also just do it this way for ANI comparisons by dumping all of the prodigal predicted genomes into one directory for all-v-all. However, then it will do comparisons within the same group, as in the same cmparisons of bins in a timepoint, so it will take a bit longer.   

### Inspecting Best Bins from Multiple Samples using ANI Comparisons 

To compare similar bins among the different binning efforts, we will use ANI comparisons to group them into similar bins, and then based off of completion/redundancy estimates, can start to pick the "best" ones. For each of the 11 time points, there are approximately 35-50 bins. We need to do pairwise ANI comparisons for all of the bins to find the similar ones. First rename all of the bins in each of the bin directories with the `makeANIcombos.sh` script. This will append the assembly name the specific bins came from for doing ANI comparisons. This will also need to be pairwise within a sample to know if there are similar organisms within the same time point. 

To schedule mass ANI comparisons on CHTC, we will use Sarah Stevens' [DAG pipeline](https://github.com/sstevens2/ani_compare_dag) for this. The important part is that ANI comparisons can only be run on the nucleotide coding regions and not tRNA/rRNA, so we have to perform gene calling with Prodigal to only get the coding sequences.

First make a file listing the directories of bins with the files ending in `.fna`, such as: 

```
/home/emcdaniel/EBPR-Bins-NUCS/3300009517
/home/emcdaniel/EBPR-Bins-NUCS/3300026282
/home/emcdaniel/EBPR-Bins-NUCS/3300026283
/home/emcdaniel/EBPR-Bins-NUCS/3300026284
/home/emcdaniel/EBPR-Bins-NUCS/3300026286
/home/emcdaniel/EBPR-Bins-NUCS/3300026287
/home/emcdaniel/EBPR-Bins-NUCS/3300026288
/home/emcdaniel/EBPR-Bins-NUCS/3300026289
/home/emcdaniel/EBPR-Bins-NUCS/3300026299
/home/emcdaniel/EBPR-Bins-NUCS/3300026302
/home/emcdaniel/EBPR-Bins-NUCS/3300026303
```

Then use Prodigal to create fna files for each bin only containing coding sequences, so each timepoint can be compared pairwise against all other bins only for the coding regions using the ANI comparison tool. Get Prodigal from [here](https://github.com/hyattpd/Prodigal). To install, do the following with an interactive job submission: 

```
git clone https://github.com/hyattpd/Prodigal.git
tar -cvf Prodigal.tar.gz Prodigal/
condor_submit -i Prodigal-CHTC-install.sub
# in interactive session
tar -xvf Prodigal.tar.gz
mkdir prodigal
cd Prodigal
make install INSTALLDIR=../prodigal/
cd ..
tar -cvf prodigal.tar.gz prodigal/
exit
```
Once you have the coding sequences for the bins, pairwise ANI comparisons can be run for all bins extracted from all time points. 

First get the ANI calculator and unpack to your home folder of CHTC and clone the ANI DAG: 

```
wget https://ani.jgi-psf.org/download_files/ANIcalculator_v1.tgz
tar -xzvf ANIcalculator_v1.tgz
git clone https://github.com/sstevens2/ani_compare_dag.git
```

Then decompress the prodigal fna ran files and move the directories into the ANI DAG folder. For each bin directory, make a file listing the bins for each time point in this manner: 

```
for a in $(awk '{print $1}' list-of-files.txt); 
    do for b in $(awk '{print $2}' list-of-files.txt); 
        do echo $a-vs-$b | xargs -n1 | sort -u | xargs; 
    done; 
done > allcombos.txt
```

Then manually remove the duplicate lines of same-vs-same, because the ANI comparisons won't like performing those and it's a waste to run ANI comparisons within the same timepoint. The runAll by default can perform this as well. Make directories of each combination by: 

```
while read line;
    do mkdir $line;
done < allcombos.txt

```

While it could be done this way, all the bins can just be dumped into one directory for all-v-all comparisons. That way less work is done for splitting things into different directories. This will also give comparisons for within a timepoint, but then we know if there are multiple "similar" organisms in a single timepoint. 

### Co-Assembly with SPAdes

So far we have only worked with single assemblies from each timepoint. This will work well for abundant organisms, especially those with a lot of strain variation (i.e. Accumulibacter). However since we are trying to extract as many high quality bins as possible, having a co-assembly will help these efforts. These single assemblies were done with SPAdes, and not the metaSPAdes version, done at the JGI. We will try to co-assemble with SPAdes as well. This will need to be done on CHTC on one of the high memory nodes. To install: 

```
wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz
tar -xzf SPAdes-3.12.0-Linux.tar.gz
```

The binaries will be in `SPAdes-3.12.0-Linux/bin/` to run the assembly and transfer the input executable when running it. 

To run the assembly, we need to use one of the high memory nodes at CHTC, which will change the submission script slightly. The input file to run the assembly is a YAML file: 

```

    [
      
      {
        orientation: "fr",  
        type: "paired-end",
        interlaced reads: [
            "2005-06-14-EBPR.qced.fastq","2007-09-17-EBPR.qced.fastq","2008-04-24-EBPR.qced.fastq","2009-02-02-EBPR.qced.fastq","2010-07-15-EBPR.qced.fastq","2010-10-04-EBPR.qced.fastq","2011-01-10-EBPR.qced.fastq","2012-01-23-EBPR.qced.fastq","2013-05-13-EBPR.qced.fastq","2013-05-23-EBPR.qced.fastq", "2013-05-28-EPBR.qced.fastq"
        ]
      }
    ]

```

Then run the `spades-assembly.sub` submission script, which will run the assembly on a high memory node. Even on a high memory node, this should take some time. Additionally, with this dataset we might have to test with metaSPAdes, although other EBPR people have gotten decent results co-assembling with normal SPAdes. 

### Check Assembly Quality 

To compare assemblies, use the program Quast for a preliminary look at the assembly statistics. The ultimate measure of assembly quality is number of reads mapping back to the assembly, but this is good for comparing different assembly methods at first. The usage for checking metagenomic assemblies is `./metaquast <contig-assemblies-to-report-on>`. Note that part of this tool uses BLAST, and can take quite a long time for two metagenomic assembly comparisons/statistics. I've only run this with 1 thread and very little memory, and it took over 2 hours to get past the BlastN step. Additionally, it calls references from the NCBI database, and then maps the contigs to the databases for assembly quality. Which seems to be a waste of time for just checking over and over again the quality of the co-assembly with different parameters, since I'm having to _de novo_ assemble everything anyways.  

Additionally the N50 stats can probably just be calculated from the `statswrapper.sh` in the BBTool suite. So that's a better way to get a quick look at the N50 statistics 

### Parsing ANI Comparisons and CheckM Estimates Across Time Points 

In order to pick identical bins between timepoints, ANI comparisons and CheckM estimates will be leveraged. From the ANI and CheckM analyses, there should be two output files: the all.cleaned ANI file and the concatenated `lineage.txt` files for all CheckM runs of the bins. Using the Rscript `EBPR-bins-ANI-vs-checkm.R` by `Rscript EBPR-bins-ANI-vs-checkm.R anifile checkmfile` and you should get an output of identical bins by the 99% ANI threshold and greater than .50 alignment fraction in one direction, and then can pick bins that way. 

### N50 Stats

In addition to picking bins by ANI/completion, in cases where it isn't clear which bin is best based off of ANI values, CheckM estimates, or genome size, use the N50 estimates. To do this, you need the `statswrapper.sh` script from the BBTools suite. Make a list of all bins, and make sure they are coded with their assembly timepoint from which they originated. Run the following: 

```
list=list-of-bins.txt
while read "file"; do /Volumes/mcmahonlab/home/emcdaniel/Software/bbmap/statswrapper.sh in="$file"; done <"$list" > all-bins-stats.txt
```

This will create a tab-delimited file of stats for each bin, and can then be added to the master file for comparing identical bins and picking the "best" one. With the N50 stats, the Rscript can be run as follows: `Rscript scripts/EBPR-bins-ANI-vs-checkm.R results/EBPR-BINS.all.ani.out.cleaned results/EBPR-bins-checkm-results.txt results/all-bins-stats.txt`. With the order of the arguments being `anifile checkmfile statsfile`. Use these combinations to pcik the "best" bin in a high matching set. 

Additionally for further quality control, the Rscript `checkm-get-less-50-complete.R` will output a liss of the .fna files/bins with less than 50% completeness to completely toss - which should've just been done before the ANI comparisons to lessen the amount of jobs and data to parse through. In a directory of all the bins together, move to `crap-bins` directory by: 

```
craplist=less-than-50-complete-bins.txt
while read "file"; do mv $file ../crap-bins/; done <"$craplist" 
```

Note: This can probably much easier be done with a bash script to identify less than 50% complete from the file and not make an intermediate file like this. Whatever, I always realize I could've used bash after I've done the thing.
 
Bins will have to be manually picked by identifying identical bins and choosing the best one based off of size, completeness, redundancy, and N50/L50 scores. Then the unique bins can be chosen by what is missing in the "identical" dataset. To be certain that you have dereplicated bins, re-run the ANI/CheckM steps on that subset and check the output. 

### Quality Check of Selected Bins and Bins from Co-Assembly 

Re-run CheckM/ANI/N50stats and the corresponding R script to get the final statistics on the selected bins to make sure there are no duplicates and they are of high enough quality going forward. The CheckM/ANI/N50 steps are repeated until you are left with a non-redundant set of high-quality bins that you think comprise the community. The ultimate check if a set of bins truly represents the community is mapping back reads from all timepoints to the set of selected bins.

### Map Metagenomic Reads to all Bins 

Now we will map all the metagenomic reads to the extracted bins to make sure we have a representative set of the whole community and each time point since most of these bins weren't extracted from a coassembly binning strategy. Since CHTC only has 2/3 servers that can pull from Gluster AND have enough disk space (~500GB) to perform mapping by timepoint, every meta to reference mapping comparison will have to be split up into individual jobs. To set this up, we will run the `makeMappingCombos.py` script on a list of references and metagenomes like when mapping the metagenomes to the assemblies. This will also make the file output name much easier since the script does that as well. **Make sure you have enough disk/directory quota space on Gluster when transferring files back**. For example, 50 bins with 10 metagenomic time points is 500 jobs/tarballed directories to move back to Gluster with the coverage stats and sorted/indexed BAM files for that time point mapped to the specific bin. 

From the directories in which the bins and QCed metagenomes reside: 

```
ls $PWD/*.fna > ~/refList.txt
ls $PWD/*.qced.fastq > ~/metagenomeList.txt
```

From the home directory, run `makeMappingCombos.py`. Submit the submission script `mapMetasToRefs.sub`. This script will queue the job by timepoint vs bin, saving in each directory that timepoint-vs-bin's coverage statistics and the sorted/indexed BAM files for refining with Anvi'o. 

### Manually Refine Bins with Anvi'o 

In order to make these bins as high of quality as possible, we will manually refine each bin in the de-duplicated set to make sure everything looks okay before scaffolding with PacBio reads. In order to do so, I will be roughly following the Anvi'o [Metagenomic Workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). Usually what Meren's lab recommends is performing a coassembly, mapping the reads back to the coassembly, and getting bins. This way you would input all sorted/indexed BAM files along with the contigs, and then also provide the bin collection saying which contigs ended up in which bins. Since I mapped each timepoint back to the individual assemblies, and de-replicated bins based on the "best" one in a particular timepoint, this scenario is a little different. I could take from the output of MetaBAT and only include the contigs that went into the bins I selected, but I forsee a couple of problems with the interactive interface since I had to manually pick bins. Instead, I will do it the opposite way, where instead of creating a contigs/profile database for the entire metagenome, and then using `anvi-refine` for each individual bin, I believe I can create contigs/profile databases for each individual bin. Use the sorted/indexed BAM files output from the `mapMetasToRefs.sub` workflow. 

Install Anvi'o on a VM/server with sudo access with `conda install -c bioconda -c conda-forge anvio=5.1.0 diamond bwa` if you use Anaconda, like I usually do (and make sure you are using python version 3.6+). Transfer all of the tarballed `ref-vs-meta` folder with the coverage/indexed BAM files from Gluster. Create a folder for each bin with: 

```
ls *-EBPR.tar.gz | awk -F- '{print $1"-"$2}' | uniq >> binsList.txt
for dir in $(cat binsList.txt); do mkdir "$dir"; done
for file in *-EBPR.tar.gz; do
    name="${file%-vs*}";
    mv $file $name/$file; 
done
```

Concatentate each bin's coverage statistics file into one to get the coverage of that bin through all timepoints. That output file will be used to visualize coverage through time. The BAI files are used for refining each individual bin by having the information of mapped reads from each timepoint to the bin. Note, when dealing with multiple bins and large BAM files, the rest of these steps are probably best done individually folder by folder to 1) Not overwhelm the VM and 2) Keep track of where you are in the refinement process. But keep the coverage files for the end of the analysis.  

```
for file in */*/*.coverage.txt; do
    name="${file%-vs*}"; 
    cat -- "$file" >> "${name}".coverage.txt; 
done
```

Rename all bins files of contigs with the file extension `.fa`. **Important Note**: Only reformat fasta files if you did so before mapping, because then the contig names for your contig database for the particular bin won't match the contig names within the mapping BAM file. For the most part, the names of the contigs and formatting of FASTA files might be ok without reformatting for these purposes. If not, reformat the fasta BIN files, and remap the metagenomic reads to each bin.

```
for file in */*.fa; do
    name="${file%.fa*}";
    anvi-script-reformat-fasta $file -o $name-fixed.fa -l 0 --simplify-names;
done
``` 

Create a contigs database for each bin. 

```
anvi-gen-contigs-database -f $BIN.fna -o $name-contigs.db
```

Run HMMs stats on each of the contigs databases to check for single copy core genes:

```
anvi-run-hmms -c $file --num-threads 15
```

We don't need to run functional annotation or import external functional annotations because we are refining the contigs based off of coverage from the metagenomic timepoints for each genome, and will work with functional annotations after the genomes are manually refined. Gene-level taxonomy from Kaiju or Centrifuge can also be imported, and might help with the manual refinement process. For now I will skip that and just base the refinement process off of checking coverage of contigs. 

Profiling the BAM files: 

```
for file in */*.sorted.bam.bai; do 
    name="${file%-vs*}"; 
    base=$(basename $file); 
    echo $base >> "${name}".samples.txt; 
done
```

This will create a file in each bin folder of the mapped reads to the specific bin for each timepoint to profile the BAM files for that bin separately from. Make sure these BAM files are already sorted/indexed with samtools, or do it within Anvi'o. The minimum number of bases on the contigs will need to be changed due to splits issues, since we are looking at individual genome bins, and not the whole metagenome. 

```
for file in *.qced.bam; do anvi-profile -i $file -c CONTIGDB -M 500; done
```

### Scaffolding with Long Reads  

- Have to find the accession numbers for raw reads
- 2018-10-03: PacBio filtered reads have been found, need to download and look at methods for hybrid assembly of short reads with long reads, or scaffolding the bins and not complete reassembly 

### Preliminary Classification

- DAG to go through classifying contigs > bins, compared to GTDB-tk output

### Classification and Phylogenetic Relationships 

- GTDB, make docker image and test-run on condor
- compare to classifications made through JGI pipeline

### Identifying Accumulibacter Bins

- Check by ANI
- Check by BLAST of _ppk1_ genes

### Functional Annotation 

- Prokka easy way, for having easy access to all files for a genome bin
- GhostKOALA for KEGG modules, metabolic pathways, using in TcT manuscript 
- 2018-10-09 Note: I had zero clue that GhostKOALA was a web server and you have to manually annotate a metagenome/genome bin ONE AT A TIME. Dagnabbit.  

### Incorporating Metatranscriptomic Datasets 

#### Filter Metatranscriptomic Read

#### rRNA Depletion 

#### Competitively map Reads to Bins 

#### Count reads and Normalize

### TbasCO Incorporation of Metatranscriptomic Reads and KEGG Annotations 

### Metabolic Pathway Prediction 

### Putative Interactions

### Depositing Genome Drafts to Public Repository 

- Open Science Framework for initial sharing with collaborators