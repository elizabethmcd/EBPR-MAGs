# Extracting Metagenome Assembled Genomes from EBPR Enrichment Reactor Time-Series 

This series of workflows demonstrate how to extract, refine, and utilize metagenome assembled genomes (MAGs) from a metagenomic time-series of EBPR reactors. Most of the workflow was constructed to be run on a high-throughput computing system, specifically [CHTC at UW-Madison](http://chtc.cs.wisc.edu/). These steps describe the filtering, assembly, refinement, classification, and annotation steps of using metagenome assembled genomes, and using them for a metatranscriptomics analysis. 

## Dependencies 

- [BBTools Suite](https://sourceforge.net/projects/bbmap/) 
- [Samtools](http://www.htslib.org/download/) 
- [MetaBAT](https://bitbucket.org/berkeleylab/metabat) 
- [CheckM](https://github.com/Ecogenomics/CheckM/wiki)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [ANI Calculator](http://enve-omics.ce.gatech.edu/ani/) 
- [SPAdes](http://cab.spbu.ru/files/release3.12.0/manual.html)
- [Anvi'o](http://merenlab.org/software/anvio/)
- [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM)
- [GTDBK-tk](https://github.com/Ecogenomics/GTDBTk)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [RaxML](https://sco.h-its.org/exelixis/software.html)
- [Prokka](https://github.com/tseemann/prokka)
- [GhostKOALA](https://www.kegg.jp/ghostkoala/)
- [KofamKOALA](https://www.genome.jp/tools/kofamkoala/)
- [antiSMASH](https://antismash.secondarymetabolites.org/#!/start)
- [kallisto](https://pachterlab.github.io/kallisto/)
- [PHANOTATE](https://github.com/deprekate/PHANOTATE)
- [multiPhATE](https://github.com/carolzhou/multiPhATE)
- [OPERA-MS](https://github.com/CSB5/OPERA-MS)
- [Canu](https://canu.readthedocs.io/en/latest/)
- [Racon](https://github.com/lbcb-sci/racon)
- [Quiver](https://github.com/PacificBiosciences/GenomicConsensus)
- [Pilon](http://software.broadinstitute.org/software/pilon/)

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

Then run the `spades-assembly.sub` submission script, which will run the assembly on a high memory node. Even on a high memory node, this should take some time. Additionally, with this dataset we might have to test with metaSPAdes, although other EBPR people have gotten decent results co-assembling with normal SPAdes. **Note:** This takes about a week to run on a high memory node at CHTC. 

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

From the directory where your raw bins are, move those into the corresponding directories: 

```
for file in *.fna; do
    name="${file%.fna*}";
    mv $file $name/$file; 
done
```

Or instead of doing the below, the submission script `Anvio-Metagenomics-Workflow.sub` uses the Docker installation of Anvi'o with CHTC to create contigs/profile databases for each bin for further refinement. Each tarball should be a bin archive, as was created above, with the `.fna` for the bin, and all the `*-EBPR.tar.gz` mapping files with the sorted/indexed BAM files from the previous step. This is all zipped together to fit on Gluster and be transferred as easily as possible. **Before doing any of this**, make sure that the Anvi'o installation versions that you have for the Docker version/VM conda installation match, or else you won't be able to visualize your profile/contig DBs. From where all the bin archives are: 

```
ls $PWD/*.tar.gz > ~/EBPR-BINS-FOLDERS.txt
```

The submission script will then queue by the specific bin, and bring back the contigs/profile DBs for visualizing in Anvi'o. All of those steps are described below if you wish to do so manually.

Concatentate each bin's coverage statistics file into one to get the coverage of that bin through all timepoints. That output file will be used to visualize coverage through time. The BAI files are used for refining each individual bin by having the information of mapped reads from each timepoint to the bin. Note, when dealing with multiple bins and large BAM files, the rest of these steps are probably best done individually folder by folder to 1) Not overwhelm the VM and 2) Keep track of where you are in the refinement process. But keep the coverage files for the end of the analysis.  

```
for file in */*.coverage.txt; do
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
for file in */*.sorted.bam; do anvi-profile -i $file -c CONTIGDB -M 500 --num-threads 14; done
```

Merge the profile databases: 

```
# Example
anvi-merge */*/PROFILE.db -o 3300009517-bin1-SAMPLES-MERGED -c 3300009517-bin.1-contigs.db
```

Note, Anvi'o doesn't like symbols other than '_' so replace names of samples/directories without '.' changed to underscores or no space. Since we are working one at a time with the single bins themselves and not the full metagenomic collection, do not import the binning results from MetaBAT. To load the interactive interface when working on a remote server, you need to login a specific way by running through an SSH tunnel. SSH like so: `ssh -L 8080:localhost:8080 meren@server.university.edu`. From where your samples merged profile database/contigs database is, run `anvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db` for whatever your samples/contigs databases are named. Then on your local browser (preferably Chrome, and that's probably your only option), open `localhost:8080`. Thus, you have a workaround from going through the full metagenomic workflow and then doing `anvi-refine` on the select bins, although you are running the `profile` and `merge` steps multiple times to manually check the bins one by one. 

When working with in the interactive interface, contigs will be hierarchically clustered, and it will be pretty easy to tell when the binner messed up and put contigs together than shouldn't be there based upon differentiall coverage profile of all the contigs. If you have an abberant contig that has a much different differential coverage profile than the other contigs in the bin, you know it needs to be removed from the bin. Once you are satisfied with the contigs you have selected in your manually curated bin, save the collection. I usually just save it as "default". Then summarize the collection with `anvi-summarize -p MERGED_PROFILE/PROFILE.db -c contigs.db -C CONCOCT -o MERGED_SUMMARY` or in our case `anvi-summarize -p MERGED_PROFILE/PROFILE.db -c contigs.db -C default -o refined_bin_name`. In the summary, there will be a `bin-name_contigs.fa` file, which will give the refined bin's contigs in FASTA format. There are a bunch of other statistics in the folder, but the most important is saving the manually refined contigs FASTA file somewhere. So now you have a manually refined bin that you can say you checked for uniform differential coverage, and not just go off of CheckM estimates. 

### Classification and Phylogenetic Relationships 

Bins were classfied with the full classify workflow through EcoGenomics' tool GTDB-tk. This was run on the WEI servers because this wouldn't install/run correctly on the VMs and a docker situation in this instance isn't ideal because of the large database needed. From the output of GTDB-tk, the putative classifications are in the `gtdbtk.bac120.classification_pplacer.tsv` file and the multi-sequence alignment FASTA file for tree making is `gtdbtk.bac120.user_msa.fasta`. Use FastTree for a quick look, but make sure to use RaxML for the final tree since this is a tree of multiple single copy markers. 

### Identifying Accumulibacter Bins

It's usually pretty clear from relative abundance measures in this case which bins are Accumulibacter. To specifically identify which clade they belong to, either perform ANI comparisons against the reference Accumulibacter genomes, as clade genomes by ANI are not similar by more than 85%. Or BLAST the _ppk1_ reference set of genes against the bins, as this is how clades are originally identified by. 

### Functional Annotation 

There are several functional annotation programs/databases I will use for these sets of bins, all for different sets of purposes. Prokka is probably the easiest way to get functional annotation and all file formats for submitting to public repositories. To get KEGG modules for interests in metabolic pathways, annotations are made manually by submitting protein sequences to [https://www.kegg.jp/ghostkoala/](). To get information about biosynthetic gene clusters (BGCs), I will use [antiSMASH](https://antismash.secondarymetabolites.org/) which is available both as an online portal for submitting jobs and through the command line, installed through conda/Docker etc. 

_Prokka:_

To run Prokka, install with `conda install -c conda-forge -c bioconda prokka`. These are all bacterial genomes, so run Prokka by: 

```
for file in *.fna; do 
    N=$(basename $file .fna);
    prokka --outdir $N --prefix $N --cpus 15 $file; 
done
```

To annotate with the KEGG database using GhostKHOALA, concatentate all proteins for all bins into one file like so: 

```
for filename in */*.faa; 
    do GENNAME=`basename ${filename%.faa}`; 
    sed "s|^>|>${GENNAME}_|" $filename; 
done > all-ebpr-prots.faa
```

_GhostKOALA and KofamKOALA:_

GhostKHOALA can take FAA files up to 300 MB in size for annotation at once, and the 58 bins' concatenated proteins file is much smaller than this. Therefore, the concatenated protein file has the genome bin the protein came from and the Prokka annotation in the header for KEGG assignments. Additionally, this concatenated set of proteins can be used to annotate with the KEGG database using KofamKOALA, which is an HMM-based way of performing annnotations. Download the pipeline and required profiles and ko_list from [here](https://www.genome.jp/tools/kofamkoala/), and run the program such as: 

```
./exec_annotation ~/EBPR/prots/all-ebpr-prots.faa;
    -p profiles/; 
    -k ko_list; 
    -o ~/EBPR/prots/all-ebpr-prots-kofamkoala-annots.txt; 
    --cpu 8;
```

What is nice about this pipeline and HMMs over the BLAST-based method that GhostKOALA annotates datasets with the KEGG database is that the output will give you confidence scores according to the threshold cutoff of each HMM and if the annotation for a certain protein is significant with a "*", from which you can parse out only the significant hits and non-duplicates (i.e. taking the annotation with the highest confidence score for a particular protein) with:

```
grep '*' all-ebpr-prots-kofamkoala-annots.txt | awk '{print $2"\t"$3"\t"$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14}' > ebpr-kofamkoala-annots-sig-mod.txt
sort -u -k1,1 ebpr-kofamkoala-annots-sig-mod.txt > ebpr-kofamkoala-annots-sig-mod-nodups.txt
awk '{print $1"\t"$2}' ebpr-kofamkoala-annots-sig-mod-nodups.txt > ebpr-kofamkoala-annots-sig-mod-nodups-ko-list.txt
```

This will give you both a tab-delimited file of significant annotations, their scores and functional annotations, and then a list of each locus tag and KEGG identifier to use for downstream purposes, and is exactly the output that GhostKOALA gives, except without really knowing quantitatively why things were annotated the way they were. 

_antiSMASH:_

To use antismash, use the generated genbank files from the Prokka output. Activate antismash with `activate source antismash` to open the virtual environment installation. A default run of antismash is `antismash <gbk file>`. I also want to run the additional parameters `--cluster-blast` and `--smcogs` to compare the idetified clusters against the antismash database, and also look for orthologous groups of the secondary metabolite clusters. For each predicted cluster, the output is a gbk file with the proteins for each part of the cluster. This will need to be fed through a genbank file parser to put into the final dataframe for annotations. When feeding the files into antismash, it will complain about the contig headers being too long and will change those. But for the purposes of giving the loci range and annotation/locus tag, this shouldn't matter. Copy all the genbank files for all bins to a central directory to run antismash as follows: 

```
for file in *.gbk; do
    name=$(basename $file .gbk);
    antismash $file --clusterblast --smcogs --knownclusterblast --outputfolder $name-antismash;
done
```

Then for each genome, download the antismash output folder, and open the `index.html` file in a browser to interactively view the results for a first pass. 

### Incorporating Metatranscriptomic Datasets 

The metatranscriptomics workflow approach I took goes through these steps: 

1. Manual removal of rRNA
2. Psuedoalign and quantify reads with kallisto
3. Parse genbank files to combine annotation sources

#### Manual removal of rRNA

1. Map reads with 80% identity to the SILVA database
2. Sort the bam file with samtools
3. Create FASTQ files from the sorted bam files with bedtools
4. Create a list of read IDs that were identified as rRNA
5. Filter out the rRNA reads with `filterbyname.sh`, part of the BBTools suite

These steps aren't super important, neither is the overall removal of rRNA, as this can be removed after mapping using the predicted ORFs/annotations. Additionally the index mapped to with kallisto could have just predicted ORFs/non-rRNA as predicted by barrnap (as part of Prokka), so you wouldn't have to waste time manually sorting out rRNA or using cumbersome software like `sortmeRNA`. 

#### Pseudoalignment and Quantification with Kallisto

First create a concatenated index of all predicted ORFs from all genomes with `kallisto index -i ebpr-orfs fastafile`. The headers for each predicted ORF should be the same as what will be in the genbank files to parse from for matching annotation purposes. A good format is `genomeID_orfID`. Then map/quantify with `kallisto quant -i index -o outdir fastqfiles`. This tool is very simple to use, and there are great libraries for parsing the outputs. 

#### Parse Annotations and Output Files

To create a master file of annotations (Prokka for example), use the `genbank-annotation-parser.py`. The R script `ebpr-transcriptomes-kallisto-analyze.R` provides steps for combining mapping files with KEGG annotations, and parsing the results kallisto mapping files to create a master raw counts table.

### Analyzing Phage Genomes

We know that phage has an impact in these bioreactors, but haven't really studied them in detail, due to mostly not having high quality phage genomes. Using the OPERA-MS hybrid assembly of the PacBio + Illumina samples sequenced for 2013-05-23, and then binned with MetaBat, I was able to pull out bins that were identified as "root" with CheckM, and are on either 5 contigs or less (2 of them are on a single and 2 contigs, respectively) and range from 200-300 kb, which is about the expected size for phage genomes. There are a couple different methods to identify and annotate phage sequences and predict their host relationships bioinformatically. 

- virSorter to predict phage sequences
- PHANNOTATE to identify phage ORFs
- multiPhAte to call functional annotations on the ORFs
- VRCA to compare TNF of phage/genomes to associate with host

Can then use the assembled genomes, annotations, and predicted relationships with more resolved sequenced samples to draw better connections

## Incorporating Long PacBio Reads

All workflows described above entail assembling short read Illlumina shotgun data into consensus bins, which are usually fragmented into lots of contigs and don't include conserved/repetative regions such as ribosomal subunits (sometimes, some of my HQ ones do). This can be alleviated by using long read sequencing technlogies such as PacBio/Oxford Nanopore. We previously sequenced one of the samples with Illumina shotgun sequencing also with Pacbio sequencing, specifically the May 28th 2013 sample. I previously tried scaffolding some of the bins recovered from this sample specifically after depreplication with the PacBio reads with very little success, and the tools for doing this are weird and poorly documented. Here, I attempt two approaches: long read assembly only  (with some Illumina polishing for error correction), hybrid assembly, and post-assembly correction and polishing. Binning will be performed by taking the differential coverage from the all three May 2013 samples, the 13th, 23rd, and 28th.  

### Hybrid Assembly of Illumina short reads and PacBio long reads

The PacBio RS II reads have a high error rate ~20%, compared to the new HiFi reads with low error rate, and should be corrected with the Illumina data before any assembly. Error correction and polishing can be done on the long reads with and without the Illumina data after the assembly as well. 

The software OPERA-MS was recently released for hybrid assembly of short + long reads for metagenomes, and is suppposedly better than IDBA-hybrid and SPAdes, and I dont' think metaSPAdes has hybrid assembly options yet. First, interleaved files have to be split up into R1 and R2 read files, which can be done with BBtools `reformat.sh`, with usage: 

```
Usage:  reformat.sh in=<file> out1=<outfile> out2=<outfile2>
```

In combination with the long reads, OPERA-MS can be run as: 

```
perl OPERA-MS.pl; 
    --short-read1 2013-05-23-EBPR.R1.fastq; 
    --short-read2 2013-05-23-EBPR.R2.fastq; 
    --long-read 2013-05-23-PACBIO-qced.fastq; 
    --out-dir RESULTS
```

A test run on one of the VMs with > 16 GB RAM finished this job in less than 24 hours with the following assembly statistics: 

```
[Thu Aug  1 16:07:55 2019]	Assembly stats
Number of contigs: 255720
Assembly size: 334183580 bp
Max contig size: 603481 bp
Contig(s) longer than 1Mbp: 0
Contig(s) longer than 500kbp: 1
Contig(s) longer than 100kbp: 218
Contig N50: 3445 bp
```

The hybrid assembly can then be used for binning purposes, repeating the mapping and binning steps as described above. There are a couple of ways this could be approached. Sometimes people will only have long reads for a couple of samples, use those for the assembly, and then have many more samples only with short reads. Then they will map all the samples with short reads to the hybrid long-read approach assembly to get differential coverage, and bin based off of that. As a preliminary test I only mapped the corresponding illumina reads to the assembly, but could try all the samples to see how binning improves. 

### Long Read Assembly and Polishing with Short Reads

The established methods of genome/metagenome assembly have been: short read assembly only, long read assembly only, or hybrid assembly of the two. Something I've learned recently is sort of a combination of the latter two, where the long reads are assembled, and the bases are polished or error corrected with the corresponding short read data to create a high-quality assembly. This involves several steps, mostly with intermediate polishing. 

First assemble the long Pacbio data with `canu`. Install with: 

```
git clone https://github.com/marbl/canu.git
cd canu/src
make -j <number of threads>
```

The basic canu run setup is `canu -p <assembly-prefix> -d <assembly-directory> genomeSize=<number> pacbio-raw <fastq/fasta file>`. There are some suggested modificiations for assembling a metagenome, where the genomeSize should be set to around 5Mb, and playing with the error rates and increasing memory. Don't run canu on "quality filtered" reads, the first step of canu will specifically correct/quality filter the PacBio reads/bases specific to the technology, whereas other QC programs are specific to short read programs. 

```
corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200
```

After the assembly, there will need to be some polishing on the raw, uncorrected raw long read assembly. For PacBio data, Quiver is suggested, which is produced specifically by PacBio. Racon has also been mentioned, and can support polishing with Illumina data. Additionally the Canu documentation recommends if you have Illumina data that you can polish with Pilon. Some resources suggest doing both polishing approaches, such as polishing the Canu long read assembly only with the long reads first with Racon, and then polishing with short read data using Pilon, such as this [ONT article](https://github.com/nanoporetech/ont-assembly-polish). This docker container is specific for nanopore reads, and I haven't yet found a similar pipeline for PacBio + Illumina reads assembly/polishing. The general steps are: 

1. Canu assembly
2. Polish long-read contigs with Racon
3. Map Illumina reads onto assembly
4. Polish contigs with short read data using Pilon

*Important notes*: Important note #1: There is some discrepancy (in emails and previous lab notes) in whether the pacbio data came from the May 23rd or the 28th sample date, for which we have both illumina and pacbio data for, so probably want to map/polish with both just to be certain. Although in the 5/28/2019 IMG collection for this sample, all the Pacbio and Illumina stuff is there, but not everybody is certain. Important note #2: It's not entirely clear from the IMG site was is Pacbio raw reads/assemblies vs Illumina reads/assemblies for this sample, and it's confusing when downloading things. Make sure when downloading raw Illumina vs Pacbio reads, the pacbio ones will say pbio in the name, whereas the Illumina raw reads file is really large. Additionally the assembly contains multiple pacbio runs (???) together, and isn't the raw data. The raw data is split into 5 individual fasta files, which canu can assemble together with wildcard `*.fasta`.  

Example command: 

```
/home/emcdaniel/Ext-Inst/canu/Linux-amd64/bin/./canu -p EBPRpac -d EBPR_long genomeSize=5m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=50 -pacbio-raw *.fasta
```

Putting the `canu` commands in my bin interfered with the native perl I had installed probably through anaconda, so had to run from within the directory. So be it. Top part of the intial report: 

```
-- In sequence store './EBPRpac.seqStore':
--   Found 1743135 reads.
--   Found 3853006791 bases (770.6 times coverage).
--
--    G=3853006791                       sum of  ||               length     num
--    NG         length     index       lengths  ||                range    seqs
--    ----- ------------ --------- ------------  ||  ------------------- -------
--    00010         3675     77331    385303969  ||       1000-1361       146635|---------------------
--    00020         2978    195845    770602413  ||       1362-1723       305202|-------------------------------------------
--    00030         2657    333382   1155904482  ||       1724-2085       449283|---------------------------------------------------------------
--    00040         2435    485109   1541202810  ||       2086-2447       366954|----------------------------------------------------
--    00050         2256    649617   1926505405  ||       2448-2809       217259|-------------------------------
--    00060         2098    826741   2311804375  ||       2810-3171       111946|----------------
--    00070         1946   1017374   2697104798  ||       3172-3533        55227|--------
--    00080         1781   1224025   3082406235  ||       3534-3895        28802|-----
--    00090         1565   1453605   3467706422  ||       3896-4257        17196|---
--    00100         1000   1743134   3853006791  ||       4258-4619        11206|--
--    001.000x             1743135   3853006791  ||       4620-4981         7924|--
--                                               ||       4982-5343         5871|-
--                                               ||       5344-5705         4394|-
```

Showing that most sequences fall within the length of 1300-2500, which seems to be decent. 