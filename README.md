# Extracting Metagenome Assembled Genomes from EBPR Enrichment Reactor Time-Series 

These series of workflows demonstrate how to extract, refine, and utilize metagenome assembled genomes (MAGs) from a metagenomic time-series of EBPR reactors. Most of the workflow was constructed to be run on a high-throughput computing system, specifically [CHTC at UW-Madison](http://chtc.cs.wisc.edu/). 

## Dependencies 

- [BBMap](https://sourceforge.net/projects/bbmap/) 
- [Samtools](http://www.htslib.org/download/) 
- [MetaBAT](https://bitbucket.org/berkeleylab/metabat) 
- ANI Calculator 
- Anvi'o

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



