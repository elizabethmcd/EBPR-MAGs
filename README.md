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

The script `makeBinningCombos.py` will manually make a text file to queue from, but you can also use bash to call assembly ID's from the BAM files in order to bin the corresponding samples based on their differential coverage of mapped reads to the particular assembly. Either way will work, and having a nice python 2.7 installation will be handy at some point. 

## To Write:

### Quality Check 

### Inspecting Best Bins from Multiple Samples using ANI Comparisons 

### Quality Check 

### Map Metagenomic Reads to all Bins 

### Manually Refine Bins with Anvi'o 

### Reassemble Bins with Long Reads 

