#! /usr/bin/python

import os, argparse, csv, re


parser = argparse.ArgumentParser(description = "Create combinations of assemblies with multiple BAM samples for calculating coverage and performing binning with MetaBAT")
parser.add_argument('BAM_DIR', metavar="BAM_DIR", help="Directory of sorted BAM files of metagenomic timepoints mapped to assemblies")
parser.add_argument("ASSEM_PATH", metavar="ASSEM_PATH", help="Path to directory of assembly files")
parser.add_argument("--outfile", default="BinningCombos.txt", help="Output file for binning combinations")

args = parser.parse_args()

BAM_DIR = args.BAM_DIR 
ASSEM_PATH = args.ASSEM_PATH
OUTFILE = open(args.outfile, "w")

# Read in list of BAM files, create list of assemblies > corresponding mapped BAM files with paths to the files

files = os.listdir(BAM_DIR)
ref_list = []
meta_list = []
for file in files:
    ref = file.split("-vs-")[0]
    if ref not in ref_list:
        ref_list.append(ref)
    meta_list.append(file)
for ref in ref_list:
    matching = [m for m in meta_list if ref in m]
    matching_with_paths = [BAM_DIR + x for x in matching]
    key = ASSEM_PATH + ref + ".a.fna"
    OUTFILE.write("{0},{1}\n".format(key, ",".join(matching_with_paths)))





    





