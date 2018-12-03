#! /usr/bin/env python

# Change locus tags names to have bin names for transcriptomic mapping from Bio import SeqIO
import os, sys
from Bio import SeqIO

gbfile=sys.argv[1]
genome=os.path.basename(gbfile).replace(".gbk", "")
final_features=[]
for record in SeqIO.parse(gbfile, "genbank"):
    for f in record.features:
        if f.type =="CDS" or f.type=="gene" and "locus_tag" in f.qualifiers:
            locus = f.qualifiers["locus_tag"][0]
            newlocus = genome + "_" + locus
            f.qualifiers['locus_tag'][0] = newlocus
        final_features.append(f)
    record.features=final_features
    with open(genome+"_modified.gbk", "w") as new_gbk:
        SeqIO.write(record, new_gbk, "genbank")
for record in SeqIO.parse(gbfile, "genbank"):
    for f in record.features:
        if f.type =="CDS" or f.type=="gene" and "locus_tag" in f.qualifiers:
            locus = f.qualifiers["locus_tag"][0]
            newlocus = genome + "_" + locus
            f.qualifiers['locus_tag'][0] = newlocus
        final_features.append(f)
    record.features=final_features       
    with open(genome+"_modified_nuc.fna", "w") as new_fasta:
        SeqIO.write(record, new_gbk, "fasta")
