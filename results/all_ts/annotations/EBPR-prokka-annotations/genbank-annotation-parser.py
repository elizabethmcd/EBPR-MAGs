#!/usr/bin/env python

import argparse, os
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

# Arguments 

parser = argparse.ArgumentParser(description = "Parse Prokka annotations from a genbank file to get a tab-delimited file of annotations and gene sizes")
parser.add_argument('gbk_file', metavar='GBK', help='Annotation file from Prokka in Genbank format')
parser.add_argument('--annotation', default='gene_annot.txt', help="Output: Functional annotation for all genes of a genome .gbk file (Default: gene_annot.txt)")

args = parser.parse_args()

# Input and output files
GBK = args.gbk_file
OUT_ANNO = open(args.annotation, "w")

# Headers
# OUT_ANNO.write("gene_id\tProkka_annotation\tsize_bp\taccession\n")

# Parse the GBK file
for record in SeqIO.parse(GBK, "genbank"):
    genome_id = os.path.basename(GBK).replace(".gbk", "").strip().splitlines()[0]
    for f in record.features:
        if f.type == 'CDS':
            locus, span, function = (f.qualifiers["locus_tag"][0], f.location, f.qualifiers["product"][0])
            beg = span.start # biopython already reads in the genbank file to start counting from 0, so don't have to subtract 1
            end = span.end 
            strand = str(span.strand)
            r = strand.replace("-1", "r")
            direction = r.replace("1", "f")
            size= (end-beg)
            try:
                gene_acc = f.qualifiers["gene"][0]
            except KeyError:
                gene_acc = ""
            gene_id = genome_id + "_" + locus
            # Write out annotation information
            OUT_ANNO.write('%s\t%s\t%d\t%s\n' % (gene_id, function, size, gene_acc))