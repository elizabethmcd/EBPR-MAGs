#! /usr/bin/env python 

import pandas as pd
import os, sys

annots = pd.read_table(sys.argv[1], sep="\t", names=['genome_name', 'locus', 'annotation'])
counts = pd.read_table(sys.argv[2], sep="\t", names=['locus', 'count'])
outf = sys.argv[3]

annot_counts = annots.merge(counts[['locus', 'count']], how="left", left_on="locus", right_on="locus")
annot_counts.to_csv(outf+"-all-bin-counts.txt", sep="\t", index=False)