#! /usr/bin/python 

import os, re, sys
import pandas as pd

def usage():
    print("Calculate coverage from depthfile")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage()
        exit()

    # Arguments
    filename = sys.argv[1]
    covFileHeader = "filename\tref\tmeta\tmetareads\tbaseCovSum\tcoveredBases\ttotalBases\tAvgCov\n"
    regCovLine = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
    ref, meta = os.path.basename(filename).split("-vs-")
    meta = meta.split(".")[0]

    # Genome size
    refLengths = open("refGenomes.len", "r")
    refLength = int()
    for line in refLengths:
        if re.search(ref, line):
            refLength = int(line.split()[1])
    
    # Metagenomic reads
    metaReadsfile = open("metaReads.txt", "r")
    metareads = int()
    for line in metaReadsFile:
        if meta in line:
            metareads = int(line.split(" ")[1].rstrip())
            break
    metaReadsFile.close()
    
    # File empty if no reads recruited, otherwise not empty
    outname = ref + ".coverage.txt"
    if os.stat(filename).st_size ==0:
        if os.path.exists(outname):
            with open(outname, "a") as outfile: 
                outfile.write(regCovLine.format(filename, ref, meta, metareads, "0", "0", refLength, "NA"))
        else: 
            with open(outname, "a") as outfile:
                outfile.write(covFileHeader)
                outfile.write(regCovLine.format(filename, ref, meta, metareads, "0", "0", refLength, "NA"))
    else:
        depth_df = pd.read_table(filename, header=None)
        covCol = len(depth_df.columns) - 1 # last column number
        if os.path.exists(outname):
            with open(outname, "a") as outfile:
                outfile.write(regCovLine.format(filename, ref, meta, metareads, depth_df[covCol].sum(), depth_df[covCol].count(), refLength, depth_df[covCol].sum() / float(refLength)))
        else: 
            with open(outname, "a") as outfile: 
                outfile.write(covFileHeader)
                outfile.write(regCovLine.format(filename, ref, meta, metareads, depth_df[covCol].sum(), depth_df[covCol].count(), refLength, depth_df[covCol].sum() / float(refLength)))