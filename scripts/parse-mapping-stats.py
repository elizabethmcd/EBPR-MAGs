#! /usr/bin/python 

import os, re, sys
import pandas as pd

def usage():
    print("Calculate coverage from depthfile")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        usage()
        exit()

    # Arguments
    depthfile = sys.argv[3]
    covFileHeader = "filename\tref\tmeta\tmetareads\tbaseCovSum\tcoveredBases\ttotalBases\tAvgCov\n"
    regCovLine = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
    ref, meta = os.path.basename(depthfile).split("-vs-")
    meta = meta.split(".")[0]

    # Genome size
    refLengths = open(sys.argv[1], "r")
    refLength = int()
    for line in refLengths:
        if re.search(ref, line):
            refLength = int(line.split()[1])
    
    # Metagenomic reads
    metaReadsfile = open(sys.argv[2], "r")
    metareads = int()
    for line in metaReadsfile:
        if meta in line:
            metareads = int(line.split(" ")[1].rstrip())
            break
    metaReadsfile.close()
    
    # File empty if no reads recruited, otherwise not empty
    outname = sys.argv[4]
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