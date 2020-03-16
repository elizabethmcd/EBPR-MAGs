#! /usr/bin/python

import sys

def usage():
    print("Count number of bases in any fasta file")

if __name__ == "__main__":
    if len(sys.argv) !=3:
        usage()
        exit()

    fastafile = open(sys.argv[1], "rU")
    fasta=fastafile.readlines()
    count = 0
    for line in fasta:
        if not line.startswith(">"):
            count = count + (len(line) - 1)
    outname = sys.argv[2]
    with open(outname, "wt") as out:
        out.write(sys.argv[1] + "\t" + str(count) + "\n")