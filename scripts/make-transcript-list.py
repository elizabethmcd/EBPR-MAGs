#! /usr/bin/env python

import os

# read in lists
with open("transcription-timepoints.txt", "r") as t:
    transList = t.read().splitlines()

with open("bin-names.txt", "r") as b:
    binList = b.read().splitlines()

# combos
with open("countCombos.txt", "w") as outf:
    for trans in transList:
        for ref in binList:
            tranbase = os.path.splitext(os.path.basename(trans))[0]
            outname = '{0}/{1}-counts.txt'.format(tranbase, ref)
            outf.write('{0}\t{1}\t{2}\n'.format(tranbase, ref, outname))