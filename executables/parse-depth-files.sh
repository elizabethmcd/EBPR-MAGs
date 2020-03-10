#! /bin/bash 

# Parse depth mapping files to get coverage results

REF=$1
REFBASE=$(basename $REF .fna)
META=$2
METABASE=$(basename $META .qced.fastq)
DEPTH=$3
DEPTHBASE=$(basename $DEPTH .qced.sorted.dept)
EXDIR=/home/GLBRCORG/emcdaniel/EBPR/coassembly/scripts
OUTDIR=/home/GLBRCORG/emcdaniel/EBPR/coassembly/metagenomes/finalBins/depthResults

# Reference length

python $EXDIR/countBases.py $REF > $OUTDIR/$REFBASE.len

# Metagenomes reads file

awk '{s++}END{print FILENAME,s/4}' $META >> $OUTDIR/$METABASE-reads.txt

# stats file

python $EXDIR/calc-mapping-stats.py $OUTDIR/$REFBASE.len $OUTDIR/$METABASE-reads.txt $DEPTH $OUTDIR/$DEPTHBASE.coverage.txt

