#! /bin/bash 

# Parse depth mapping files to get coverage results

PYTHONPATH=/home/GLBRCORG/emcdaniel/anaconda3/bin/python3
REF=$1
REFBASE=$(basename $REF .fna)
META=$2
METABASE=$(basename $META .qced.fastq)
DEPTH=$3
DEPTHBASE=$(basename $DEPTH .qced.sorted.depth)
EXDIR=/home/GLBRCORG/emcdaniel/EBPR/coassembly/scripts
OUTDIR=/home/GLBRCORG/emcdaniel/EBPR/coassembly/metagenomes/finalBins/depthResults

# Reference length

python $EXDIR/countBases.py $REF $OUTDIR/$DEPTHBASE.len

# Metagenomes reads file

awk '{s++}END{print FILENAME,s/4}' $META >> $OUTDIR/$DEPTHBASE-reads.txt

# stats file

$PYTHONPATH $EXDIR/calc-mapping-stats.py $OUTDIR/$DEPTHBASE.len $OUTDIR/$DEPTHBASE-reads.txt $DEPTH $OUTDIR/$DEPTHBASE.coverage.txt
