#! /bin/bash 

tar -xvf prodigal.tar.gz

base=$(basename $1)
mkdir $base-NUCS
cd $base
for file in *.fna; do
    ../prodigal/./prodigal -i $file -o ../$base-NUCS/$file-genes.fna;
done

cd ..

tar -cvf $base-NUCS.tar.gz $base-NUCS/