#! /bin/bash

folder=$(basename $PWD -EBPR-bins)
for file in *.fna; do
    bin=$(basename $file .fna);
    mv $file ../$folder/$folder-$bin.fna;
done