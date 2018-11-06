#! /bin/bash
 
 # split large metatranscriptome fastqs into smaller subsets with BBMap to filter out rRNA faster with sortmeRNA

# setup 
tar -xzvf BBMap_38.07.tar.gz
cp $1 .
tar -xzvf $1 
file=$(basename $1 .tar.gz)
name=$(basename $1 .qced.fastq.tar.gz)

sed -i '/^$/d' $file
maxreads=$((`wc -l < $name` / 8 - 1))
startpoints=$(seq 0 500000 $maxreads)

for num in $startpoints;
    do endpoint=$(($num + 499999));
    bbmap/getreads.sh in=$file id=$num-$endpoint out=$name-$endpoint.fastq overwrite=T;
done 

rm $file 
mkdir $name-splitfiles 
mv $name*.fastq $name-splitfiles
tar -czvf $name-splitfiles.tar.gz $name-splitfiles/
cp -r $name-splitfiles/ /mnt/gluster/emcdaniel/