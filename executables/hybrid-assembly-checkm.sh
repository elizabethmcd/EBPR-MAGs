#!/bin/bash

mkdir input
cp /mnt/gluster/emcdaniel/EBPR-Mapping-Results/hybrid-assembly-bins/*.fa .
mkdir hybrid-checkm
for filename in *.fa; do cp $filename ${filename%.fa}.fna; done
for filename in *.fna; do cp $filename input/; done

# Running lineage specific marker set

checkm lineage_wf input output

tax=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /marker lineage/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
length=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Genome size/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
complete=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Completeness/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
contamination=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Contamination/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')

mv output/ hybrid-checkm
tar -czf hybrid-checkm.tgz hybrid-checkm/
mv hybrid-checkm.tgz /mnt/gluster/emcdaniel/EBPR-Mapping/.