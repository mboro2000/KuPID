#!/bin/bash

#arg1 = main directory

mkdir $1/Experiment2
mkdir $1/Experiment2/minimap2_output

cd "$1/minimap2" || { echo "Directory not found!"; exit 1; }

./minimap2 -ax splice --MD $1/reference_data/GRCh38.primary_assembly.genome.fa $1/novel_isoform_data/all.ccs.fa -t 3 -o $1/Experiment2/minimap2_output/all.ccs.sam

samtools view -S -b $1/Experiment2/minimap2_output/all.ccs.sam > $1/minimap2_output/all.ccs.bam
samtools sort $1/Experiment2/minimap2_output/all.ccs.bam -o $1/Experiment2/minimap2_output/all.ccs.sorted.bam
samtools index $1/Experiment2/minimap2_output/all.ccs.sorted.bam
bedtools bamtobed -bed12 -i $1/Experiment2/minimap2_output/all.ccs.bam > $1/Experiment2/minimap2_output/all.ccs.bed12
