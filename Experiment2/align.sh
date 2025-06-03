#!/bin/bash

#arg1 = main directory

mkdir $1/Experiment2
mkdir $1/Experiment2/minimap2_output

cd "$1/minimap2" || { echo "Directory not found!"; exit 1; }

./minimap2 -ax splice --MD $1/reference_data/GRCh38.primary_assembly.genome.fa $1/novel_isoform_data/all.ccs.fa -t 3 -o $1/minimap2_output/all.ccs.sam
