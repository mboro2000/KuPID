#!/bin/bash

#arg1 = main directory

cd "$1/minimap2" || { echo "Directory not found!"; exit 1; }

./minimap2 -ax splice --MD ~/human_p14_chr_1_22.fa all_chr1_22.ccs.fa -t 3 -o 4_13_ccs_chr1_22_all_t3.sam
