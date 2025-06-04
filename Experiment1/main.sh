#!/bin/bash

mkdir $1/Experiment1

./sample_files.sh $1/novel_isoform_data/novel.ccs.bam $1/novel_isoform_data/annotated.ccs.bam $1/Experiment1

mkdir $1/Experiment1/minimap2_output
./align_samples.sh $1
#./analyze.sh $1
