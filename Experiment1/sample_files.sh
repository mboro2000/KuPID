#!/bin/bash

#arg1: ccs reads of novel transcripts
#arg2: ccs reads of annotated transcripts
#arg3: output directory (/usr1/mborowia/Experiment1)

# Print a welcome message
mkdir novel_samples
python 5_29_novel_sample.py $1 20 "$3/novel_samples"
python 5_29_novel_sample.py $1 40 "$3/novel_samples"
python 5_29_novel_sample.py $1 60 "$3/novel_samples"
python 5_29_novel_sample.py $1 80 "$3/novel_samples"
python 5_29_novel_sample.py $1 100 "$3/novel_samples"

mkdir ratio_samples
python 5_29_ratio_sample.py $1 $2 0 "$3/ratio_samples"
python 5_29_ratio_sample.py $1 $2 1 "$3/ratio_samples"
python 5_29_ratio_sample.py $1 $2 2 "$3/ratio_samples"
python 5_29_ratio_sample.py $1 $2 5 "$3/ratio_samples"
python 5_29_ratio_sample.py $1 $2 10 "$3/ratio_samples"
