#!/bin/bash

#arg1: ccs reads of novel transcripts
#arg2: ccs reads of annotated transcripts
#arg3: main directory (/usr1/mborowia)

# Print a welcome message
mkdir $3/Experiment1/novel_samples
for s in 20 40 60 80 100;do
python $3/Experiment1/novel_sample.py $1 $s "$3/Experiment1/novel_samples"
done

mkdir $3/Experiment1/ratio_samples
for a in 0 1 2 5 10;do
python $3/Experiment1/ratio_sample.py $1 $2 $s "$3/Experiment1/ratio_samples"
done
