#!/bin/bash

mkdir Experiment1/minimap2_output
cd "/usr1/mborowia/minimap2" || { echo "Directory not found!"; exit 1; }

for s in 20 40 60 80 100;do
./minimap2 -ax splice --MD /usr1/mborowia/human_p14_chr_1_22.fa /usr1/mborowia/Experiment1/novel_samples/novel_$s%.css.fa -t 3 -o /usr1/mborowia/Experiment1/minim>#       samtools view -S -b /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sam > /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bam
samtools sort /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bam -o /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam
samtools index /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam
done

for a in 0 1 2 5 10;do
./minimap2 -ax splice --MD /usr1/mborowia/human_p14_chr_1_22.fa /usr1/mborowia/Experiment1/ratio_samples/1:$a.css.fa -t 3 -o /usr1/mborowia/Experiment1/minimap2_ou>        samtools view -S -b /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sam > /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.bam
samtools sort /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.bam -o /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam
samtools index /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam
done
