#!/bin/bash
mkdir /usr1/mborowia/Experiment1

for s in 20 40 60 80 100;do
./minimap2 -ax splice --MD /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa /usr1/mborowia/Experiment1/novel_samples/novel_$s%.css.fa -t 3 -o /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sam      
	
samtools view -S -b /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sam > /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bam
samtools sort /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bam -o /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam
samtools index /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam
bedtools bamtobed -bed12 -i /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bam > /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bed12

flair correct -q /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bed12 -f /usr1/mborowia/reference_data/gencode.v48.annotation.gtf -g /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa --output /usr1/mborowia/Experiment1/flair/corrected/novel_$s%
flair collapse -g  /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa -q /usr1/mborowia/Experiment1/flair/corrected/novel_$s%_all_corrected.bed -r /usr1/mborowia/Experiment1/novel_samples/novel_$s%.css.fa --gtf /usr1/mborowia/reference_data/gencode.v48.annotation.gtf --output /usr1/mborowia/Experiment1/flair/novel_$s%
mv /usr1/mborowia/Experiment1/flair/novel_$s%.isoforms.gtf /usr1/mborowia/Experiment1/flair/novel_$s%_guided.gtf

cd "/usr1/mborowia/stringtie2" || { echo "Directory not found!"; exit 1; }
			
./stringtie -L /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam -G /usr1/mborowia/reference_data/gencode.v48.annotation.gtf -o /usr1/mborowia/Experiment1/stringtie2/novel_$s%_guided.gtf
./stringtie -L /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam -o /usr1/mborowia/Experiment1/stringtie2/novel_$s%_unguided.gtf

cd "/usr1/mborowia/IsoQuant" || { echo "Directory not found!"; exit 1; }
			
python isoquant.py --reference /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa --bam /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam --data_type pacbio_ccs --genedb /usr1/mborowia/reference_data/gencode.v48.annotation.gtf --threads 3 -o /usr1/mborowia/Experiment1/IsoQuant/novel_$s%_guided
python isoquant.py --reference /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa --bam /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.sorted.bam --data_type pacbio_ccs -o /usr1/mborowia/Experiment1/IsoQuant/novel_$s%_unguided

mv /usr1/mborowia/Experiment1/IsoQuant/novel_$s%_guided/OUT/OUT.transcript_models.gtf /usr1/mborowia/Experiment1/$ID/novel_$s%_guided.gtf
mv /usr1/mborowia/Experiment1/IsoQuant/novel_$s%_unguided/OUT/OUT.transcript_models.gtf /usr1/mborowia/Experiment1/$ID/novel_$s%_unguided.gtf

for ID in flair stringtie2 IsoQuant;do

mkdir /usr1/mborowia/Experiment1/$ID/gffcompare
mkdir /usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_ref
mkdir /usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_novel

for mode in guided unguided;do

cd "/usr1/mborowia/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r /usr1/mborowia/reference_data/gencode.v48.annotation.gtf /usr1/mborowia/Experiment1/$ID/novel_$s%_$mode.gtf -o /usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_ref/novel_$s%_$mode

cd "/usr1/mborowia" || { echo "Directory not found!"; exit 1; }
python 5_29_ID_predicted_novel.py "/usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_ref/novel_$s%_$mode.tracking" "/usr1/mborowia/Experiment1/$ID/novel_$s%_$mode.gtf" "/usr1/mborowia/Experiment1/$ID/novel_$s%_$mode.predicted_novel.gtf"
	
cd "/usr1/mborowia/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r /usr1/mborowia/4_10_yasim_data/4_10_25_human_p14_chr1_22.novel.gtf /usr1/mborowia/Experiment1/$ID/novel_$s%_$mode.predicted_novel.gtf -o /usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_novel/novel_$s%_$mode

done
done
done


for a in 0 1 2 5 10;do

./minimap2 -ax splice --MD /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa /usr1/mborowia/Experiment1/ratio_samples/ratio_1:$a.css.fa -t 3 -o /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sam      
	
samtools view -S -b /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sam > /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.bam
samtools sort /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.bam -o /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam
samtools index /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam
bedtools bamtobed -bed12 -i /usr1/mborowia/Experiment1/minimap2_output/novel_$s%.bam > /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.bed12

flair correct -q /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.bed12 -f /usr1/mborowia/reference_data/gencode.v48.annotation.gtf -g /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa --output /usr1/mborowia/Experiment1/flair/corrected/ratio_1:$a.
flair collapse -g /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa -q /usr1/mborowia/Experiment1/flair/corrected/ratio_1:$a._all_corrected.bed -r /usr1/mborowia/Experiment1/ratio_samples/1:$a.css.fa --gtf /usr1/mborowia/reference_data/gencode.v48.annotation.gtf --output /usr1/mborowia/Experiment1/flair/ratio_1:$a
mv /usr1/mborowia/Experiment1/flair/ratio_1:$a.isoforms.gtf /usr1/mborowia/Experiment1/flair/ratio_1:$a.guided.gtf

cd "/usr1/mborowia/stringtie2" || { echo "Directory not found!"; exit 1; }
			
./stringtie -L /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam -G /usr1/mborowia/reference_data/gencode.v48.annotation.gtf -o /usr1/mborowia/Experiment1/stringtie2/ratio_1:$a.guided.gtf
./stringtie -L /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam -o /usr1/mborowia/Experiment1/stringtie2/ratio_1:$a.unguided.gtf

cd "/usr1/mborowia/IsoQuant" || { echo "Directory not found!"; exit 1; }
			
python isoquant.py --reference /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa --bam /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam --data_type pacbio_ccs --genedb /usr1/mborowia/reference_data/gencode.v48.annotation.gtf --threads 3 -o /usr1/mborowia/Experiment1/IsoQuant/ratio_1:$a.guided
python isoquant.py --reference /usr1/mborowia/reference_data/GRCh38.primary_assembly.genome.fa --bam /usr1/mborowia/Experiment1/minimap2_output/ratio_1:$a.sorted.bam --data_type pacbio_ccs -o /usr1/mborowia/Experiment1/IsoQuant/ratio_1:$a.unguided

mv /usr1/mborowia/Experiment1/IsoQuant/ratio_1:$a.guided/OUT/OUT.transcript_models.gtf /usr1/mborowia/Experiment1/IsoQuant/ratio_1:$a.guided.gtf
mv /usr1/mborowia/Experiment1/IsoQuant/ratio_1:$a.unguided/OUT/OUT.transcript_models.gtf usr1/mborowia/Experiment1/IsoQuant/ratio_1:$a.unguided.gtf

for ID in flair stringtie2 IsoQuant;do
for mode in guided unguided;do
cd "/usr1/mborowia/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r /usr1/mborowia/reference_data/gencode.v48.annotation.gtf /usr1/mborowia/Experiment1/$ID/ratio_1:$a.$mode.gtf -o /usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_ref/ratio_1:$a.$mode

cd "/usr1/mborowia" || { echo "Directory not found!"; exit 1; }
python 5_29_ID_predicted_novel.py "/usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_ref/ratio_1:$a.$mode.tracking" "/usr1/mborowia/Experiment1/$ID/ratio_1:$a.$mode.gtf" "/usr1/mborowia/Experiment1/$ID/ratio_1:$a.$mode.predicted_novel.gtf"
	
cd "/usr1/mborowia/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r /usr1/mborowia/4_10_yasim_data/4_10_25_human_p14_chr1_22.novel.gtf /usr1/mborowia/Experiment1/$ID/ratio_1:$a.$mode.predicted_novel.gtf -o /usr1/mborowia/Experiment1/$ID/gffcompare/comp_to_novel/ratio_1:$a.$mode
done
done
done
