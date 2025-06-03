#!/bin/bash

#arg1 = main directory

for ID in flair stringtie2 IsoQuant;do
mkdir $1/Experiment1/$ID/gffcompare
mkdir $1/Experiment1/$ID/gffcompare/comp_to_ref
mkdir $1/Experiment1/$ID/gffcompare/comp_to_novel
done

for s in 20 40 60 80 100;do

flair correct -q $1/Experiment1/minimap2_output/novel_$s%.bed12 -f $1/reference_data/gencode.v48.annotation.gtf -g $1/reference_data/GRCh38.primary_assembly.genome.fa --output $1/Experiment1/flair/corrected/novel_$s%
flair collapse -g  $1/reference_data/GRCh38.primary_assembly.genome.fa -q $1/Experiment1/flair/corrected/novel_$s%_all_corrected.bed -r $1/Experiment1/novel_samples/novel_$s%.css.fa --gtf $1/reference_data/gencode.v48.annotation.gtf --output $1/Experiment1/flair/novel_$s%
mv $1/Experiment1/flair/novel_$s%.isoforms.gtf $1/Experiment1/flair/novel_$s%_guided.gtf

cd "$1/stringtie2" || { echo "Directory not found!"; exit 1; }
			
./stringtie -L $1/Experiment1/minimap2_output/novel_$s%.sorted.bam -G $1/reference_data/gencode.v48.annotation.gtf -o $1/Experiment1/stringtie2/novel_$s%_guided.gtf
./stringtie -L $1/Experiment1/minimap2_output/novel_$s%.sorted.bam -o $1/Experiment1/stringtie2/novel_$s%_unguided.gtf

cd "$1/IsoQuant" || { echo "Directory not found!"; exit 1; }
			
python isoquant.py --reference $1/reference_data/GRCh38.primary_assembly.genome.fa --bam $1/Experiment1/minimap2_output/novel_$s%.sorted.bam --data_type pacbio_ccs --genedb $1/reference_data/gencode.v48.annotation.gtf --threads 3 -o $1/Experiment1/IsoQuant/novel_$s%_guided
python isoquant.py --reference $1/reference_data/GRCh38.primary_assembly.genome.fa --bam $1/Experiment1/minimap2_output/novel_$s%.sorted.bam --data_type pacbio_ccs -o $1/Experiment1/IsoQuant/novel_$s%_unguided

mv $1/Experiment1/IsoQuant/novel_$s%_guided/OUT/OUT.transcript_models.gtf $1/Experiment1/IsoQuant/novel_$s%_guided.gtf
mv $1/Experiment1/IsoQuant/novel_$s%_unguided/OUT/OUT.transcript_models.gtf $1/Experiment1/IsoQuant/novel_$s%_unguided.gtf

for ID in flair stringtie2 IsoQuant;do
for mode in guided unguided;do

cd "$1/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r $1/reference_data/gencode.v48.annotation.gtf $1/Experiment1/$ID/novel_$s%_$mode.gtf -o $1/Experiment1/$ID/gffcompare/comp_to_ref/novel_$s%_$mode

cd "$1" || { echo "Directory not found!"; exit 1; }
python 5_29_ID_predicted_novel.py "$1/Experiment1/$ID/gffcompare/comp_to_ref/novel_$s%_$mode.tracking" "$1/Experiment1/$ID/novel_$s%_$mode.gtf" "$1/Experiment1/$ID/novel_$s%_$mode.predicted_novel.gtf"
	
cd "$1/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r $1/4_10_yasim_data/4_10_25_human_p14_chr1_22.novel.gtf $1/Experiment1/$ID/novel_$s%_$mode.predicted_novel.gtf -o $1/Experiment1/$ID/gffcompare/comp_to_novel/novel_$s%_$mode

done
done
done


for a in 0 1 2 5 10;do

flair correct -q $1/Experiment1/minimap2_output/ratio_1:$a.bed12 -f $1/reference_data/gencode.v48.annotation.gtf -g $1/reference_data/GRCh38.primary_assembly.genome.fa --output $1/Experiment1/flair/corrected/ratio_1:$a.
flair collapse -g $1/reference_data/GRCh38.primary_assembly.genome.fa -q $1/Experiment1/flair/corrected/ratio_1:$a._all_corrected.bed -r $1/Experiment1/ratio_samples/1:$a.css.fa --gtf $1/reference_data/gencode.v48.annotation.gtf --output $1/Experiment1/flair/ratio_1:$a
mv $1/Experiment1/flair/ratio_1:$a.isoforms.gtf $1/Experiment1/flair/ratio_1:$a.guided.gtf

cd "$1/stringtie2" || { echo "Directory not found!"; exit 1; }
			
./stringtie -L $1/Experiment1/minimap2_output/ratio_1:$a.sorted.bam -G $1/reference_data/gencode.v48.annotation.gtf -o $1/Experiment1/stringtie2/ratio_1:$a.guided.gtf
./stringtie -L $1/Experiment1/minimap2_output/ratio_1:$a.sorted.bam -o $1/Experiment1/stringtie2/ratio_1:$a.unguided.gtf

cd "$1/IsoQuant" || { echo "Directory not found!"; exit 1; }
			
python isoquant.py --reference $1/reference_data/GRCh38.primary_assembly.genome.fa --bam $1/Experiment1/minimap2_output/ratio_1:$a.sorted.bam --data_type pacbio_ccs --genedb $1/reference_data/gencode.v48.annotation.gtf --threads 3 -o $1/Experiment1/IsoQuant/ratio_1:$a.guided
python isoquant.py --reference $1/reference_data/GRCh38.primary_assembly.genome.fa --bam $1/Experiment1/minimap2_output/ratio_1:$a.sorted.bam --data_type pacbio_ccs -o $1/Experiment1/IsoQuant/ratio_1:$a.unguided

mv $1/Experiment1/IsoQuant/ratio_1:$a.guided/OUT/OUT.transcript_models.gtf $1/Experiment1/IsoQuant/ratio_1:$a.guided.gtf
mv $1/Experiment1/IsoQuant/ratio_1:$a.unguided/OUT/OUT.transcript_models.gtf usr1/mborowia/Experiment1/IsoQuant/ratio_1:$a.unguided.gtf

for ID in flair stringtie2 IsoQuant;do
for mode in guided unguided;do
cd "$1/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r $1/reference_data/gencode.v48.annotation.gtf $1/Experiment1/$ID/ratio_1:$a.$mode.gtf -o $1/Experiment1/$ID/gffcompare/comp_to_ref/ratio_1:$a.$mode

cd "$1" || { echo "Directory not found!"; exit 1; }
python 5_29_ID_predicted_novel.py "$1/Experiment1/$ID/gffcompare/comp_to_ref/ratio_1:$a.$mode.tracking" "$1/Experiment1/$ID/ratio_1:$a.$mode.gtf" "$1/Experiment1/$ID/ratio_1:$a.$mode.predicted_novel.gtf"
	
cd "$1/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r $1/4_10_yasim_data/4_10_25_human_p14_chr1_22.novel.gtf $1/Experiment1/$ID/ratio_1:$a.$mode.predicted_novel.gtf -o $1/Experiment1/$ID/gffcompare/comp_to_novel/ratio_1:$a.$mode
done
done
done
