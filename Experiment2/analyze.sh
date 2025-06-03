for ID in flair stringtie2 IsoQuant;do
mkdir $1/Experiment2/$ID/gffcompare
mkdir $1/Experiment2/$ID/gffcompare/comp_to_ref
mkdir $1/Experiment2/$ID/gffcompare/comp_to_novel
mkdir $1/Experiment2/$ID/gffcompare/comp_to_as
done

flair correct -q $1/Experiment2/minimap2_output/all.bed12 -f $1/reference_data/gencode.v48.annotation.gtf -g $1/reference_data/GRCh38.primary_assembly.genome.fa --output $1/Experiment2/flair/corrected/all
flair collapse -g  $1/reference_data/GRCh38.primary_assembly.genome.fa -q $1/Experiment2/flair/corrected/all_all_corrected.bed -r $1/novel_isoform_data/all.ccs.fa --gtf $1/reference_data/gencode.v48.annotation.gtf --output $1/Experiment2/flair/all
mv $1/Experiment2/flair/all.isoforms.gtf $1/Experiment2/flair/all.guided.gtf

cd "$1/stringtie2" || { echo "Directory not found!"; exit 1; }
			
./stringtie -L $1/Experiment2/minimap2_output/all.sorted.bam -G $1/reference_data/gencode.v48.annotation.gtf -o $1/Experiment2/stringtie2/all.guided.gtf
./stringtie -L $1/Experiment2/minimap2_output/all.sorted.bam -o $1/Experiment2/stringtie2/all.unguided.gtf

cd "$1/IsoQuant" || { echo "Directory not found!"; exit 1; }
			
python isoquant.py --reference $1/reference_data/GRCh38.primary_assembly.genome.fa --bam $1/Experiment2/minimap2_output/all.sorted.bam --data_type pacbio_ccs --genedb $1/reference_data/gencode.v48.annotation.gtf --threads 3 -o $1/Experiment2/IsoQuant/all.guided
sleep(30)
python isoquant.py --reference $1/reference_data/GRCh38.primary_assembly.genome.fa --bam $1/Experiment2/minimap2_output/all.sorted.bam --data_type pacbio_ccs -o $1/Experiment1/IsoQuant/all.unguided

mv $1/Experiment2/IsoQuant/all.guided/OUT/OUT.transcript_models.gtf $1/Experiment2/IsoQuant/all.guided.gtf
mv $1/Experiment2/IsoQuant/all.unguided/OUT/OUT.transcript_models.gtf $1/Experiment2/IsoQuant/all.unguided.gtf

for ID in flair stringtie2 IsoQuant;do
for mode in guided unguided;do

cd "$1/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r $1/reference_data/gencode.v48.annotation.gtf $1/Experiment2/$ID/all.$mode.gtf -o $1/Experiment2/$ID/gffcompare/comp_to_ref/all.$mode
./gffcompare -r $1/reference_data/gencode.v48.annotation.gtf $1/Experiment2/$ID/all.$mode.gtf -o $1/Experiment2/$ID/gffcompare/comp_to_ref/all.$mode

cd "$1" || { echo "Directory not found!"; exit 1; }
python 5_29_ID_predicted_novel.py "$1/Experiment2/$ID/gffcompare/comp_to_ref/all.$mode.tracking" "$1/Experiment2/$ID/all.$mode.gtf" "$1/Experiment2/$ID/all.$mode.predicted_novel.gtf"
	
cd "$1/gffcompare" || { echo "Directory not found!"; exit 1; }
./gffcompare -r $1/4_10_yasim_data/4_10_25_human_p14_chr1_22.novel.gtf $1/Experiment1/$ID/novel_$s%_$mode.predicted_novel.gtf -o $1/Experiment2/$ID/gffcompare/comp_to_novel/all.$mode
done
done
