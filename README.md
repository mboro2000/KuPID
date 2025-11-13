<a id="readme-top"></a>

<br />
<div align="left">
  <h1 align="left">KuPID: Kmer-based Upstream Preprocessing for Isoform Discovery</h1>
</div>

<!-- intro -->
## Introduction

Eukaryotic genes can encode multiple protein isoforms based on alternative splicing of exonic and intronic regions. More than 95\% of human genes have been found to undergo alternative splicing. Alternative isoforms have been implicated in a wide-range of biological processes, such cell differentiation, stress response, and tumorigenesis. Most modern novel isoform discovery methods function by identifying and assembling exon splice junctions from an RNAseq sample. However, splice junctions can only be accurately annotated with time-intensive dynamic programming alignment.

KuPID is a preprocessing method designed for RNAseq analysis of long transcript reads. When given an RNAseq sample, KuPID will filter out the reads expressed from reference (annotated) isoforms. The remaining reads can be submitted to the downstream alignment and isoform discovery software of the user's choice. KuPID operates by applying kmer sketching methods to quickly pseudo-align RNAseq reads to a reference transcriptome. The KuPID-processed samples require less time for downstream alignment, and routinely discover more true novel isoforms.

<!-- install -->
## Install

KuPID is currently available to build from source as a rust crate. 

Requirements:
1. rust programming language and associated tools such as cargo are required and assumed to be in PATH.
2. A c compiler (e.g. GCC)
3. make

```
git clone https://github.com/mboro2000/KuPID.git
cd KuPID

# If default rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo

# If ~/.cargo doesn't exist use below commands instead
#cargo build --release
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Input

To function, KuPID requires the following inputs:
1. RNAseq sample (-i); should be formatted as a .fasta of the sequenced reads
2. Reference transcriptome (-r); should be formatted of a .fasta of the reference isoforms

## Sample Tutorial

To demonstrate KuPID's abilities, we've provided an RNAseq sample of PacBio HiFi reads sequenced from chr1 of the human genome. In addition, we've provided an annotation file of the novel isoforms present in the sample.

To complete the downstream analysis, users should have access to the following software:
1. minimap2
2. stringtie2
3. gffcompare

<h4 align="left">Download reference data</h4>

```
cd KuPID/reference_data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
gzip -d gencode.v48.transcripts.fa.gz
gzip -d gencode.v48.annotation.gtf.gz
gzip -d Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
python edit_heading.py Homo_sapiens.GRCh38.dna.chromosome.1.fa "chr1"
```

<h4 align="left">Apply KuPID</h4>

```
cd KuPID/src
for mode in discovery quantify;do
cargo run -- -o "small_sample" -r "~/KuPID/reference_data/gencode.v48.transcripts.fa" -i "~/KuPID/sample/small_sample.ccs.fasta" -k 22 -s 0.1 -e 0.002 -n 30 -b 16 -m 3 -B 0.98 -t 1 -c 1.5 -l 5 -z 100 --mode $mode
done
```

<h4 align="left">Downstream Alignment and Isoform Discovery</h4>

```
cd KuPID/sample
mkdir minimap2_output
#Align non-processed reads
./minimap2 -ax splice:hq -uf --MD ~/KuPID/reference_data/Homo_sapiens.GRCh38.dna.chromosome.1.fa ~/KuPID/sample/small_sample.ccs.fasta -t 3 -o ~/KuPID/sample/minimap2_output/small_sample.sam
samtools view -S -b ~/KuPID/sample/minimap2_output/small_sample.sam > ~/KuPID/sample/minimap2_output/small_sample.bam
samtools sort ~/KuPID/sample/minimap2_output/small_sample.bam -o ~/KuPID/sample/minimap2_output/small_sample.sorted.bam
samtools index ~/KuPID/sample/minimap2_output/small_sample.sorted.bam

#Align KuPID-processed reads
for mode in discovery quantify;do
./minimap2 -ax splice:hq -uf --MD ~/KuPID/reference_data/Homo_sapiens.GRCh38.dna.chromosome.1.fa ~/KuPID/src/small_sample.$mode.fa -t 1 -o ~/KuPID/sample/minimap2_output/small_sample.$mode.sam
samtools view -S -b ~/KuPID/sample/minimap2_output/small_sample.$mode.sam > ~/KuPID/sample/minimap2_output/small_sample.$mode.bam
samtools sort ~/KuPID/sample/minimap2_output/small_sample.$mode.bam -o ~/KuPID/sample/minimap2_output/small_sample.$mode.sorted.bam
samtools index ~/KuPID/sample/minimap2_output/small_sample.$mode.sorted.bam
done

#Apply stringtie2
./stringtie -L ~/KuPID/sample/minimap2_output/small_sample.sorted.bam -G ~/KuPID/reference_data/gencode.v48.annotation.gtf -o ~/KuPID/sample/stringtie2.small_sample.gtf
for mode in discovery quantify;do
./stringtie -L ~/KuPID/sample/minimap2_output/small_sample.$mode.sorted.bam -G ~/KuPID/reference_data/gencode.v48.annotation.gtf -o ~/KuPID/sample/stringtie2.small_sample.$mode.gtf
done

```

<h4 align="left">Scale Quantification Results (Required for KuPID-quantify)</h4>

```
### scale quantification results for KuPID reads
python ~/scale_quantification_results.py -a ~/KuPID/sample/stringtie2.small_sample.quantify.gtf --method stringtie2 --output ~/KuPID/sample/stringtie2.small_sample.KuPID.scaled_tpm.csv --scale ~/KuPID/src/small_sample.scale_factors.csv -l 5 -p KuPID
### adjust quantifcation results of non-processed reads to only show abundances of known reference transcripts
python ~/scale_quantification_results.py -a ~/KuPID/sample/stringtie2.small_sample.quantify.gtf --method stringtie2 --output ~/KuPID/sample/stringtie2.small_sample.tpm.csv --scale ~/KuPID/src/small_sample.scale_factors.csv -l 5 -p None

# -a: initial abundance results reported by the chosen quantification method
# --method: ID method used to quantify transcripts. KuPID can currently be paired with IsoQuant, flair, or stringtie2
# -l: number of reads sampled from each annotated isoform during KuPID-quantify
```

<h4 align="left">Analyze Discovery Results</h4>

```
cd KuPID/sample
mkdir gffcompare
for trial in stringtie2.small_sample stringtie2.small_sample.discovery stringtie2.small_sample.quantify;do
./gffcompare -r ~/KuPID/reference_data/gencode.v48.annotation.gtf ~/KuPID/sample/$trial.gtf -o ~/KuPID/sample/gffcompare/$trial.comp_to_ref
## Extract the transcript models that are predicted to be novel
python ~/KuPID/sample/predicted_novel.py "/KuPID/sample/gffcompare/$trial.comp_to_ref.tracking" "~/KuPID/sample/$trial.gtf" "~/KuPID/sample/$trial.predicted_novel.gtf"
./gffcompare -r ~/KuPID/sample/small_sample.novel.gtf ~/KuPID/sample/$trial.predicted_novel.gtf -o ~/KuPID/sample/gffcompare/$trial.comp_to_novel
done

### gffcompare will report the precision and accuracy of novel transcript models assembled by stringtie2. To view this output, use:
for trial in stringtie2.small_sample stringtie2.small_sample.discovery stringtie2.small_sample.quantify;do
cat ~/KuPID/sample/gffcompare/$trial.comp_to_novel
done
```
