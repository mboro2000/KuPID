<a id="readme-top"></a>

<br />
<div align="left">
  <h1 align="left">KuPID: Kmer-based Upstream Preprocessing for Isoform Discovery</h1>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#introduction">Introduction</a> </li>
    <li><a href="#install">Install</a></li>
    <li><a href="#tutorial">Sample Tutorial</a></li>
  </ol>
</details>

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

## Sample Tutorial

To demonstrate KuPID's abilities, we've provided an RNAseq sample of PacBio HiFi reads sequenced from 3500 genes in the human genome. In addition, we've provided an annotation file of the novel isoforms present in the sample.

To complete the downstream analysis, users should have access to the following software:
1. minimap2
2. stringtie2
3. gffcompare

<h4 align="left">Download reference data</h4>

```
mkdir reference_data
cd reference_data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
gzip -d gencode.v48.transcripts.fa.gz
gzip -d gencode.v48.annotation.gtf.gz
```

<h4 align="left">Apply KuPID</h4>

```
cd KuPID/src
for mode in discovery quantify;do
cargo run -- -o "small_sample" -r "~/reference_data/gencode.v48.transcripts.fa" -i "~/sample/5000_genes.ccs.fasta" -k 22 -s 0.1 -e 0.002 -n 30 -b 16 -m 3 -B 0.98 -t 1 -c 1.5 -l 5 -g 100 --mode $mode
done
```

<h4 align="left">Downstream Alignment and Isoform Discovery</h4>

```
#Align non-processed reads
./minimap2 -ax splice:hq -uf --MD ~/reference_data/GRCh38.chr1-22.fa ~/sample/small_sample.ccs.fasta -t 3 -o ~/sample/minimap2_output/small_sample.sam
samtools view -S -b ~/sample/minimap2_output/small_sample.sam > ~/sample/minimap2_output/small_sample.bam
samtools sort ~/sample/minimap2_output/small_sample.bam -o ~/sample/minimap2_output/small_sample.sorted.bam
samtools index ~/sample/minimap2_output/small_sample.sorted.bam

#Align KuPID-processed reads
for mode in discovery quantify;do
./minimap2 -ax splice:hq -uf --MD ~/reference_data/GRCh38.chr1-22.fa ~/KuPID/src/small_sample.$mode.fa -t 1 -o /usr1/mborowia/11_6_KuPID/minimap2_output/small_sample.$mode.sam
samtools view -S -b ~/sample/minimap2_output/small_sample.$mode.sam > ~/sample/minimap2_output/small_sample.$mode.bam
samtools sort ~/sample/minimap2_output/small_sample.$mode.bam -o ~/sample/minimap2_output/small_sample.$mode.sorted.bam
samtools index ~/sample/minimap2_output/small_sample.$mode.sorted.bam
done

#Apply stringtie2
./stringtie -L ~/sample/minimap2_output/small_sample.sorted.bam -G ~/reference_data/gencode.v48.annotation.gtf -o ~/sample/stringtie2.small_sample.gtf
for mode in discovery quantify;do
./stringtie -L ~/sample/minimap2_output/small_sample.$mode.sorted.bam -G ~/reference_data/gencode.v48.annotation.gtf -o ~/sample/stringtie2.small_sample.$mode.gtf
done

```

<h4 align="left">Analyze Discovery Results</h4>


## Output
