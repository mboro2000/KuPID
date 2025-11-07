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
    <li><a href="#quickstart">Quickstart</a></li>
    <li><a href="#output">Output</a></li>
    <li><a href="#citation">Citation</a></li>
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
<h3 align="left">Downloads reference data</h3>
```
mkdir reference_data
cd reference_data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
gzip -d gencode.v48.transcripts.fa.gz
gzip -d gencode.v48.annotation.gtf.gz
```
## Output
