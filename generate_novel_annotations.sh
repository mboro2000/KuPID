#!/bin/bash

#arg1 = main directory

mkdir $1/reference_data
cd $1/reference_data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz

#release 48 of human genome via GenCode

mkdir $1/novel_isoform_data
python -m yasim generate_as_events -f $1/reference_data/GRCh38.primary_assembly.genome.fa -g $1/reference_data/gencode.v48.annotation.gtf -o $1/novel_isoform_data/gencode.v48.as.gtf -c 2
python -m yasim generate_gene_depth -g $1/novel_isoform_data/gencode.v48.as.gtf -o $1/novel_isoform_data/gencode.v48.gene_depth.tsv -d 20
python -m yasim generate_isoform_depth -g $1/novel_isoform_data/gencode.v48.as.gtf -d $1/novel_isoform_data/gencode.v48.gene_depth.tsv -o $1/novel_isoform_data/gencode.v48.isoform_depth.tsv
python -m labw_utils.bioutils transcribe -f $1/reference_data/GRCh38.primary_assembly.genome.fa -g $1/novel_isoform_data/gencode.v48.as.gtf -o $1/novel_isoform_data/gencode.v48.as.fa

#install fastaparser
python make_transcript_file.py $1/novel_isoform_data/gencode.v48.isoform_depth.tsv $1/novel_isoform_data/gencode.v48.as.fa $1/novel_isoform_data/
