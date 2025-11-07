use clap::{Parser};

#[derive(Parser)]
#[clap(author, version, about = "Preprocessing script to select reads from an RNAseq sample that are likely transcribed from novel isoforms", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    //#[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional. Default: fastq file produces a sample sketch (*.sylsp) while fasta files are combined into a database (*.syldb).")]
    //pub files: Vec<String>,
    //from sylph github
    #[clap(short='o',long="output path", default_value = "./", help_heading = "OUTPUT", help = "Output fasta file of the chosen novel read candidates")]
    pub output: String,
    #[clap(short='r',long="reference", default_value = "", help = "Input fasta file of reference transcripts")]
    pub reference: String,
    #[clap(short='i',long="input sample", default_value = "", help = "Input fasta of RNAseq sample reads")]
    pub input: String,

    #[clap(short='t',long="threads", default_value_t = 3, help = "Set number of threads to use")]
    pub t: i32,
    #[clap(short='k',long="kmer length", default_value_t = 12, help = "Set kmer length used for sketching")]
    pub k: i32,
    #[clap(short='e',long="Error rate", default_value_t = 0.1, help = "Set the average error rate of the input reads")]
    pub e: f64,
    #[clap(short='a',long="Hash seed value", default_value_t = 11, help = "Set the seed value for the FracMinHash function")]
    pub a: i64,
    #[clap(short='s',long="Sketch density", default_value_t = 0.1, help = "Set the density of the FracMinHash function")]
    pub s: f64,
    #[clap(short='b',long="Max density of kmer map bins", default_value_t = 12, help = "Set the maximum num. of isoforms that can map to a given kmer")]
    pub b: usize,
    #[clap(short='n',long="Gap between matches", default_value_t = 30, help = "Set the maximum gap allowed in an optimal chain without AS")]
    pub n: i32,  
    #[clap(short='m',long="Band width", default_value_t = 30, help = "Set the band width of the chaining procedure's dynamic programming table")]
    pub m: i32, 
    #[clap(short='B',long="Maximum novel similarity", default_value_t = 0.98, help = "Set the maximum similarity score allowed for a novel read")]
    pub B: f32,
    #[clap(short='c',long="Candidate set size", default_value_t = 1.5, help = "Set the scale factor for the number of ATSS candidates")]
    pub c: f64
}
