use std::collections::HashMap;
use std::env;
use clap::Parser;
use crate::types::*;
use std::sync::{Arc, RwLock};
use std::path::Path;
use std::fs;
use std::process;

use KuPID::types;
use KuPID::sketch;
use KuPID::map;
use KuPID::cmdline;

fn invalid_range(x:f64, name:&str){
    if (x < 0.0) || (x > 1.0){
        eprintln!("Invalid range for parameter {}. Please select a value between 0 and 1.", name);
            std::process::exit(1);
    }
}

fn main() {      
    let cli = cmdline::Cli::parse();
    let ref_file = cli.reference;
    let sample_file = cli.input;
    let output = cli.output;
    let t = cli.t;
    let k = cli.k;
    let s = cli.s;
    let e = cli.e;
    let b = cli.b;
    let m = cli.m;
    let c = cli.c;
    let B = cli.B;
    let n = cli.n;
    let d = cli.d;
    let mode = cli.mode;
    let l = cli.l;
    let buffer = cli.f;

    invalid_range(B as f64, "B");
    invalid_range(e, "e");
    invalid_range(s, "s");

    if (mode != "quantify") && (mode != "discovery"){
        eprintln!("Invalid mode selected. Please choose either 'discovery' or 'quantify'");
            std::process::exit(1);
    }
    if k as f64 * e >= 1.0{
        eprintln!("Kmer length chosen is non-ideal. Given your error rate, all of your kmers are expected to be erroneous. We recommend selecting k such that k*e < 1.");
            std::process::exit(1);
    }

    env::set_var("RUST_BACKTRACE", "1");
    
    //Read in and convert the reference transcriptome into kmer sketches
    let (ref_data, transcript_names, ref_chunks) = sketch::read_input(t, &ref_file, true);
    
    let ref_data_shared: Arc<RwLock<Vec<HashMap<i32, String>>>> = Arc::new(RwLock::new(ref_data));
    let ref_sketch:Arc<RwLock<HashMap<i32, Sketch>>>  = sketch::get_sketches(Arc::clone(&ref_data_shared), k, s, t);  
    
    let ref_chunks_shared: Arc<RwLock<HashMap<i32, usize>>> = Arc::new(RwLock::new(ref_chunks));
    //Read in RNAseq sample
    let (sample_data, query_names, _query_chunks) = sketch::read_input(t, &sample_file, true);
    
    let sample_data_shared = Arc::new(RwLock::new(sample_data));
    //Psuedo-align RNAseq reads to reference transcriptome    
    let selected = map::find_ref_matches(output.clone(), ref_sketch, Arc::clone(&ref_data_shared), Arc::clone(&ref_chunks_shared), Arc::clone(&sample_data_shared), b, t, n, k, s, m, e as f32, B, c, mode.clone(), l, transcript_names, buffer, d as f32);       
    //Output the filtered reads
    let output_file = output.clone() + "." + &mode + ".fa";
    let output_path = Path::new(&output_file);
    let mut output_line = "".to_string();
    let sample = sample_data_shared.read().unwrap();

    for (id, chunk) in selected.read().unwrap().iter(){
        let query_name = query_names.get(id).expect("msg");
        output_line.push_str(&(">".to_string() + &query_name.to_string() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
    }
    fs::write(output_path, output_line);     
    }

