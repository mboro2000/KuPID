use std::collections::HashMap;
use std::time::Instant;
use std::env;
use clap::Parser;
use crate::types::*;
use std::sync::{Arc, RwLock};
use std::path::Path;
use std::fs;

use KuPID::types;
use KuPID::sketch;
use KuPID::map;
use KuPID::cmdline;

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
    let mode = cli.mode;
    let l = cli.l;
    let z = cli.z;

    env::set_var("RUST_BACKTRACE", "1");
    
    //Read in and convert the reference transcriptome into kmer sketches
    let ref_data = sketch::read_input(t, &ref_file);
    let ref_sketch:Arc<RwLock<HashMap<String, Sketch>>>   = sketch::get_sketches(ref_data, k, s, t);  
    //Read in RNAseq sample
    let sample_data = sketch::read_input(t, &sample_file);
    let sample_data_shared = Arc::new(RwLock::new(sample_data));
    //Psuedo-align RNAseq reads to reference transcriptome    
    let selected = map::find_ref_matches(output.clone(), ref_sketch, Arc::clone(&sample_data_shared), b, t, n, k, s, m, e as f32, B, c, mode.clone(), l, z);       
    //Output the filtered reads
    let output_file = output.clone() + "." + &mode + ".fa";
    let output_path = Path::new(&output_file);
    let mut output_line = "".to_string();
    let sample = sample_data_shared.read().unwrap();
    for (id, chunk) in selected.read().unwrap().iter(){
        output_line.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
    }
    fs::write(output_path, output_line); 
    
    }

