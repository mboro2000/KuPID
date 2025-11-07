use std::collections::HashMap;
use std::ptr::read;
use std::time::Instant;
use std::env;
use clap::Parser;
use crate::types::*;
use std::sync::{Arc, RwLock};
use std::path::Path;
use std::fs;

use isoforms::types;
use isoforms::sketch;
use isoforms::map;
//use isoforms::filter;
use isoforms::cmdline;

use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

fn main() {      

    let cli = cmdline::Cli::parse();

    let ref_file = cli.reference;
    let sample_file = cli.input;
    let output = cli.output;
    let t = cli.t;
    let k = cli.k;
    let s = cli.s;
    let e = cli.e;
    let a = cli.a;
    let b = cli.b;
    let m = cli.m;
    let c = cli.c;
    let B = cli.B;
    let n = cli.n;

    env::set_var("RUST_BACKTRACE", "1");

    let now = Instant::now();
    
    let (sample_data, total_reads) = sketch::read_input(t, &sample_file);
    let sample_data_shared = Arc::new(RwLock::new(sample_data));

    let ref_data = sketch::read_input(t, &ref_file).0;
    
    let ref_sketch:Arc<RwLock<HashMap<String, Sketch>>>   = sketch::get_sketches(sketch::read_input(t, &ref_file).0, k, a, s, t);  
    
    let (AS_shared, ATSS, AS_novel, AS_annot, ATSS_novel, ATSS_annot) = map::find_ref_matches(output.clone(), ref_sketch, Arc::clone(&sample_data_shared), b, t, n, k, a, s, m, total_reads, e as f32, B, c);       
    let sample = sample_data_shared.read().unwrap();

    let AS_file = output.clone() + ".AS.fa";
    let ATSS_file = output.clone() + ".ATSS.fa";

    let candidates_path_AS = Path::new(&AS_file);
    let mut candidates_outline_AS = "".to_string();
    let candidates_path_ATSS = Path::new(&ATSS_file);
    let mut candidates_outline_ATSS = "".to_string();
    
         
    for (id, chunk) in AS_shared.read().unwrap().iter(){
        candidates_outline_AS.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
        candidates_outline_ATSS.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
    }
    for (id, chunk) in ATSS.iter(){
        candidates_outline_ATSS.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
    }

        
    fs::write(candidates_path_AS, candidates_outline_AS);
    fs::write(candidates_path_ATSS, candidates_outline_ATSS);              

    let elapsed = now.elapsed();    

    let stats_file = output.clone() + ".stats.csv";
    let stats_path = Path::new(&stats_file);
    let mut stats_outline = "AS novel,AS annot,ATSS novel,ATSS annot,runtime\n".to_string();
    stats_outline.push_str(&(AS_novel.to_string() + "," + &AS_annot.to_string() + "," + &ATSS_novel.to_string() + "," + &ATSS_annot.to_string() + ","  + &elapsed.as_secs_f32().to_string()));


    println!("{}", AS_novel);
    println!("{}", AS_annot);
    println!("{}", ATSS_novel);
    println!("{}", ATSS_annot);

    fs::write(stats_path, stats_outline);
    
    }
