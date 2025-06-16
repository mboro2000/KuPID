use std::collections::HashMap;
use std::time::Instant;
use std::env;
use clap::Parser;
use std::sync::Mutex;
use crate::types::*;
use std::sync::{Arc, RwLock};
use std::path::Path;
use std::fs;
use std::f32;

use isoforms::types;
use isoforms::sketch;
use isoforms::map;
use isoforms::filter;
use isoforms::cmdline;
   
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
    let q = cli.q;

    let nb = (10.0 / (1.0 - k as f64*e)).floor() as i32;
    let p = k as f64*e * (2.0 + (k as f64 * e - 1.0) * (1.0 - s).powi(k - 1) - k as f64*e);

    let C_prob = 0.99;
    let mut C = 0.0;
    let mut check = 0;

    while C < C_prob{
        C += p.powi(check) * (1.0-p);
        check += 1;
    }

    //println!("p is {}", p);
    //println!("C is {}, check is {}", C, check);

    env::set_var("RUST_BACKTRACE", "1");

    let chks = [check].to_vec();
    //let chks = [3, 5, 7, 10].to_vec();

    for ck in chks.iter(){

        let now = Instant::now();
    
        let sample_data = sketch::read_input(t, &sample_file);
        let mut sample_data_shared = Arc::new(RwLock::new(sample_data));    
        let ref_sketch:Arc<RwLock<HashMap<String, Sketch>>>   = sketch::get_sketches(sketch::read_input(t, &ref_file), k, a, s, t);
        let (ref_map, crit_1_shared) = map::find_ref_matches(ref_sketch, Arc::clone(&sample_data_shared), b, t, m, k, a, s, *ck);  
        let num_exon_gaps = crit_1_shared.read().unwrap().len();
        let criteria_2 = filter::filter(Arc::clone(&sample_data_shared), Arc::clone(&ref_map), nb, q, c, num_exon_gaps);
        //let out = "novel_candidates_ck".to_string() + &ck.to_string() + ".fa";
        let candidates_path = Path::new(&output);
        let mut candidates_outline = "".to_string();

        let sample = sample_data_shared.read().unwrap();
    
        for (id, chunk) in crit_1_shared.read().unwrap().iter(){
            candidates_outline.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
        }
        for (id, chunk) in criteria_2.iter(){
            candidates_outline.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
        }

        fs::write(candidates_path, candidates_outline);

        let elapsed = now.elapsed();    
        println!("Elapsed: {:.2?}", elapsed); 
    }
