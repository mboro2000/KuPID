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

    //println!("Number of bins is {}", nb);

    let C_prob = 0.99;
    let mut C = 0.0;
    let mut check = 0;

    while C < C_prob{
        C += p.powi(check) * (1.0-p);
        //println!("{}", C);
        check += 1;
    }

    //println!("p is {}", p);
    //println!("C is {}, check is {}", C, check);

    /*
    let mut sample:HashMap<i64, Vec<i32>> = HashMap::new();
    let mut ref_hash:HashMap<i64, Vec<i32>> = HashMap::new();

    //Ref: a b c d e f a c
    //Q2:  a e c d - - - - - - - - - e f a c

    ref_hash.entry(1).or_insert([1, 7].to_vec());
    ref_hash.entry(2).or_insert([2].to_vec());
    ref_hash.entry(3).or_insert([3, 8].to_vec());
    ref_hash.entry(4).or_insert([4].to_vec());
    ref_hash.entry(5).or_insert([5].to_vec());
    ref_hash.entry(6).or_insert([6].to_vec());

    sample.entry(1).or_insert([1, 16].to_vec());
    sample.entry(3).or_insert([3, 17].to_vec());
    sample.entry(4).or_insert([4].to_vec());
    sample.entry(5).or_insert([2, 14].to_vec());
    sample.entry(6).or_insert([15].to_vec());

    let (chain, gap, score) = map::kmer_chain(&sample, &ref_hash, 3);
    println!("{}, {}", chain, gap);


    println!("{:#?}", anchors);
    //k=27, s=0.1, e = 0.0025, m = 3, b = 12, a = 11, c = 0.1, t = 3, q = 0.5
    */
    env::set_var("RUST_BACKTRACE", "1");

    let chks = [check].to_vec();
    //let chks = [3, 5, 7, 10].to_vec();

    for ck in chks.iter(){

        let now = Instant::now();
    
        let sample_data = sketch::read_input(t, &sample_file);
        let mut sample_data_shared = Arc::new(RwLock::new(sample_data));
        //println!("Read in the RNAseq reads");
    
        let ref_sketch:Arc<RwLock<HashMap<String, Sketch>>>   = sketch::get_sketches(sketch::read_input(t, &ref_file), k, a, s, t);
        //println!("Finished the unlock; Read and Sketched the reference reads");
        let elapsed = now.elapsed();    
        //println!("Elapsed: {:.2?}", elapsed); 
        let (ref_map, crit_1_shared) = map::find_ref_matches(ref_sketch, Arc::clone(&sample_data_shared), b, t, m, k, a, s, *ck);
        //println!("Finished Mapping");
        let elapsed = now.elapsed();    
        //println!("Elapsed: {:.2?}", elapsed); 

        //println!("Num. of exon gaps is {}", crit_1_shared.read().unwrap().len());
        let num_exon_gaps = crit_1_shared.read().unwrap().len();

    
        let criteria_2 = filter::filter(Arc::clone(&sample_data_shared), Arc::clone(&ref_map), nb, q, c, num_exon_gaps);

        let out = "novel_candidates_ck".to_string() + &ck.to_string() + ".fa";
        let candidates_path = Path::new(&output);
        let mut candidates_outline = "".to_string();

        let sample = sample_data_shared.read().unwrap();
    
        for (id, chunk) in crit_1_shared.read().unwrap().iter(){
            candidates_outline.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
        }
        for (id, chunk) in criteria_2.iter(){
            candidates_outline.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
        }

        //candidates_outline.push_str(&(samp_match.sample_id.clone() + "\n" + &sample_data[samp_match.chunk].get(&samp_match.sample_id).expect("msg") + "\n"));
        //candidates_outline.push_str(&(read.sample_id.clone() + "\n" + &sample_data[read.chunk].get(&read.sample_id).expect("msg") + "\n"));
        fs::write(candidates_path, candidates_outline);

        let elapsed = now.elapsed();    
        println!("Elapsed: {:.2?}", elapsed); 
    }

    
        
    /*
    
    let qs:Vec<f32> = [0.0, 0.1, 0.25, 0.5, 1.0].to_vec();
    //let qs:Vec<f32> = [0.5].to_vec();

    for q in qs.iter(){
        println!("Filter w/ q = {}", q);
        let criteria_2 = filter::filter(Arc::clone(&sample_data_shared), Arc::clone(&ref_map), nb, *q, c, num_exon_gaps);
    
        let out = "novel_candidates_q".to_string() + &q.to_string() + "_t2.fa";
        let candidates_path = Path::new(&out);
        let mut candidates_outline = "".to_string();

        let sample = sample_data_shared.read().unwrap();
    
        for (id, chunk) in crit_1_shared.read().unwrap().iter(){
            candidates_outline.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
        }
        for (id, chunk) in criteria_2.iter(){
            candidates_outline.push_str(&(">".to_string() + &id.clone() + "\n" + &sample[*chunk].get(id).expect("msg") + "\n"));
        }

        //candidates_outline.push_str(&(samp_match.sample_id.clone() + "\n" + &sample_data[samp_match.chunk].get(&samp_match.sample_id).expect("msg") + "\n"));
        //candidates_outline.push_str(&(read.sample_id.clone() + "\n" + &sample_data[read.chunk].get(&read.sample_id).expect("msg") + "\n"));
        fs::write(candidates_path, candidates_outline);
     */    
         
}
