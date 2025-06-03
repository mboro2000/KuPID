use std::collections::HashMap;
use std::time::Instant;
use std::env;
use clap::Parser;

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

    let nb = 10;

    println!("{}, {}, {}", ref_file, s, m);

    //k=27, s=0.1, e = 0.02, m = 3, b = 12, a = 11, c = 0.1, t = 3, q = 0.5

    env::set_var("RUST_BACKTRACE", "1");

    let now = Instant::now();

    

    let sample_data = sketch::read_input(t, &sample_file);
    println!("Read in the RNAseq reads");
    println!("Total of {} RNAseq reads", sample_data.len());


    let ref_sketch:HashMap<String, types::Sketch>  = sketch::get_sketches(sketch::read_input(t, &ref_file), k, a, s, t);
    println!("Finished the unlock; Read and Sketched the reference reads");
    let elapsed = now.elapsed();    
    println!("Elapsed: {:.2?}", elapsed); 
    let ref_map:HashMap<String, Vec<types::Match>> = map::find_ref_matches(ref_sketch, &sample_data, b, t, m, k, a, s);
    println!("Finished Mapping");
    let elapsed = now.elapsed();    
    println!("Elapsed: {:.2?}", elapsed); 
    filter::filter(&sample_data, ref_map, output, nb, m, q, c);
    
    let elapsed = now.elapsed();    
    println!("Elapsed: {:.2?}", elapsed);         
}
