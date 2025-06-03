
use std::collections::HashMap;
use std::sync::Mutex;
use bio::io::fasta;
use std::str;
use std::cmp;
use crate::types::*;
use std::thread;
use std::sync::{Arc, RwLock};
use num::integer::div_ceil;

pub fn read_input(t:i32, file:&String) -> (Vec<HashMap<String, String>>){
    let mut seq_data:Vec<HashMap<String, String>> = Vec::new();    
    let read_results = fasta::Reader::from_file(file);

    let mut r = 0;
    for i in 0 .. t{
        let mut read_chunk:HashMap<String, String> = HashMap::new();
        seq_data.push(read_chunk);
    }

    for reader in read_results{
        for result in reader.records(){
        
            let record = result.expect("Error during fasta record parsing");
            let seq = str::from_utf8(record.seq()).unwrap().to_string();
            seq_data[(r % t) as usize].entry(record.id().to_string()).or_insert(seq);
            r += 1;
        }
    }  
    seq_data
}

pub fn get_sketches(read_data:Vec<HashMap<String, String>>, k:i32, a:i64, s:f64, t:i32) -> HashMap<String, Sketch>{
   
    let mut seq_sketches:HashMap<String, Sketch> = HashMap::new();
    let mut handles = vec![];    
    let mut seq_sketches_shared = Arc::new(Mutex::new(seq_sketches));
    let mut read_data_shared = Arc::new(RwLock::new(read_data));

    for i in 0..t{
        let seq_sketches_shared = Arc::clone(&seq_sketches_shared); 
        let read_data_shared = Arc::clone(&read_data_shared); 

        let handle = thread::spawn(move || {
            let guard = read_data_shared.read().unwrap();
            let seq_chunk = (guard).get(i as usize).expect("msg");

            let mut r_count = 0;

            for (id, seq) in seq_chunk{       
                let (num_kmers, kmer_set) = FracMinHash(&seq, k, a,s);
                let mut seq_sketches = seq_sketches_shared.lock().unwrap();
                seq_sketches.entry(id.clone()).or_insert(build_Sketch(id.clone(), num_kmers, kmer_set));
                drop(seq_sketches);
                   
                r_count += 1;
                if r_count % 1000 == 0{
                    println!("Input read {}", r_count);
                }
            }
        });
        handles.push(handle);
    }

    for i in handles{
        i.join().unwrap();
    } 

    println!("Beginning to unlock");
    let seq_sketches = seq_sketches_shared.lock().unwrap().clone();
    seq_sketches
}

pub fn map_ATCG(item:char, mut label:i64) -> i64{

    if 'A' == item || 'a' == item{
        label += 0;
    }
    else if 'C' == item || 'c' == item {
        label += 1;
    }
    else if 'G' == item || 'g' == item {
        label += 2;
    }
    else{
        label += 3;
    }
    label
}


//Returns a vector encoding of minimizers selected using the FracMinHash method
//table = <score of minimizer m, set of ordinal positions where m is located in the sketch>
pub fn FracMinHash(seq:&str, k:i32, a:i64, s:f64) -> (i32, HashMap<i64, Vec<i32>>) {

    let max:i64 = ((1 << (2*k)) - 1) >> 3;
    let Hs:f64 = max as f64 * s;

    let rounds = (seq.len() as i32 - k + 1);

    let mut v:Vec<i64> = Vec::new();
    let mut label:i64 = 0;
    let mut kmers:HashMap<i64, Vec<i32>> = HashMap::new();
    let mut ind = 0;

    let s = &seq[0 .. cmp::min(k as usize, seq.len())];

    for item in s.chars() {
        label  = label << 2;
        label = map_ATCG(item, label);            
    }  

    if seq.len() > k as usize{
        let mod_score = (a*label) & max;
    
        if mod_score as f64 <= Hs{
            kmers.entry(label).or_insert(Vec::new()).push(ind);
            ind += 1;
        }

        let mut exit_char = seq[ .. (seq.len() as i32 - k) as usize].chars();
        let mut enter_char = seq[k as usize..].chars();

               
        for i in 0 .. seq.len() - k as usize{
            let p0 = map_ATCG(exit_char.next().expect("msg"), 0);
            label -= (p0 << ((k-1) << 1));
            label = label << 2;
            label = map_ATCG(enter_char.next().expect("msg"), label);
                     
            let mod_score = (a*label) & max;
            if mod_score as f64 <= Hs{
                kmers.entry(label).or_insert(Vec::new()).push(ind);
                ind += 1;
            }           
        }
    }
    
    (ind+1, kmers)
    
}
