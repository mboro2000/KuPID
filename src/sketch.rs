
use std::collections::HashMap;
use bio::io::fasta;
use std::str;
use std::cmp;
use crate::types::*;
use std::thread;
use std::sync::{Arc, RwLock};

//Read in and process input files for multi-threading
pub fn read_input(t:i32, file:&String) -> Vec<HashMap<String, String>>{
    let mut seq_data:Vec<HashMap<String, String>> = Vec::new();    
    let read_results = fasta::Reader::from_file(file);
    for _i in 0 .. t{
        let read_chunk:HashMap<String, String> = HashMap::new();
        seq_data.push(read_chunk);
    }
    let mut r = 0;
    for reader in read_results{
        for result in reader.records(){
            r += 1;
            let record = result.expect("Error during fasta record parsing. Please submit a viable fasta file");
            let seq = str::from_utf8(record.seq()).unwrap().to_string();
            seq_data[(r % t) as usize].entry(record.id().to_string()).or_insert(seq);
        }
    }  
    seq_data
}

//Generate the kmer sketches for a set of nucleotide strings
pub fn get_sketches(read_data:Vec<HashMap<String, String>>, k:i32, s:f64, t:i32) -> Arc<RwLock<HashMap<String, Sketch>>> {   
    let sketches:HashMap<String, Sketch> = HashMap::new();
    let mut handles = vec![];    
    let sketches_shared = Arc::new(RwLock::new(sketches));
    let read_data_shared = Arc::new(RwLock::new(read_data));
    for i in 0..t{
        let sketches_shared = Arc::clone(&sketches_shared); 
        let read_data_shared = Arc::clone(&read_data_shared); 
        let handle = thread::spawn(move || {
            let guard = read_data_shared.read().unwrap();
            let seq_chunk = (guard).get(i as usize).expect("msg");
            for (id, seq) in seq_chunk{       
                let (num_kmers, sketch) = frac_min_hash(&seq, k,s);
                let mut sketches = sketches_shared.write().unwrap();
                sketches.entry(id.clone()).or_insert(build_sketch(num_kmers, sketch, seq.len()));
                drop(sketches);
            }
        });
        handles.push(handle);
    }
    for i in handles{
        i.join().unwrap();
    } 
    sketches_shared
}

//Encodes nucleotides as numerical values
pub fn map_atcg(item:char, mut label:i64) -> i64{
    if 'A' == item || 'a' == item{label += 0;}
    else if 'C' == item || 'c' == item {label += 1;}
    else if 'G' == item || 'g' == item {label += 2;}
    else{label += 3;}
    label
}

//Applies invertible hash function to the nucleotide encodings
pub fn integer_hash(r:i64, m:i128) -> i64{
    let mut x = r as i128;
    x = (!x + (x<<21)) & m;
    x = x ^ x>>24;
    x = (x + (x<<3) + (x<<8)) & m;
    x = x ^ x>>14;
    x = (x+(x<<2) + (x<<4)) & m;
    x = x^x>>28;
    x=(x+(x<<31)) & m;
    x as i64
}

// Returns a table of minimizers selected using the frac_min_hash method
// Returns: total:i32 - number of kmers that were included in the sketch
//          sketch:HashMap<i64, Vec<i32>> - <numerical hash value of kmer: starting positions of the kmer in the sequence>
pub fn frac_min_hash(seq:&str, k:i32, s:f64) -> (i32, HashMap<i64, Vec<i32>>) {
    let max:i64 = (1 << (2*k)) - 1;
    let hs:f64 = max as f64 * s;
    let mut label:i64 = 0;
    let mut sketch:HashMap<i64, Vec<i32>> = HashMap::new();
    let mut total = 0;
    let s = &seq[0 .. cmp::min(k as usize, seq.len())];

    for item in s.chars() {
        label  = label << 2;
        label = map_atcg(item, label);            
    }  

    if seq.len() > k as usize{
        let mod_score = integer_hash(label, max as i128);    
        if mod_score as f64 <= hs{
            sketch.entry(mod_score).or_insert(Vec::new()).push(0);
            total += 1;
        }
        let mut exit_char = seq[ .. (seq.len() as i32 - k) as usize].chars();
        let mut enter_char = seq[k as usize..].chars();                              
        for i in 0 .. seq.len() - k as usize{
            let p0 = map_atcg(exit_char.next().expect("msg"), 0);
            label -= p0 << ((k-1) << 1);
            label = label << 2;
            label = map_atcg(enter_char.next().expect("msg"), label);                     
            let mod_score = integer_hash(label, max as i128);
            if mod_score as f64 <= hs{
                sketch.entry(mod_score).or_insert(Vec::new()).push(i as i32 + 1);
                total += 1;
            }           
        }
    }
    (total, sketch)    
}