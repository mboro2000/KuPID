use std::collections::HashMap;
use std::collections::HashSet;
use std::process;
use bio::io::fasta;
use std::str;
use crate::types::*;
use std::thread;
use std::sync::{Arc, RwLock};

//Read in and process input files for multi-threading
pub fn read_input(t:i32, file:&String, encode:bool) -> (Vec<HashMap<i32, String>>, HashMap<i32, String>, HashMap<i32, usize>){
    let mut seq_data:Vec<HashMap<i32, String>> = Vec::new();    
    
    let mut reader = match fasta::Reader::from_file(file) {
        Ok(r) => r,
        Err(e) => {
            println!("File not found: {}", file);
            std::process::exit(1);
        }
    };

    let mut read_encoding:HashMap<i32, String> = HashMap::new();
    let mut chunk_assign:HashMap<i32, usize> = HashMap::new();

    for _i in 0 .. t{
        let read_chunk:HashMap<i32, String> = HashMap::new();
        seq_data.push(read_chunk);
    }
    let mut r = 0;
    
    for result in reader.records(){
        r += 1;
        let record: fasta::Record = result.expect("Error during fasta record parsing. Please submit a viable fasta file");
        let seq = str::from_utf8(record.seq()).unwrap().to_string();
        seq_data[(r % t) as usize].entry(r.clone()).or_insert(seq);
        if encode{
            read_encoding.entry(r.clone()).or_insert(record.id().to_string());
        }
        chunk_assign.entry(r.clone()).or_insert((r % t) as usize);
    }  
    (seq_data, read_encoding, chunk_assign)
}

//Generate the kmer sketches for a set of nucleotide strings
pub fn get_sketches(read_data_shared:Arc<RwLock<Vec<HashMap<i32, String>>>>, k:i32, s:f64, t:i32) -> Arc<RwLock<HashMap<i32, Sketch>>> {  
    let max:i64 = (1 << (2*k)) - 1;
    let hs = (max as f64 * s) as i64; 
    let sketches:HashMap<i32, Sketch> = HashMap::new();
    let mut handles = vec![];    
    let sketches_shared = Arc::new(RwLock::new(sketches));
    //let read_data_shared = Arc::new(RwLock::new(read_data));
    for i in 0..t{
        let sketches_shared = Arc::clone(&sketches_shared); 
        let read_data_shared = Arc::clone(&read_data_shared); 
        let handle = thread::spawn(move || {
            let guard = read_data_shared.read().unwrap();
            let seq_chunk = (guard).get(i as usize).expect("msg");
            for (id, seq) in seq_chunk{       
                let (num_kmers, sketch) = frac_min_hash(&seq, k,max, hs);
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
pub fn integer_hash(r:i64, m:i64) -> i64{
    let mut x = r;
    x = (!x.wrapping_add(x<<21)) & m;
    x = x ^ x>>24;
    x = (x.wrapping_add(x<<3).wrapping_add(x<<8)) & m;
    x = x ^ x>>14;
    x = (x.wrapping_add(x<<2).wrapping_add(x<<4)) & m;
    x = x^x>>28;
    x=(x.wrapping_add(x<<31)) & m;
    x
}

pub fn get_unmatched(seq:&str, k:i32, end:&String, pos:i32, max:i64, l:i32) -> HashSet<i64>{
    let mut label:i64 = 0;
    let mut i = 0;
    let mut kmers:HashSet<i64> = HashSet::new();

    if *end == "5'"{
        for item in seq.chars(){
            label  = label << 2;
            label = map_atcg(item, label);
            i += 1;

            if i > pos{
                break;
            }
            if i >= k{
                label = label & (max >> 2);
                let mod_score = integer_hash(label, max);
                kmers.insert(mod_score);

                if i >= pos - l{
                    let mod_score = integer_hash(label, max);
                    kmers.insert(mod_score);
                }
            }
        }
    }

    if *end == "3'"{
        for item in seq.chars(){
            label  = label << 2;
            label = map_atcg(item, label);
            i += 1;

            
            if i > pos + l{
                break;
            }
            
            if i >= k{
                label = label & (max >> 2);

                if i > pos{
                    let mod_score = integer_hash(label, max);
                    kmers.insert(mod_score);
                }
            }

        }
    }
    kmers
}


// Returns a table of minimizers selected using the frac_min_hash method
// Returns: total:i32 - number of kmers that were included in the sketch
//          sketch:HashMap<i64, Vec<i32>> - <numerical hash value of kmer: starting positions of the kmer in the sequence>
pub fn frac_min_hash(seq:&str, k:i32, max:i64, hs:i64) -> (i32, HashMap<i64, Vec<i32>>) {
    let mut label:i64 = 0;
    let mut sketch:HashMap<i64, Vec<i32>> = HashMap::new();
    let mut total = 0;
    let mut i = 0;

    for item in seq.chars(){
        label  = label << 2;
        label = map_atcg(item, label);
        i += 1;   

        if i >= k{
            let mod_score = integer_hash(label, max);
            if mod_score <= hs{
                sketch.entry(mod_score).or_insert(Vec::new()).push(i as i32);
                total += 1;
            }
            label = label & (max >> 2);
        }
    }
    (total, sketch)    
}
