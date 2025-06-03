use crate::sketch;
use crate::sketch::FracMinHash;
use crate::types;
use crate::types::build_Match;
use crate::types::Sketch;
use std::collections::HashMap;
use std::collections::HashSet;
use std::thread;
use std::sync::{Mutex, Arc, RwLock};
use num::integer::div_ceil;
use std::time::Instant;

pub fn create_kmer_map(ref_sketch:&HashMap<String, types::Sketch>, b:usize ) -> HashMap<i64, HashSet<String>>{

    let mut kmer_map:HashMap<i64, HashSet<String>> = HashMap::new();
    let mut to_remove:HashSet<i64> = HashSet::new();

    for (id, sketch) in ref_sketch.iter(){
        for (kmer, pos) in sketch.kmers.iter(){
            if kmer_map.entry(*kmer).or_insert(HashSet::new()).len() < b{
                kmer_map.entry(*kmer).or_insert(HashSet::new()).insert(sketch.id.clone());
            }
            else{
                to_remove.insert(*kmer);
            }
        }
    }
    
    for kmer in to_remove.iter(){
        kmer_map.remove_entry(kmer);
    }

    println!("Made kmer map");

    kmer_map
}

pub fn find_ref_matches(ref_sketches:HashMap<String, types::Sketch>, sample_data:&Vec<HashMap<String, String>>, b:usize, t: i32, m:i32, k:i32, a:i64, s:f64) -> HashMap<String, Vec<types::Match>> {

    let now = Instant::now();

    let kmer_map = create_kmer_map(&ref_sketches, b);
    let elapsed = now.elapsed();    
    println!("Elapsed: {:.2?}", elapsed); 

    let mut ref_map:HashMap<String, Vec<types::Match>> = HashMap::new(); 

    let mut handles = vec![];    
    let mut sample_chunks_shared = Arc::new(RwLock::new(sample_data.clone()));
    let mut ref_sketches_shared = Arc::new(RwLock::new(ref_sketches));
    let mut kmer_map_shared = Arc::new(RwLock::new(kmer_map));
    let mut ref_map_shared = Arc::new(Mutex::new(ref_map));

    for i in 0..t as usize{
        let sample_chunks_shared = Arc::clone(&sample_chunks_shared); 
        let kmer_map_shared = Arc::clone(&kmer_map_shared);
        let ref_sketches_shared = Arc::clone(&ref_sketches_shared);
        let ref_map_shared = Arc::clone(&ref_map_shared);

        let handle = thread::spawn(move || {
            let guard = sample_chunks_shared.read().unwrap();
            let sample_chunk = (guard).get(i).expect("msg");
            let kmer_map = kmer_map_shared.read().unwrap();      
            let ref_sketches = ref_sketches_shared.read().unwrap();


            let mut r_count = 0;

            for (id, seq) in sample_chunk.iter(){

                let (size, sketch) = FracMinHash(&seq, k, a,s);

                r_count += 1;
                if r_count % 1000 == 0{
                    println!("Map Input read {}", r_count);
                }

                let mut opt_similarity = -1.0 * f32::INFINITY;            
                let mut opt_ref_match = "".to_string();
                let mut opt_gap = 0;
                
                let mut ref_matches:HashMap<String, i32> = HashMap::new();
                let mut max_matches = 0;

                for kmer in sketch.keys(){
                    if kmer_map.contains_key(kmer){
                        for ref_id in kmer_map.get(kmer).expect("msg"){
                            *ref_matches.entry(ref_id.clone()).or_insert(0) += 1;
                            let m = *ref_matches.entry(ref_id.clone()).or_insert(0);
                            if  m > max_matches{
                                max_matches = m;
                            }
                        }
                    }
                }

                for (ref_id, num_matches) in ref_matches.iter(){
                    if *num_matches == max_matches{
                        let r_sketch = ref_sketches.get(ref_id).expect("msg");
                        let (chain_score, gap) = kmer_chain(&(sketch), m, &r_sketch.kmers);

                        let sim_score = chain_score / (size as f32 + r_sketch.size as f32 - max_matches as f32);
                        if sim_score > opt_similarity{
                            opt_similarity = sim_score;
                            opt_ref_match = ref_id.clone();
                            opt_gap = gap;
                        }
                    }
                    let mut ref_map = ref_map_shared.lock().unwrap(); 
                    ref_map.entry(opt_ref_match.clone()).or_insert(Vec::new()).push(build_Match(id.clone(), opt_similarity, opt_gap, i));
                    drop(ref_map)
                }
            }

        });
        handles.push(handle);
    }

    for i in handles{
        i.join().unwrap();
    } 
    let ref_map = ref_map_shared.lock().unwrap().clone();
    ref_map
}

pub fn kmer_chain(sample_kmers:&HashMap<i64, Vec<i32>>, m:i32, ref_kmers:&HashMap<i64, Vec<i32>>) -> (f32, i32){
    
    let mut g:Vec<(f32, i32)> = Vec::new();  
    let mut anchors:Vec<(i32, i32)> = Vec::new();

    for (kmer, pos) in sample_kmers{
        if ref_kmers.contains_key(kmer){
            for i in pos{
                for j in ref_kmers.get(kmer).expect("msg"){
                    anchors.push((*i, *j));
                }
            }            
        }
    }       
    anchors.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
      
    for i in 0.. anchors.len(){
        let mut optimal_gap = 0;
        if i == 0{
            g.push((1.0, 0));
        }
        else{
            let mut max_val = -1.0 * std::f32::INFINITY;
            let mut unique_idx = 0;
            let mut l = 0;

            while unique_idx < m{
                let j:i32 = i as i32 - l - 1;
                l += 1;

                if j >= 0{
                    let j = j as usize;
                    if anchors[j].0 != anchors[j+1].0{
                        unique_idx += 1;
                    }
                    if (anchors[j].1 != anchors[i].1) && (anchors[j].0 != anchors[i].0){
                        let mut val = g[j as usize].0 + 1.0;
                        let gap = ((anchors[j].1 - anchors[i].1) - (anchors[j].0 - anchors[i].0)).abs();
                        if gap > m {
                            val -= ((gap) * gap) as f32;
                        }                        
                        if val > max_val{
                            max_val = val;
                            optimal_gap = std::cmp::max(gap, g[j as usize].1 as i32);
                        }
                    }
                }
                else{
                    unique_idx = m;
                }
            }
            g.push((max_val, optimal_gap));
        }
    }    
    g.pop().expect("msg")
}
