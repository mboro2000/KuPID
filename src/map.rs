use crate::sketch::frac_min_hash;
use crate::types::*;
use crate::chain;
use std::collections::HashMap;
use std::collections::HashSet;
use std::thread;
use std::sync::{Arc, RwLock};
use std::path::Path;
use std::fs;
use rand;
use rand::seq::IndexedRandom;

pub fn create_kmer_map(ref_sketch:&HashMap<String, Sketch>, b:usize ) -> HashMap<i64, HashSet<String>>{
    let mut kmer_map:HashMap<i64, HashSet<String>> = HashMap::new();
    let mut to_remove:HashSet<i64> = HashSet::new();

    for (id, sketch) in ref_sketch.iter(){
        for kmer in sketch.kmers.keys(){
            if kmer_map.entry(*kmer).or_insert(HashSet::new()).len() >= b{
                to_remove.insert(*kmer);
            }
            else{
                kmer_map.entry(*kmer).or_insert(HashSet::new()).insert(id.clone());
            }
        }
    }
    for kmer in to_remove.iter(){
        kmer_map.remove_entry(kmer);
    }
    kmer_map
}

pub fn find_ref_matches(output:String, ref_sketches_shared: Arc<RwLock<HashMap<String, Sketch>>>, sample_chunks_shared:Arc<RwLock<Vec<HashMap<String, String>>>>, b:usize, t: i32, n:i32, k:i32, s:f64, m:i32, e:f32, B:f32, c:f64, mode:String, l:usize, g:i32)
 -> (Arc<RwLock<HashMap<String, usize>>>) {
    let ref_map:HashMap<String, Vec<Match>> = HashMap::new(); 
    let selected:HashMap<String, usize> = HashMap::new();
    let mut handles = vec![];   
    let as_reads = 0;
    let novel_exon_reads = 0;

    let ref_sketches = ref_sketches_shared.read().unwrap();
    let kmer_map = create_kmer_map(&ref_sketches, b);
    drop(ref_sketches);

    let kmer_map_shared = Arc::new(RwLock::new(kmer_map));
    let ref_map_shared = Arc::new(RwLock::new(ref_map));
    let selected_shared = Arc::new(RwLock::new(selected));
    let as_reads_shared = Arc::new(RwLock::new(as_reads));
    let novel_exon_shared = Arc::new(RwLock::new(novel_exon_reads));

    for i in 0..t as usize{
        let sample_chunks_shared = Arc::clone(&sample_chunks_shared); 
        let kmer_map_shared = Arc::clone(&kmer_map_shared);
        let ref_sketches_shared = Arc::clone(&ref_sketches_shared);
        let ref_map_shared = Arc::clone(&ref_map_shared);
        let selected_shared = Arc::clone(&selected_shared);
        let as_reads_shared = Arc::clone(&as_reads_shared);
        let novel_exon_shared = Arc::clone(&novel_exon_shared);
        
        let handle = thread::spawn(move || {
            let guard = sample_chunks_shared.read().unwrap();
            let sample_chunk = (guard).get(i).expect("msg");
            let kmer_map = kmer_map_shared.read().unwrap();      
            let ref_sketches = ref_sketches_shared.read().unwrap();
            
            for (id, seq) in sample_chunk.iter(){
                let mut tail_len = 0;
                let mut chars = seq.chars().rev();
                let mut c = chars.next().expect("msg");
                while c=='A' || c=='a'{
                    tail_len += 1;
                    c = chars.next().expect("msg");
                }
                let adjusted_len = seq.len() as i32 - tail_len;
                let (size, sketch) = frac_min_hash(&seq, k,s);
                
                let mut opt_chain = init_chain();                
                let mut ref_matches:HashMap<String, i32> = HashMap::new();
                let mut max_matches = 0;

                for kmer in sketch.keys(){
                    if kmer_map.contains_key(kmer){
                        for ref_id in kmer_map.get(kmer).expect("msg"){
                            *ref_matches.entry(ref_id.clone()).or_insert(0) += 1;
                            let mt = *ref_matches.entry(ref_id.clone()).or_insert(0);
                            if  mt > max_matches{
                                max_matches = mt;
                            }
                        }
                    }
                }

                for (ref_id, num_matches) in ref_matches.iter(){
                    if (*num_matches == max_matches) && (max_matches > 0){
                        let r_sketch = ref_sketches.get(ref_id).expect("msg");
                        let (chain_score, gap, (first_pos_q, first_pos_r), (last_pos_q, last_pos_r)) = chain::kmer_chain(&(sketch),  &r_sketch.kmers, m);
                        let mut ref_gap = 0;
                        let gap_5p = std::cmp::max(0,first_pos_q - first_pos_r);
                        let gap_3p = std::cmp::max(0,(adjusted_len - last_pos_q) - (r_sketch.seq_len as i32 - last_pos_r));
                        if std::cmp::max(first_pos_r - first_pos_q, (r_sketch.seq_len as i32 - last_pos_r) - (adjusted_len - last_pos_q)) > 0{
                            ref_gap = std::cmp::max(first_pos_q, adjusted_len - last_pos_q - k + 1);
                        }
                        
                        let sim_score = chain_score as f32 / ((size + r_sketch.size) as f32 - chain_score as f32);
                        if sim_score > opt_chain.similarity{
                            opt_chain = build_chain(chain_score, ref_id.clone(), gap, std::cmp::max(gap_3p, gap_5p), ref_gap, sim_score);
                        }
                    }              
                }

                if opt_chain.gap > n{                    
                    let mut selected = selected_shared.write().unwrap();
                    selected.entry(id.clone()).or_insert(i);
                    *as_reads_shared.write().unwrap() += 1;
                }
                else if  opt_chain.query_gap > n{
                    let mut selected = selected_shared.write().unwrap();
                    selected.entry(id.clone()).or_insert(i);
                    *novel_exon_shared.write().unwrap() += 1;          
                }
                else if opt_chain.ref_gap > g{
                    let mut selected = selected_shared.write().unwrap();
                    selected.entry(id.clone()).or_insert(i);
                    *novel_exon_shared.write().unwrap() += 1; 
                }             
                else{
                    let mut ref_map = ref_map_shared.write().unwrap(); 
                    ref_map.entry(opt_chain.ref_match.clone()).or_insert(Vec::new()).push(build_match(id.clone(), opt_chain.similarity, i));
                }      
            }
        });
        handles.push(handle);
    }

    for i in handles{
        i.join().unwrap();
    } 

    let ref_map = ref_map_shared.write().unwrap();
    if mode == "quantify".to_string(){
        let mut selected = selected_shared.write().unwrap();
        let scale_file = output.clone() + ".scale_factors.csv";
        let scale_factor = Path::new(&scale_file);
        let mut scale_outline = "".to_string();
        scale_outline.push_str("Transcript,Group Count,Scale\n");

        for (iso, matches) in ref_map.iter(){
            if matches.len() > l{
                let mut rng = rand::rng();
                let chosen = matches.choose_multiple(&mut rng, l);
                let scale = matches.len() as f64 / chosen.len() as f64;
                for read in chosen{
                    selected.entry(read.sample_id.clone()).or_insert(read.chunk);
                }               
                scale_outline.push_str(&(iso.clone()+ "," + &matches.len().to_string() + "," + &scale.to_string() + "\n"));
            }
            else{
                for read in matches.iter(){
                    selected.entry(read.sample_id.clone()).or_insert(read.chunk);
                }
                let scale = matches.len() as f64 / l as f64;
                scale_outline.push_str(&(iso.clone()+ "," + &matches.len().to_string() + "," + &scale.to_string() + "\n"));
            } 
        }
        fs::write(scale_factor, scale_outline);  
    }

    if mode == "discovery".to_string(){
        let mut selected = selected_shared.write().unwrap();
        let mut avg_sq:Vec<i32> = Vec::new();
        let mut group_avg_sq:HashMap<String, i32> = HashMap::new();
        let sim_ceiling = B * (1.0 - k as f32 *e);

        for iso in ref_map.keys(){
            let mut scores:Vec<f32> = Vec::new();
            let num_reads = ref_map.get(iso).expect("msg").len();
            for read in ref_map.get(iso).expect("msg").iter(){
                scores.push(read.similarity);             
            }   
            let score_avg = (1000.0 * scores.iter().sum::<f32>() / scores.len() as f32) as i32;
            if score_avg >= 0{
                group_avg_sq.entry(iso.clone()).or_insert(score_avg*score_avg);
                for _i in 0..num_reads{
                    avg_sq.push(score_avg*score_avg);
                }  
            }
        }

        avg_sq.sort();
        let mut cand_to_select = (*as_reads_shared.read().unwrap() as f64 * c) as usize;
        if (*novel_exon_shared.read().unwrap()) > cand_to_select{
            cand_to_select = 0;
        }
        else{
            cand_to_select -= (*novel_exon_shared.read().unwrap()) as usize;
        }

        if cand_to_select > 0{
            if avg_sq.len() > 0{
                let mut avg_sq_threshold = avg_sq.pop().expect("msg");
                if cand_to_select < avg_sq.len(){
                    avg_sq_threshold = avg_sq[cand_to_select];
                }
                for (iso, avg_sq) in group_avg_sq.into_iter(){
                    let matches = ref_map.get(&iso).expect("msg");            
                    if avg_sq <= avg_sq_threshold{
                        for read in matches.iter(){
                            if read.similarity < sim_ceiling{        
                                selected.entry(read.sample_id.clone()).or_insert(read.chunk);                                                                                           
                            }
                        }
                    }
                }        
            }
        }
    }
selected_shared
 }
