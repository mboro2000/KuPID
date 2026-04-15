use crate::sketch;
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

pub fn create_kmer_map(ref_sketch:&HashMap<i32, Sketch>, b:usize ) -> HashMap<i64, Vec<i32>>{

    let mut kmer_map:HashMap<i64, Vec<i32>> = HashMap::new();
    let mut to_remove:HashSet<i64> = HashSet::new();

    for (id, sketch) in ref_sketch.iter(){
        for kmer in sketch.kmers.keys(){
            if kmer_map.entry(*kmer).or_insert(Vec::new()).len() >= b{
                to_remove.insert(*kmer);
            }
            else{
                kmer_map.entry(*kmer).or_insert(Vec::new()).push(id.clone());
            }
        }
    }
    for kmer in to_remove.iter(){
        kmer_map.remove_entry(kmer);
    }
    kmer_map
}

pub fn compare_exons(ref_seq:&str, seq:&str, unmatched:(i32, i32, i32), max:i64, k:i32, e:f32, end:&String, novel:&mut bool) -> bool{
    let r_unmatch = sketch::get_unmatched(ref_seq, k, end, unmatched.1, max, unmatched.0);
    let q_unmatch = sketch::get_unmatched(seq, k, end, unmatched.2, max, unmatched.0);

    let inter = q_unmatch.intersection(&r_unmatch);
    if (inter.count() as f32) < q_unmatch.len() as f32 * (1.0 - k as f32 * e){
        *novel = true;
    }
    *novel
}

pub fn find_ref_matches(output:String, ref_sketches_shared: Arc<RwLock<HashMap<i32, Sketch>>>, ref_data_shared:Arc<RwLock<Vec<HashMap<i32, String>>>>, ref_chunks_shared:Arc<RwLock<HashMap<i32, usize>>>, sample_data_shared:Arc<RwLock<Vec<HashMap<i32, String>>>>,
     b:usize, t: i32, n:i32, k:i32, s:f64, m:i32, e:f32, B:f32, c:f64, mode:String, l:usize, transcript_names:HashMap<i32, String>, buffer:i32, d:f32)
 -> Arc<RwLock<HashMap<i32, usize>>> {
    let ref_map:HashMap<i32, (usize, Vec<Match>)> = HashMap::new(); 
    let selected:HashMap<i32, usize> = HashMap::new();
    let mut handles = vec![];   
    let as_reads = 0;
    let novel_exon_reads = 0;
    let max:i64 = (1 << (2*k)) - 1;
    let hs = (max as f64 * s) as i64; 

    let ref_sketches = ref_sketches_shared.read().unwrap();
    let kmer_map = create_kmer_map(&ref_sketches, b);
    drop(ref_sketches);

    let kmer_map_shared = Arc::new(RwLock::new(kmer_map));
    let ref_map_shared = Arc::new(RwLock::new(ref_map));
    let selected_shared = Arc::new(RwLock::new(selected));
    let as_reads_shared = Arc::new(RwLock::new(as_reads));
    let novel_exon_shared = Arc::new(RwLock::new(novel_exon_reads));

    for i in 0..t as usize{
        let sample_data_shared = Arc::clone(&sample_data_shared); 
        let kmer_map_shared = Arc::clone(&kmer_map_shared);
        let ref_sketches_shared = Arc::clone(&ref_sketches_shared);
        let ref_data_shared = Arc::clone(&ref_data_shared);
        let ref_chunks_shared = Arc::clone(&ref_chunks_shared);
        let ref_map_shared = Arc::clone(&ref_map_shared);
        let selected_shared = Arc::clone(&selected_shared);
        let as_reads_shared = Arc::clone(&as_reads_shared);
        let novel_exon_shared = Arc::clone(&novel_exon_shared);

        let handle = thread::spawn(move || {

            let guard = sample_data_shared.read().unwrap();
            let sample_chunk = (guard).get(i).expect("msg");
            let kmer_map = kmer_map_shared.read().unwrap();      
            let ref_sketches = ref_sketches_shared.read().unwrap();
            
            for (id, seq) in sample_chunk.iter(){
                
                if seq.len() >= k as usize{
                    let adjusted_len = seq.len() as i32;
                    let rc_sequence = seq.chars()
                        .rev()
                        .map(|c| match c {
                            'A' => 'T',
                            'T' => 'A',
                            'C' => 'G',
                            'G' => 'C',
                            'a' => 't',
                            't' => 'a',
                            'c' => 'g',
                            'g' => 'c',
                            _ => c,
                        }).collect::<String>();

                    let (size, sketch) = sketch::frac_min_hash(&seq, k,max, hs);

                    //// currently, the sketch stores the end point of the kmer as the position
                    
                    let mut opt_chain = init_chain();             
                    let mut forward_matches:HashMap<i32, i32> = HashMap::new();
                    let mut rc_matches:HashMap<i32, i32> = HashMap::new();
                    let mut max_matches = 0;
                    let mut forward_max = 0;
                    let mut rc_max = 0;
                    let default:Vec<i32> = Vec::new();

                    let best_size: i32;
                    let best_seq: &str;
                    let best_sketch: HashMap<i64, Vec<i32>>;
                    let ref_matches: HashMap<i32, i32>;

                    for kmer in sketch.keys(){
                        let genes = kmer_map.get(kmer).unwrap_or_else(|| &default);
                        for ref_id in genes.iter(){
                            *forward_matches.entry(ref_id.clone()).or_insert(0) += 1;
                        }
                    }
                    for mt in forward_matches.values(){
                        if *mt > forward_max{
                            forward_max = *mt;
                        }
                    }
                            
                    if max_matches < ((seq.len() - k as usize) as f32 * (1.0 - k as f32 * e) * s as f32) as i32{    
                        let (rc_size, rc_sketch) = sketch::frac_min_hash(&rc_sequence, k,max, hs);

                        for kmer in rc_sketch.keys(){
                            let genes = kmer_map.get(kmer).unwrap_or_else(|| &default);
                            for ref_id in genes.iter(){
                                *rc_matches.entry(ref_id.clone()).or_insert(0) += 1;
                            }
                        }
                        for mt in rc_matches.values(){
                            if *mt > rc_max{
                                rc_max = *mt;
                            }
                        }

                        if rc_max > forward_max{
                            ref_matches = rc_matches;
                            best_sketch = rc_sketch;
                            best_size = rc_size;
                            best_seq = &rc_sequence;
                            max_matches = rc_max;
                        }
                        else{
                            ref_matches = forward_matches;
                            best_sketch = sketch;
                            best_size = size;
                            best_seq = seq;
                            max_matches = forward_max;
                        }
                    }
                    else{
                        ref_matches = forward_matches;
                        best_sketch = sketch;
                        best_size = size;
                        best_seq = seq;
                        max_matches = forward_max;
                    }
                    
                    for (ref_id, num_matches) in ref_matches.iter(){
                        if (*num_matches == max_matches) && (max_matches > 1){
                            let r_sketch = ref_sketches.get(ref_id).expect("msg");
                            let (chain_score, gap, (first_pos_q, first_pos_r), (last_pos_q, last_pos_r), num_repeat) = chain::kmer_chain_unique(&(best_sketch),  &r_sketch.kmers, m);
                            let overhang_5p = first_pos_q - first_pos_r;
                            let overhang_3p = (adjusted_len - last_pos_q) - (r_sketch.seq_len as i32 - last_pos_r);
                            
                            let unmatched_5p:(i32,i32,i32) = (first_pos_q-k, first_pos_r, first_pos_q);
                            let unmatched_3p:(i32,i32,i32) = (adjusted_len - last_pos_q, last_pos_r, last_pos_q);

                            let sim_score = chain_score as f32 / ((best_size + r_sketch.size) as f32 - chain_score as f32 - num_repeat as f32);
                            if sim_score > opt_chain.similarity{
                                opt_chain = build_chain(chain_score, ref_id.clone(), gap, overhang_3p, overhang_5p, unmatched_5p, unmatched_3p, sim_score);
                            }
                        }              
                    }

                    if opt_chain.score > 0{

                        if opt_chain.gap >= n{                    
                            let mut selected = selected_shared.write().unwrap();
                            selected.entry(id.clone()).or_insert(i);
                            *as_reads_shared.write().unwrap() += 1;
                        }
                        
                        if std::cmp::max(opt_chain.query_overhang_3p, opt_chain.query_overhang_5p) >= n + buffer{
                            
                            let mut selected = selected_shared.write().unwrap();
                            selected.entry(id.clone()).or_insert(i);
                            *novel_exon_shared.write().unwrap() += 1;       
                        }
                        else if (opt_chain.unmatched_5p.0 >= n + buffer) || (opt_chain.unmatched_3p.0 >= n + buffer){
                            let ref_chunk = ref_chunks_shared.read().unwrap();
                            let ref_data = ref_data_shared.read().unwrap();
                            let ref_seq = ref_data[*ref_chunk.get(&opt_chain.ref_match).expect("msg")].get(&opt_chain.ref_match).expect("msg");

                            let mut novel_exon = false;

                            if (opt_chain.query_overhang_5p < 0) && (opt_chain.unmatched_5p.0 >= n + buffer){
                                novel_exon = compare_exons(ref_seq, best_seq, opt_chain.unmatched_5p, max, k, e, &"5'".to_string(), &mut novel_exon);
                            }
                            if (opt_chain.query_overhang_3p < 0) && (opt_chain.unmatched_3p.0 >= n + buffer){
                                novel_exon = compare_exons(ref_seq, best_seq, opt_chain.unmatched_3p, max, k, e, &"3'".to_string(), &mut novel_exon);
                            }

                            if novel_exon{
                                let mut selected = selected_shared.write().unwrap();
                                selected.entry(id.clone()).or_insert(i);
                                *novel_exon_shared.write().unwrap() += 1;
                            }
                            else{
                                let mut ref_map = ref_map_shared.write().unwrap(); 
                                ref_map.entry(opt_chain.ref_match.clone()).or_insert((ref_seq.len(), Vec::new())).1.push(build_match(id.clone(), opt_chain.similarity, i));
                            }  
                        }   
                        else{
                            let ref_chunk = ref_chunks_shared.read().unwrap();
                            let ref_data = ref_data_shared.read().unwrap();
                            let ref_seq = ref_data[*ref_chunk.get(&opt_chain.ref_match).expect("msg")].get(&opt_chain.ref_match).expect("msg");
                            let mut ref_map = ref_map_shared.write().unwrap(); 
                            ref_map.entry(opt_chain.ref_match.clone()).or_insert((ref_seq.len(), Vec::new())).1.push(build_match(id.clone(), opt_chain.similarity, i));
                        }
                    }
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

        for (iso, (_r_len, matches)) in ref_map.iter(){
            if matches.len() > l{
                let mut rng = rand::rng();
                let chosen = matches.choose_multiple(&mut rng, l);
                let scale = matches.len() as f64 / chosen.len() as f64;
                for read in chosen{
                    selected.entry(read.sample_id.clone()).or_insert(read.chunk);
                }               
                scale_outline.push_str(&(iso.to_string()+ "," + &matches.len().to_string() + "," + &scale.to_string() + "\n"));
            }
            else{
                for read in matches.iter(){
                    selected.entry(read.sample_id.clone()).or_insert(read.chunk);
                }
                let scale = matches.len() as f64 / l as f64;
                let ref_id = transcript_names.get(iso).expect("msg");
                scale_outline.push_str(&(ref_id.to_string()+ "," + &matches.len().to_string() + "," + &scale.to_string() + "\n"));
            } 
        }
        fs::write(scale_factor, scale_outline);  
    }

    if mode == "discovery".to_string(){
        let mut selected = selected_shared.write().unwrap();
        let mut avg_sq:Vec<i32> = Vec::new();
        let mut group_avg_sq:HashMap<i32, (i32, HashSet<usize>)> = HashMap::new();
        let sim_ceiling = B * (1.0 - k as f32 *e);

        for iso in ref_map.keys(){
            let mut scores:Vec<(usize, i32)> = Vec::new();
            let mut scores_filtered:Vec<i32> = Vec::new();
            let (r_len, reads) = ref_map.get(iso).expect("msg");
            let exon_sim = (1000.0 * d / *r_len as f32) as i32;
            let num_reads = reads.len();
            let mut to_omit:HashSet<usize> = HashSet::new();
            let mut i = 0;
            for read in reads.iter(){
                scores.push((i, (read.similarity * 1000.0) as i32));     
                i += 1;        
            }   

            scores.sort_by(            |a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
            if scores.len() > 1{
                for i in 0..scores.len(){
                    let score = scores[i];
                    let mut prev = false;
                    let mut next = false;
                    if i > 0{
                        if scores[i-1].1 > (score.1 - exon_sim){
                            prev = true;
                        }
                    }
                    if i < scores.len() - 1{
                        if scores[i+1].1 < score.1 + exon_sim{
                            next = true;
                        }
                    }

                    if prev || next{
                        scores_filtered.push(score.1)
                    }
                    else{
                        to_omit.insert(score.0);
                    }
                }
            }


            let mut score_avg = 0;
            if scores_filtered.len() > 0{
                score_avg = (scores_filtered.iter().sum::<i32>() / scores_filtered.len() as i32) as i32;
            }
            if score_avg >= 0{
                group_avg_sq.entry(iso.clone()).or_insert((score_avg*score_avg, to_omit));
                for _i in 0..num_reads{
                    avg_sq.push(score_avg*score_avg);
                }  
            }
        }

        avg_sq.sort();
        
        let mut cand_to_select = (c as f64  * (*as_reads_shared.read().unwrap() as f64)) as usize;
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
                for (iso, (avg_sq, to_omit)) in group_avg_sq.into_iter(){
                    let matches = &ref_map.get(&iso).expect("msg").1;            
                    if avg_sq <= avg_sq_threshold{
                        let mut i = 0;
                        for read in matches.iter(){
                            if read.similarity < sim_ceiling{        
                                if !to_omit.contains(&i){
                                    selected.entry(read.sample_id.clone()).or_insert(read.chunk);  
                                }                                                                                    
                            }
                            i += 1;
                        }
                    }
                }        
            }
        }
    }

selected_shared
 }
