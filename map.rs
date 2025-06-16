use crate::sketch;
use crate::sketch::FracMinHash;
use crate::types::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::i32::MAX;
use std::thread;
use std::sync::{Mutex, Arc, RwLock};
use num::integer::div_ceil;
use std::time::Instant;
use std::cmp;

pub fn create_kmer_map(ref_sketch:&HashMap<String, Sketch>, b:usize ) -> HashMap<i64, HashSet<String>>{

    let mut kmer_map:HashMap<i64, HashSet<String>> = HashMap::new();
    let mut to_remove:HashSet<i64> = HashSet::new();

    let mut count = 0;

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
            //if kmer_map.entry(*kmer).or_insert(HashSet::new()).len() < b{
                //kmer_map.entry(*kmer).or_insert(HashSet::new()).insert(sketch.id.clone());
            //else{
                //to_remove.insert(*kmer);}}
    
    for kmer in to_remove.iter(){
        kmer_map.remove_entry(kmer);
    }

    //println!("Made kmer map");

    kmer_map
}

pub fn find_ref_matches(ref_sketches_shared: Arc<RwLock<HashMap<String, Sketch>>>, sample_chunks_shared:Arc<RwLock<Vec<HashMap<String, String>>>>, b:usize, t: i32, m:i32, k:i32, a:i64, s:f64, check:i32)
 -> (Arc<RwLock<HashMap<String, Vec<Match>>>>, Arc<RwLock<HashMap<String, usize>>>) {

    let mut ref_map:HashMap<String, Vec<Match>> = HashMap::new(); 
    let mut criteria_1:HashMap<String, usize> = HashMap::new();
    let mut handles = vec![];    
    //let mut sample_chunks_shared = Arc::new(RwLock::new(sample_data.clone()));
    //let mut ref_sketches_shared = Arc::new(RwLock::new(ref_sketches));  

    let now = Instant::now();

    let ref_sketches = ref_sketches_shared.read().unwrap();
    //println!("Beginning to make the kmer map");
    let kmer_map = create_kmer_map(&ref_sketches, b);
    drop(ref_sketches);
    let elapsed = now.elapsed();    
    //println!("Elapsed: {:.2?}", elapsed); 
    let mut crit_1_novel = 0;
    let mut crit_1_annot = 0;

    let mut kmer_map_shared = Arc::new(RwLock::new(kmer_map));
    let mut ref_map_shared = Arc::new(RwLock::new(ref_map));
    let mut crit_1_shared = Arc::new(RwLock::new(criteria_1));

    let mut crit_novel_shared = Arc::new(RwLock::new(crit_1_novel));
    let mut crit_annot_shared = Arc::new(RwLock::new(crit_1_annot));

    for i in 0..t as usize{
        let sample_chunks_shared = Arc::clone(&sample_chunks_shared); 
        let kmer_map_shared = Arc::clone(&kmer_map_shared);
        let ref_sketches_shared = Arc::clone(&ref_sketches_shared);
        let ref_map_shared = Arc::clone(&ref_map_shared);
        let crit_1_shared = Arc::clone(&crit_1_shared);

        let mut crit_novel_shared = Arc::clone(&crit_novel_shared);
        let mut crit_annot_shared = Arc::clone(&crit_annot_shared);

        let handle = thread::spawn(move || {
            let guard = sample_chunks_shared.read().unwrap();
            let sample_chunk = (guard).get(i).expect("msg");
            let kmer_map = kmer_map_shared.read().unwrap();      
            let ref_sketches = ref_sketches_shared.read().unwrap();
            
            let mut r_count = 0;

            for (id, seq) in sample_chunk.iter(){

                let (size, sketch) = FracMinHash(&seq, k, a,s);

                let mut opt_similarity = -1.0 * f32::INFINITY;            
                let mut opt_ref_match = "".to_string();
                let mut opt_gap = 0;
                
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
                    if *num_matches == max_matches{
                        let r_sketch = ref_sketches.get(ref_id).expect("msg");
                        let (chain_score, gap, val) = kmer_chain(&(sketch),  &r_sketch.kmers, check);

                        let sim_score = chain_score / (size as f32 + r_sketch.size as f32 - max_matches as f32);
                        if sim_score > opt_similarity{
                            opt_similarity = sim_score;
                            opt_ref_match = ref_id.clone();
                            opt_gap = gap;
                        }
                    }              
                }

                r_count += 1;
                if r_count % 20000 == 0{
                    //println!("Map Input read {}", r_count);
                    //println!("Opt score: {}", opt_similarity);
                }
                
                if opt_gap > m{
                    let mut crit_1 = crit_1_shared.write().unwrap();
                    crit_1.entry(id.clone()).or_insert(i);
                    if id.contains("novel"){
                        *crit_novel_shared.write().unwrap() += 1;
                    }
                    else{
                        *crit_annot_shared.write().unwrap() += 1;
                    }
                }
                else{
                    let mut ref_map = ref_map_shared.write().unwrap(); 
                    ref_map.entry(opt_ref_match.clone()).or_insert(Vec::new()).push(build_Match(id.clone(), (10000.0 * opt_similarity) as i32, opt_gap, i));
                    drop(ref_map)
                }      
            }

        });
        handles.push(handle);
    }

    for i in handles{
        i.join().unwrap();
    } 
    //println!("Finished Joining");
    //let ref_map = ref_map_shared.lock().unwrap().clone();
    //println!("Beginning to pass back");

    //println!("Novel w/ gaps: {}", crit_novel_shared.read().unwrap());
    //println!("Annot w/ gaps: {}", crit_annot_shared.read().unwrap());

    (ref_map_shared, crit_1_shared)
}

pub fn kmer_chain(sample_kmers:&HashMap<i64, Vec<i32>>, ref_kmers:&HashMap<i64, Vec<i32>>, check:i32) -> (f32, i32, f32){
    
    let mut g:Vec<(f32, i32, f32)> = Vec::new();  
    let mut anchors:Vec<(i32, i32)> = Vec::new();

    //println!("Check is {}", check);

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

    //println!("Anchors length: {}", anchors.len());
      
    for i in 0.. anchors.len(){
        let mut optimal_chain = 1.0;
        let mut optimal_gap = 0;
        let mut max_score = 1.0;
        //let mut max_score = k as f32;
        if i == 0{
            
            g.push((optimal_chain, optimal_gap, max_score));
        }
        else{
            //let mut optimal_chain = g[(i-1) as usize].0;
            let mut unique_idx = 0;
            let mut l = 0;

            while unique_idx < check{
            //while unique_idx < 3{
            //while unique_idx < m{
                let j:i32 = i as i32 - l - 1;
                l += 1;

                if j >= 0{
                    let j = j as usize;

                    if anchors[j].0 != anchors[j+1].0{
                        unique_idx += 1;
                    }
                    if (anchors[j].1 < anchors[i].1) && (anchors[j].0 < anchors[i].0){
                        //let mut chain = g[j as usize].0 + 1.0;
                        let gap = ((anchors[j].1 - anchors[i].1) - (anchors[j].0 - anchors[i].0)).abs();
                        let chain = g[j as usize].0 + 1.0;
                        let mut score = g[j as usize].2;
                        //let new_bases = cmp::min(k, cmp::min(anchors[i].1 - anchors[j].1, anchors[i].0 - anchors[j].0)) as f32;
                        //score += new_bases;

                        /*
                         
                        if gap >= m{
                            //score -= cmp::min(0.01 * (k * gap) as f32, (gap as f32).ln() / (2.0_f32).ln());
                            //score -= (gap as f32).ln() / (2.0_f32).ln();
                            score -= gap as f32;
                        }
                        else{
                            //score -= (0.01 * (k * gap) as f32 + 0.5 * (gap as f32).ln() / (2.0_f32).ln());
                            score -= (gap * gap) as f32;
                        }
                        */
                        //if score >= max_score{
                        if chain > optimal_chain{
                            max_score = score;
                            optimal_gap = std::cmp::max(gap, g[j as usize].1 as i32);
                            //optimal_chain += 1.0;
                            optimal_chain = chain;
                        }

                        //println!("{}", gap);
                        //let gap = cmp::max((anchors[j].1 - anchors[i].1).abs(), (anchors[j].0 - anchors[i].0).abs());
                        //println!("{}, {}, {}", (anchors[j].1 - anchors[i].1).abs(), (anchors[j].0 - anchors[i].0).abs(), gap);
                        
                        /*
                        if gap > m {
                            //println!("gap is {}", gap);

                            val -= ((gap) * gap) as f32;
                            
                        }                        
                        if val > max_val{
                            
                            max_val = val;
                            optimal_gap = std::cmp::max(gap, g[j as usize].1 as i32);
                        }
                         */
                        
                    }
                }
                else{
                    unique_idx = check;
                }
            }
            g.push((optimal_chain, optimal_gap, max_score));
        }
    }    
    g.pop().expect("msg")
}
