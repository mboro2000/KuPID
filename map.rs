use crate::sketch::frac_min_hash;
use crate::types::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::thread;
use std::sync::{Arc, RwLock};
use std::time::Instant;
use std::path::Path;
use std::fs;
use std::cmp;

pub fn create_kmer_map(ref_sketch:&HashMap<String, Sketch>, b:usize ) -> HashMap<i64, HashSet<String>>{

    let mut kmer_map:HashMap<i64, HashSet<String>> = HashMap::new();
    let mut to_remove:HashSet<i64> = HashSet::new();

    for (id, sketch) in ref_sketch.iter(){
        for kmer in sketch.kmers.keys(){
            if kmer_map.entry(*kmer).or_insert(HashSet::new()).len() < b{
                kmer_map.entry(*kmer).or_insert(HashSet::new()).insert(id.clone());
            }
            else{
                to_remove.insert(*kmer);
            }
        }
    }
    
    for kmer in to_remove.iter(){
        kmer_map.remove_entry(kmer);
    }

    kmer_map
}

pub fn find_ref_matches(output:String, ref_sketches_shared: Arc<RwLock<HashMap<String, Sketch>>>, sample_chunks_shared:Arc<RwLock<Vec<HashMap<String, String>>>>, b:usize, t: i32, n:i32, k:i32, a:i64, s:f64, m:i32, total_reads:i32, e:f32, B:f32, d:f32, c:f64)
 -> (Arc<RwLock<HashMap<String, usize>>>, HashMap<String, usize>, i32, i32, i32, i32, i32, i32, i32, i32, HashMap<String, (String, usize)>) {

    let ref_map:HashMap<String, Vec<Match>> = HashMap::new(); 
    let AS:HashMap<String, usize> = HashMap::new();

    let mut large_diff:HashMap<String, (String, usize)> = HashMap::new();

    let ref_sketches = ref_sketches_shared.read().unwrap();
    let kmer_map = create_kmer_map(&ref_sketches, b);
    drop(ref_sketches);
    let mut handles = vec![];   
    let as_novel = 0;
    let as_annot = 0;

    let kmer_map_shared = Arc::new(RwLock::new(kmer_map));
    let ref_map_shared = Arc::new(RwLock::new(ref_map));
    let as_shared = Arc::new(RwLock::new(AS));

    let as_novel_shared = Arc::new(RwLock::new(as_novel));
    let as_annot_shared = Arc::new(RwLock::new(as_annot));

    for i in 0..t as usize{
        let sample_chunks_shared = Arc::clone(&sample_chunks_shared); 
        let kmer_map_shared = Arc::clone(&kmer_map_shared);
        let ref_sketches_shared = Arc::clone(&ref_sketches_shared);
        let ref_map_shared = Arc::clone(&ref_map_shared);
        let as_shared = Arc::clone(&as_shared);

        let as_novel_shared = Arc::clone(&as_novel_shared);
        let as_annot_shared = Arc::clone(&as_annot_shared);

        let handle = thread::spawn(move || {
            let guard = sample_chunks_shared.read().unwrap();
            let sample_chunk = (guard).get(i).expect("msg");
            let kmer_map = kmer_map_shared.read().unwrap();      
            let ref_sketches = ref_sketches_shared.read().unwrap();
            
            let r_count = 0;

            for (id, seq) in sample_chunk.iter(){

                let (size, sketch) = frac_min_hash(&seq, k, a,s);

                let mut opt_similarity = -1.0 * f32::INFINITY;            
                let mut opt_ref_match = "".to_string();
                let mut opt_ref_len = 0;
                let mut opt_gap = 0;
                let mut opt_gap_5 = 0;
                let mut opt_gap_3 = 0;
                
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
                        let (chain_score, gap, first_pos_q, first_pos_r, last_pos_q, last_pos_r) = kmer_chain(&(sketch),  &r_sketch.kmers, m);

                        let gap_5 = std::cmp(first_pos_q, first_pos_r);
                        let gap_3 = std::cmp(seq.len() as i32 - last_pos_q, opt_ref_len as i32 - last_pos_r);

                        let sim_score = chain_score / ((size + r_sketch.size) as f32 - chain_score);
                        if sim_score > opt_similarity{
                            opt_similarity = sim_score;
                            opt_ref_match = ref_id.clone();
                            opt_ref_len = r_sketch.seq_len;
                            opt_gap = gap;
                        }
                    }              
                }
                
                if opt_gap > n{
                    let mut AS = as_shared.write().unwrap();
                    AS.entry(id.clone()).or_insert(i);
                    if id.contains("novel"){
                        *as_novel_shared.write().unwrap() += 1;
                    }
                    else{
                        *as_annot_shared.write().unwrap() += 1;
                    }
                }
                else{
                    let mut ref_map = ref_map_shared.write().unwrap(); 
                    let cov_q = seq.len() as f32 / opt_ref_len as f32;
                    let cov_r = opt_ref_len as f32 / seq.len() as f32;
                    let mut cov = 0.0;
                    if cov_q < cov_r{
                        cov = cov_q;
                    }
                    else{
                        cov = cov_r;
                    }
                    ref_map.entry(opt_ref_match.clone()).or_insert(Vec::new()).push(build_match(id.clone(), opt_similarity, i, cov));
                    drop(ref_map)
                }      
            }

        });
        handles.push(handle);
    }

    for i in handles{
        i.join().unwrap();
    } 

    let mut ref_map = ref_map_shared.write().unwrap();

    println!("Novel w/ gaps: {}", as_novel_shared.read().unwrap());
    println!("Annot w/ gaps: {}", as_annot_shared.read().unwrap());

    let num_AS = as_shared.read().unwrap().len();
    let mut avg_sq:Vec<i32> = Vec::new();
    let mut group_avg_sq:HashMap<String, i32> = HashMap::new();

    /*
    let avg_sq_file = output.clone() + ".group_averages.csv";
    let avg_sq_path = Path::new(&avg_sq_file);
    let mut avg_sq_outline = "average_sqs,num_annot,num_novel\n".to_string();


     */   

    let cov_file = output.clone() + ".coverage_all.csv";
    let cov_path = Path::new(&cov_file);
    let mut cov_outline = "score,coverage,type\n".to_string();

    let sim_ceiling = B * (1.0 - k as f32 *e);
    let max_diff = d * (1.0 - k as f32*e);

    println!("Max diff is {}", max_diff);

    for iso in ref_map.keys(){
        let mut scores:Vec<f32> = Vec::new();
        let num_reads = ref_map.get(iso).expect("msg").len();

        let mut num_incl = 0;

        for read in ref_map.get(iso).expect("msg").iter(){
            let diff = (read.coverage * (1.0 - k as f32 *e)) - read.similarity;
            scores.push(read.similarity);  

            
            if read.sample_id.contains("_novel"){
                cov_outline.push_str(&(read.similarity.to_string() + "," + &(read.coverage.to_string()) + ",novel\n"));
            }
            else{
                cov_outline.push_str(&(read.similarity.to_string() + "," + &(read.coverage.to_string()) + ",annot\n"));
            }
             
        }   
        let score_avg = (1000.0 * scores.iter().sum::<f32>() / scores.len() as f32) as i32;

        if score_avg >= 0{
            group_avg_sq.entry(iso.clone()).or_insert(score_avg*score_avg);
            for _i in 0..num_reads{
                avg_sq.push(score_avg*score_avg);
            }  
        }
    }

    let mut ATSS:HashMap<String, usize> = HashMap::new();
    let mut ATSS_novel = 0;
    let mut ATSS_annot = 0;

    avg_sq.sort();
    let cand_to_select = (c * num_AS as f64) as usize;

    println!("Cand to select: {}", cand_to_select);

    let mut exclude_d_novel = 0;
    let mut exclude_B_novel = 0;
    let mut exclude_d_annot = 0;
    let mut exclude_B_annot = 0;

    if cand_to_select > 0{
        let mut avg_sq_threshold = avg_sq.pop().expect("msg");

        if cand_to_select < avg_sq.len(){
            avg_sq_threshold = avg_sq[cand_to_select];
        }

        for (iso, avg_sq) in group_avg_sq.into_iter(){
            let mut group_annot = 0;
            let mut group_novel = 0;

            let matches = ref_map.get(&iso).expect("msg");

            /*
            for read in matches.iter(){
                if read.sample_id.contains("novel"){
                group_novel += 1;
                }
                else{
                    group_annot += 1;
                }
            }

            avg_sq_outline.push_str(&(avg_sq.to_string() + "," +  &group_annot.to_string() + "," + &group_novel.to_string() + "\n"));
             */
            
            if avg_sq <= avg_sq_threshold{
                for read in matches.iter(){
                    if read.similarity < sim_ceiling{

                        /*
                        if read.sample_id.contains("_novel"){
                            cov_outline_2.push_str(&(read.similarity.to_string() + "," + &(read.coverage.to_string()) + ",novel\n"));
                        }
                        else{
                            cov_outline_2.push_str(&(read.similarity.to_string() + "," + &(read.coverage.to_string()) + ",annot\n"));
                        }
                         */
                        

                        let diff = (read.coverage * (1.0 - k as f32 *e)) - read.similarity;

                        if diff >= 0.8{
                            println!("Ref is {}", iso);

                            large_diff.entry(read.sample_id.clone()).or_insert((iso.clone(), read.chunk));
                        }

                        if diff <= max_diff{
                            ATSS.entry(read.sample_id.clone()).or_insert(read.chunk);  

                            if read.sample_id.contains("novel"){
                                ATSS_novel += 1;
                            }
                            else{
                                ATSS_annot += 1;
                            }
                        }
                        else{
                            if read.sample_id.contains("novel"){
                                exclude_d_novel += 1;
                            }
                            else{
                                exclude_d_annot += 1;
                            }
                        }                        
                    }
                    else {

                        if read.sample_id.contains("novel"){
                                exclude_B_novel += 1;
                            }
                            else{
                                exclude_B_annot += 1;
                            }
                    }
                }
            }
        }        
    }
    
    //fs::write(avg_sq_path, avg_sq_outline);
    //fs::write(stats_path, stats_outline);

    fs::write(cov_path, cov_outline);
    //fs::write(cov_path_2, cov_outline_2);

let as_novel = as_novel_shared.read().unwrap().clone();
let as_annot = as_annot_shared.read().unwrap().clone(); 
(as_shared, ATSS, as_novel, as_annot, ATSS_novel, ATSS_annot, exclude_d_novel, exclude_d_annot, exclude_B_novel, exclude_B_annot, large_diff)
 }

pub fn kmer_chain(sample_kmers:&HashMap<i64, Vec<i32>>, ref_kmers:&HashMap<i64, Vec<i32>>, m:i32) -> (f32, i32, i32, i32, i32, i32){
    
    let mut g:Vec<(f32, i32, i32, i32, i32, i32)> = Vec::new();  
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
        let mut optimal_chain = 1.0;
        let mut optimal_gap = 0;
        let mut first_pos_q = anchors[i].0;
        let mut last_pos_q = anchors[i].0;
        let mut first_pos_r = anchors[i].1;
        let mut last_pos_r = anchors[i].1;


        if i == 0{            
            g.push((optimal_chain, optimal_gap, first_pos_q, first_pos_r, last_pos_q, last_pos_r));
        }
        else{
            let mut unique_idx = 0;
            let mut l = 0;

            while unique_idx < m{
                let j:i32 = i as i32 - l - 1;
                l += 1;

                if j >= 0{
                    let j = j as usize;

                    if anchors[j].0 != anchors[j+1].0 && anchors[j].1 != anchors[j+1].1{
                        unique_idx += 1;
                    }
                    if (anchors[j].1 < anchors[i].1) && (anchors[j].0 < anchors[i].0){
                        let gap = ((anchors[j].1 - anchors[i].1) - (anchors[j].0 - anchors[i].0)).abs();
                        let chain = g[j as usize].0 + 1.0;
                        
                        if chain > optimal_chain{
                            optimal_gap = std::cmp::max(gap, g[j as usize].1 as i32);
                            optimal_chain = chain;

                            first_pos_q = g[j as usize].2;
                            first_pos_r = g[j as usize].3;
                        }                       
                    }
                }
                else{
                    unique_idx = m;
                }
            }
            g.push((optimal_chain, optimal_gap, first_pos_q, first_pos_r, last_pos_q, last_pos_r));
        }
    }    
    g.pop().expect("msg")
}
