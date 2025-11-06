use clap::builder::StringValueParser;
use rand::seq::IndexedRandom;
use rand::seq::IteratorRandom;

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
use rand;

pub fn create_kmer_map(ref_sketch:&HashMap<String, Sketch>, b:usize ) -> (HashMap<i64, HashSet<String>>){

    let mut kmer_map:HashMap<i64, HashSet<String>> = HashMap::new();
    let mut to_remove:HashSet<i64> = HashSet::new();

    for (id, sketch) in ref_sketch.iter(){
        for kmer in sketch.kmers.keys(){
            kmer_map.entry(*kmer).or_insert(HashSet::new()).insert(id.clone());
            if kmer_map.entry(*kmer).or_insert(HashSet::new()).len() >= b{
                to_remove.insert(*kmer);
            }
        }
    }
    
    for kmer in to_remove.iter(){
        kmer_map.remove_entry(kmer);
    }

    (kmer_map)
}

pub fn find_ref_matches(output:String, ref_sketches_shared: Arc<RwLock<HashMap<String, Sketch>>>, sample_chunks_shared:Arc<RwLock<Vec<HashMap<String, String>>>>, b:usize, t: i32, n:i32, k:i32, a:i64, s:f64, m:i32, total_reads:i32, e:f32, B:f32, c:f64, mode:String, l:usize)
 -> (Arc<RwLock<HashMap<String, usize>>>, HashMap<String, usize>) {

    let ref_map:HashMap<String, Vec<Match>> = HashMap::new(); 
    let AS:HashMap<String, usize> = HashMap::new();

    let ref_sketches = ref_sketches_shared.read().unwrap();
    let (kmer_map) = create_kmer_map(&ref_sketches, b);
    drop(ref_sketches);
    let mut handles = vec![];   
    let as_novel = 0;
    let as_annot = 0;

    let exon_novel = 0;
    let no_match_novel = 0;
    let exon_annot = 0;
    let r_gap_annot = 0;
    let q_gap_annot = 0;
    let no_match_annot = 0;
    
    let kmer_map_shared = Arc::new(RwLock::new(kmer_map));
    let ref_map_shared = Arc::new(RwLock::new(ref_map));
    let as_shared = Arc::new(RwLock::new(AS));
    
    let as_novel_shared = Arc::new(RwLock::new(as_novel));
    let as_annot_shared = Arc::new(RwLock::new(as_annot));
    let exon_novel_shared = Arc::new(RwLock::new(exon_novel));
    let exon_annot_shared = Arc::new(RwLock::new(exon_annot));
    let r_gap_annot_shared = Arc::new(RwLock::new(r_gap_annot));
    let q_gap_annot_shared = Arc::new(RwLock::new(q_gap_annot));
    let no_match_annot_shared = Arc::new(RwLock::new(no_match_annot));
    let no_match_novel_shared = Arc::new(RwLock::new(no_match_novel));

    let end_gap_file = output.clone() + ".end_gap_reads.csv";
    let end_gap_path = Path::new(&end_gap_file);
    let mut end_gap_outline = "id,opt ref gap,query gap from 5p,query gap from 3p,internal gap,optimal score\n".to_string();
    let end_gap_outline_shared = Arc::new(RwLock::new(end_gap_outline));

    let ceiling_file = output.clone() + ".score_ceiling.csv";
    let ceiling_path = Path::new(&ceiling_file);
    let mut ceiling_outline = "read,score\n".to_string();
    let ceiling_outline_shared = Arc::new(RwLock::new(ceiling_outline));

    for i in 0..t as usize{
        let sample_chunks_shared = Arc::clone(&sample_chunks_shared); 
        let kmer_map_shared = Arc::clone(&kmer_map_shared);
        let ref_sketches_shared = Arc::clone(&ref_sketches_shared);
        let ref_map_shared = Arc::clone(&ref_map_shared);
        let as_shared = Arc::clone(&as_shared);

        let as_novel_shared = Arc::clone(&as_novel_shared);
        let as_annot_shared = Arc::clone(&as_annot_shared);

        let exon_novel_shared = Arc::clone(&exon_novel_shared);
        let exon_annot_shared = Arc::clone(&exon_annot_shared);
        let r_gap_annot_shared = Arc::clone(&r_gap_annot_shared);
        let q_gap_annot_shared = Arc::clone(&q_gap_annot_shared);
        let no_match_annot_shared = Arc::clone(&no_match_annot_shared);
        let no_match_novel_shared = Arc::clone(&no_match_novel_shared);
        let end_gap_outline_shared = Arc::clone(&end_gap_outline_shared);
        let ceiling_outline_shared = Arc::clone(&ceiling_outline_shared);
        

        let handle = thread::spawn(move || {
            let guard = sample_chunks_shared.read().unwrap();
            let sample_chunk = (guard).get(i).expect("msg");
            let kmer_map = kmer_map_shared.read().unwrap();      
            let ref_sketches = ref_sketches_shared.read().unwrap();
            
            let r_count = 0;

            for (id, seq) in sample_chunk.iter(){

                let mut tail_len = 0;
                let mut chars = seq.chars().rev();
                let mut c = chars.next().expect("msg");
                while c=='A' || c=='a'{
                    tail_len += 1;
                    c = chars.next().expect("msg");
                }
                let adjusted_len = seq.len() as i32 - tail_len;

                let (size, sketch) = frac_min_hash(&seq, k, a,s);
                let mut no_match = 0;

                let mut opt_similarity = -1.0 * f32::INFINITY;   
                let mut opt_num_matches = 0;         
                let mut opt_ref_match = "".to_string();
                let mut opt_ref_len = 0;
                let mut opt_gap = 0;
                let mut opt_query_gap = 0;
                let mut opt_chain = 0;
                let mut opt_ref_gap = 0;
                
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
                        let (chain_score, gap, (first_pos_q, first_pos_r), (last_pos_q, last_pos_r)) = kmer_chain(&(sketch),  &r_sketch.kmers, m);

                        let gap_5p = std::cmp::max(0,first_pos_q - first_pos_r);
                        let gap_3p = std::cmp::max(0,(adjusted_len - last_pos_q) - (r_sketch.seq_len as i32 - last_pos_r));

                        let mut ref_gap = 0;
                        if std::cmp::max(first_pos_r - first_pos_q, ((r_sketch.seq_len as i32 - last_pos_r) - (adjusted_len - last_pos_q))) > 0{
                            ref_gap = std::cmp::max(first_pos_q, (adjusted_len - last_pos_q - k + 1));
                        }
                        
                        let sim_score = chain_score as f32 / ((size + r_sketch.size) as f32 - chain_score as f32);
                        if sim_score > opt_similarity{
                            opt_similarity = sim_score;
                            opt_ref_match = ref_id.clone();
                            opt_ref_len = r_sketch.seq_len;
                            opt_gap = gap;
                            opt_chain = chain_score;
                            opt_query_gap = std::cmp::max(gap_3p, gap_5p);
                            opt_ref_gap = ref_gap;
                            opt_num_matches = chain_score as i32;

                        }
                    }              
                }

                
                let mut end_gap_outline = end_gap_outline_shared.write().unwrap();
                end_gap_outline.push_str(&(id.to_string() + "," + &(opt_ref_gap).to_string() + "," + &(opt_query_gap).to_string() + "," + &(opt_gap).to_string() + "," + &(opt_similarity).to_string() + "\n"));
                drop(end_gap_outline);
                

                let mut ceiling_outline = ceiling_outline_shared.write().unwrap();
                ceiling_outline.push_str(&(id.to_string() + "," + &(opt_similarity.to_string()) + "\n"));

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
                
                else if opt_num_matches == 0{
                    if adjusted_len >= 250{
                        let mut AS = as_shared.write().unwrap();
                        if id.contains("novel"){
                            *no_match_novel_shared.write().unwrap() += 1;
                        }
                        else{
                            *no_match_annot_shared.write().unwrap() += 1;
                        }
                    }
                }
                
                else if  opt_query_gap > n{
                    let mut AS = as_shared.write().unwrap();
                    AS.entry(id.clone()).or_insert(i);
                    if id.contains("novel"){
                        *exon_novel_shared.write().unwrap() += 1;   
                    }
                    else{
                        *q_gap_annot_shared.write().unwrap() += 1;
                    }             
                }

                else if opt_ref_gap > 100{
                    let mut AS = as_shared.write().unwrap();
                    AS.entry(id.clone()).or_insert(i);
                    if id.contains("novel"){
                        *exon_novel_shared.write().unwrap() += 1;
                    }
                    else{
                        *r_gap_annot_shared.write().unwrap() += 1;
                    }  
                }
                                
                else{
                    let mut ref_map = ref_map_shared.write().unwrap(); 
                    ref_map.entry(opt_ref_match.clone()).or_insert(Vec::new()).push(build_match(id.clone(), opt_similarity, i, opt_chain as f32, adjusted_len as usize));
                    drop(ref_map);
                }      
            }
        });
        handles.push(handle);
    }

    for i in handles{
        i.join().unwrap();
    } 

    let mut ref_map = ref_map_shared.write().unwrap();
    let mut additional:HashMap<String, usize> = HashMap::new();

    println!("Annot AS: {}", as_annot_shared.read().unwrap());
    println!("Annot ref gap: {}", r_gap_annot_shared.read().unwrap());
    println!("Annot query gap: {}", q_gap_annot_shared.read().unwrap());
    println!("Novel no match: {}", no_match_novel_shared.read().unwrap());
    println!("Annot no match: {}", no_match_annot_shared.read().unwrap());

    if mode == "quantify".to_string(){

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
                    additional.entry(read.sample_id.clone()).or_insert(read.chunk);
                }               
                scale_outline.push_str(&(iso.clone()+ "," + &matches.len().to_string() + "," + &scale.to_string() + "\n"));
            }
            else{
                for read in matches.iter(){
                    additional.entry(read.sample_id.clone()).or_insert(read.chunk);
                }
                let scale = matches.len() as f64 / l as f64;
                scale_outline.push_str(&(iso.clone()+ "," + &matches.len().to_string() + "," + &scale.to_string() + "\n"));
            } 
        }
        fs::write(scale_factor, scale_outline);  
    }

    if mode == "discovery".to_string(){

        println!("Novel w/ gaps: {}", as_novel_shared.read().unwrap());
        println!("Annot w/ gaps: {}", as_annot_shared.read().unwrap());
        let num_AS = *as_novel_shared.read().unwrap() + *as_annot_shared.read().unwrap();

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

        let mut ATSS_novel = *exon_novel_shared.read().unwrap();
        let mut ATSS_annot = *exon_annot_shared.read().unwrap();

        avg_sq.sort();
        let mut cand_to_select = (c * num_AS as f64) as usize;
        if (ATSS_annot+ATSS_novel) > cand_to_select{
            cand_to_select = 0;
        }
        else{
            cand_to_select -= (ATSS_annot+ATSS_novel) as usize;
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
                                additional.entry(read.sample_id.clone()).or_insert(read.chunk);           
                                if read.sample_id.contains("novel"){
                                    ATSS_novel += 1;
                                }
                                else{
                                    ATSS_annot += 1;
                                }                                                                                  
                            }
                        }
                    }
                }        
            }
        }
        println!("ATSS novel: {}", ATSS_novel);
        println!("ATSS annot: {}", ATSS_annot);
    }

(as_shared, additional)
 }

pub fn kmer_chain(sample_kmers:&HashMap<i64, Vec<i32>>, ref_kmers:&HashMap<i64, Vec<i32>>, m:i32) -> (i32, i32, (i32, i32), (i32, i32)){
    
    let mut g_add:Vec<(i32, i32, i32, i32)> = Vec::new();  
    let mut g_not:Vec<(i32, i32, i32, i32)> = Vec::new();  
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
      
    for i in 0.. anchors.len() as i32{
        let mut opt_chain_add = 1;
        let mut opt_gap_add = 0;
        let mut first_add = i;
        let mut last_add = i;

        let mut opt_chain_not = 0;
        let mut opt_gap_not = 0;
        let mut first_not:i32 = -1;
        let mut last_not:i32 = -1;

        if i == 0{            
            g_add.push((opt_chain_add, opt_gap_add, first_add, last_add));
            g_not.push((opt_chain_not, opt_gap_not, first_not, last_not));
        }
        else{
            
            let mut unique_idx = 0;
            let mut l = 0;

            while unique_idx < m{
                let j = i - l - 1;
                l += 1;

                if j >= 0{
                    let j = j as usize;

                    if anchors[j].0 != anchors[j+1].0 && anchors[j].1 != anchors[j+1].1{
                        unique_idx += 1;
                    }

                    //if anchor i is not added to the optimal chain
                    if g_not[j].0 >= opt_chain_not{
                        opt_chain_not = g_not[j].0;
                        opt_gap_not = g_not[j].1;
                        first_not = g_not[j].2;
                        last_not = g_not[j].3;                        
                    }  
                    if g_add[j].0 >= opt_chain_not{
                        opt_chain_not = g_add[j].0;
                        opt_gap_not = g_add[j].1;
                        first_not = g_add[j].2;
                        last_not = g_add[j].3;                        
                    }                

                    //if anchor i is added to the optimal chain
                    let prev = g_add[j].3;
                    if prev != -1{
                        let prev = prev as usize;
                        if (anchors[prev].1 < anchors[i as usize].1) && (anchors[prev].0 < anchors[i as usize].0){
                            let gap = ((anchors[prev].1 - anchors[i as usize].1) - (anchors[prev].0 - anchors[i as usize].0)).abs();
                            let w_anchor = g_add[j].0 + 1;                        
                            
                            if w_anchor >= opt_chain_add{
                                opt_chain_add = w_anchor;
                                opt_gap_add = std::cmp::max(gap, g_add[j].1 as i32);                     
                                first_add = g_add[j].2;
                                last_add = i;                         
                            }                       
                        }
                    }
                    
                    let prev = g_not[j].3;
                    if prev != -1{
                        let prev = prev as usize;
                        if (anchors[prev].1 < anchors[i as usize].1) && (anchors[prev].0 < anchors[i as usize].0){
                            let gap = ((anchors[prev].1 - anchors[i as usize].1) - (anchors[prev].0 - anchors[i as usize].0)).abs();
                            let w_anchor = g_not[j].0 + 1;                        
                            
                            if w_anchor >= opt_chain_add{
                                opt_chain_add = w_anchor;
                                opt_gap_add = std::cmp::max(gap, g_not[j].1 as i32);                     
                                first_add = g_not[j].2;
                                last_add = i;                         
                            }                       
                        }
                    }
                }
                else{
                    unique_idx = m;
                }
            }
            
            g_add.push((opt_chain_add, opt_gap_add, first_add, last_add));
            g_not.push((opt_chain_not, opt_gap_not, first_not, last_not));
        }
    }    
    let ret_add = g_add.pop().expect("msg");
    let ret_not = g_not.pop().expect("msg");

    if ret_add.0 > ret_not.0{
        let first_anchor = anchors[ret_add.2 as usize];
        let last_anchor = anchors[ret_add.3 as usize];
        (ret_add.0, ret_add.1, first_anchor, last_anchor)
    }
    else{
        let first_anchor = anchors[ret_not.2 as usize];
        let last_anchor = anchors[ret_not.3 as usize];
        (ret_not.0, ret_not.1, first_anchor, last_anchor)
    }

}
