use std::collections::HashMap;
use std::sync::Mutex;
use crate::types::*;
use std::sync::{Arc, RwLock};
use crate::types;
use std::path::Path;
use std::fs;

pub fn filter(sample_data_shared:Arc<RwLock<Vec<HashMap<String, String>>>>, ref_map_shared:Arc<RwLock<HashMap<String, Vec<Match>>>>, nb:i32, bin_s:f32, s:f64, num_exon_gaps:usize)
-> HashMap<String, usize>{

    let sample_data = sample_data_shared.read().unwrap();
    let ref_map = ref_map_shared.read().unwrap();

    let mut total_reads = 0;
    for chunk in sample_data.iter(){
        total_reads += chunk.len();
    }
    //println!("Total reads is {}", total_reads);
    let mut sim_scores:Vec<i32> = Vec::new();
    let mut all_sim_scores:Vec<i32> = Vec::new();
    let mut subsampled_reads:Vec<types::Match> = Vec::new();

    //println!("Num. of Isoforms in Map: {}", ref_map.len());
    let mut i_count = 0;

    let mut criteria_2:HashMap<String, usize> = HashMap::new();

    for (ref_id, sample_matches) in ref_map.iter(){
        
        let mut bins:Vec<types::Bin> = Vec::new();
        for i in 1 .. nb+1{     
            bins.push(types::build_Bin(((i as f32) / nb as f32), Vec::new(), false, 2.0, 0, 0, Vec::new()))
        }

        for read in sample_matches.iter(){
            for i in 1 .. (nb+1) as usize{
                if (read.similarity as f32 / 10000.0) >= ((i as f32 - 1.0) / nb as f32) && (read.similarity as f32 / 10000.0) < ((i as f32) / nb as f32){
                            
                    bins[i-1].num_total += 1;
                    bins[i-1].scores.push(read.similarity);
                    bins[i-1].matches.push(read.clone());
                }   
            }
        }
            
        for bin in bins{
            if bin.num_total > 0{
                let mut bin_scores = bin.scores;
                bin_scores.sort();
                let bin_sample = ((bin_scores.len() as f32 * bin_s).ceil() as usize);
                let mut selected = 0;
                let mut subsample = -1;
                if bin_s > 0.0{
                    subsample = bin_scores[bin_sample - 1];
                }

                for read in bin.matches{
                    all_sim_scores.push(read.similarity);
                    if read.similarity <= subsample{   
                        sim_scores.push(read.similarity);
                        subsampled_reads.push(read.clone());
                        selected += 1;
                    }
                }
            }
        }
    }

    let cand_to_select = (s * total_reads as f64) - num_exon_gaps as f64;

    sim_scores.sort();
    all_sim_scores.sort();

    let mut threshold_collapse = -1;
    let threshold_no_collapse = all_sim_scores[cand_to_select as usize];

    if bin_s > 0.0{
        if cand_to_select >= sim_scores.len() as f64{
            threshold_collapse = *sim_scores.last().expect("msg");
        }
        else{
            threshold_collapse = sim_scores[cand_to_select as usize];
        }
    }
    
    //println!("Collapse threshold is {}", threshold_collapse  as f32 / 10000.0);
    //println!("No Collapse threshold is {}", threshold_no_collapse  as f32 / 10000.0);

    let mut crit_2_novel = 0;
    let mut crit_2_annot = 0;

    for read in subsampled_reads{
        if read.similarity <= threshold_collapse{
            if read.sample_id.contains("novel"){
                crit_2_novel += 1;
            }
            else{
                crit_2_annot += 1;
            }
            criteria_2.entry(read.sample_id).or_insert(read.chunk);  
        }
    }
    //println!("Novel w/ 3'/5' AS: {}", crit_2_novel);
    //println!("Annot w/ 3'/5' AS: {}", crit_2_annot);

    criteria_2
}
