/*
use std::collections::HashMap;
use std::sync::Mutex;
use crate::types::*;
use std::sync::{Arc, RwLock};
use crate::types;
use std::path::Path;
use std::fs;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use rand::seq::{IndexedRandom, IteratorRandom};
//use average::Variance;
//use math::stats;

pub fn filter(sample_data_shared:Arc<RwLock<Vec<HashMap<String, String>>>>, ref_map_shared:Arc<RwLock<HashMap<String, Vec<Match>>>>, s:f64, num_exon_gaps:usize, v:f64)
-> (HashMap<String, usize>){
//->(HashMap<String, usize>, HashMap<String, usize>, HashMap<String, usize>, HashMap<String, usize>, HashMap<String, usize>, HashMap<String, usize>, HashMap<String, usize>){

    println!("Beginning to filter");

    let f = "reads_per_isoform.csv".to_string();
    let f_path = Path::new(&f);
    let mut line = "isoform,novel,annotated,variance\n".to_string();

    let s_file = "similarity_scores_ge_5_le_10.csv".to_string();
    let s_path = Path::new(&s_file);
    let mut s_line = "similarity,coverage\n".to_string();

    let s_40 = "similarity_scores_var_40.csv".to_string();
    let s_40_path = Path::new(&s_40);
    let mut s_40_line = "type,similarity,coverage\n".to_string();

    let diff = "diff_btw_coverage_score.csv".to_string();
    let diff_path = Path::new(&diff);
    let mut diff_line = "type,difference\n".to_string();


    println!("{}", f);
    println!("{}", line);

    let sample_data = sample_data_shared.read().unwrap();
    let ref_map = ref_map_shared.read().unwrap();

    let mut iso_ge_5 = 0;
    let mut read_ge_5 = 0;

    let mut novel_le_all = 0;
    let mut annot_le_all = 0;

    let mut novel_le = 0;
    let mut annot_le = 0;

    let mut novel_chosen = 0;
    let mut annot_chosen = 0;

    //let mut num_exon_gaps = 0;
    let mut total_reads = 0;
    for chunk in sample_data.iter(){
        total_reads += chunk.len();
    }
    //println!("Total reads is {}", total_reads);
    let mut sim_scores:Vec<i32> = Vec::new();
    let mut all_sim_scores:Vec<i32> = Vec::new();
    

    let mut sim_scores_10:Vec<i32> = Vec::new();
    let mut subsampled_reads_10:Vec<types::Match> = Vec::new();
    let mut sim_scores_15:Vec<i32> = Vec::new();
    let mut subsampled_reads_20:Vec<types::Match> = Vec::new();
    let mut sim_scores_20:Vec<i32> = Vec::new();
    let mut subsampled_reads_30:Vec<types::Match> = Vec::new();
    let mut sim_scores_25:Vec<i32> = Vec::new();
    let mut subsampled_reads_40:Vec<types::Match> = Vec::new();
    let mut sim_scores_30:Vec<i32> = Vec::new();
    let mut subsampled_reads_50:Vec<types::Match> = Vec::new();
    let mut sim_scores_35:Vec<i32> = Vec::new();
    let mut subsampled_reads_60:Vec<types::Match> = Vec::new();
    let mut subsampled_reads_70:Vec<types::Match> = Vec::new();

    let f_10 = "similarity_scores_random_diff_10.csv".to_string();
    let f_path_10 = Path::new(&f_10);
    let mut line_10 = "type,score,coverage\n".to_string();
    let f_15 = "similarity_scores_random_diff_15.csv".to_string();
    let f_path_15 = Path::new(&f_15);
    let mut line_15 = "type,score,coverage\n".to_string();
    let f_20 = "similarity_scores_random_diff_20.csv".to_string();
    let f_path_20 = Path::new(&f_20);
    let mut line_20 = "type,score,coverage\n".to_string();
    let f_25 = "similarity_scores_random_diff_25.csv".to_string();
    let f_path_25 = Path::new(&f_25);
    let mut line_25 = "type,score,coverage\n".to_string();
    let f_30 = "similarity_scores_random_diff_30.csv".to_string();
    let f_path_30 = Path::new(&f_30);
    let mut line_30 = "type,score,coverage\n".to_string();
    let f_60 = "similarity_scores_random_diff_60.csv".to_string();
    let f_60 = Path::new(&f_60);
    let mut line_35 = "type,score,coverage\n".to_string();

    //println!("Num. of Isoforms in Map: {}", ref_map.len());
    let mut i_count = 0;
    let mut subsampled_reads:Vec<types::Match> = Vec::new();


    let mut criteria_2:HashMap<String, usize> = HashMap::new();
    let mut cov_40:HashMap<String, usize> = HashMap::new();

    /*
    let mut novel_40 = 0;
    let mut annot_40 = 0;

    let mut novel_30 = 0;
    let mut annot_30 = 0;
    let mut novel_50 = 0;
    let mut annot_50 = 0;
    let mut novel_60 = 0;
    let mut annot_60 = 0;
    let mut novel_70 = 0;
    let mut annot_70 = 0;

    let seed: u64 = 42;
    let mut rng = StdRng::seed_from_u64(seed);

    let mut novel_cov = 0;
    let mut annot_cov = 0;

    for (ref_id, sample_matches) in ref_map.iter(){
        let mut isoform_reads:Vec<types::Match> = Vec::new();
        let mut novel = 0;
        let mut annot = 0;
        let mut scores:Vec<f64> = Vec::new();
        let mut scores_sq:Vec<f64> = Vec::new();

        for read in sample_matches.iter(){
            scores.push(read.similarity as f64 / 10000.0);
            scores_sq.push((read.similarity as f64 / 10000.0) * (read.similarity as f64 / 10000.0));

            let d = read.coverage - (read.similarity as f32 / 10000.0);
            /*
            if read.sample_id.contains("_novel"){
                novel += 1;
                diff_line.push_str(&("novel,".to_string() + &(d.to_string()) + "\n"));
            }
            else{
                annot += 1;
                diff_line.push_str(&("annotated,".to_string() + &(d.to_string()) + "\n"));
            }
             */
            
        }

        let score_avg:f64 = scores.iter().sum::<f64>() / scores.len() as f64;
        let var:f64 = (scores_sq.iter().sum::<f64>() / scores_sq.iter().sum::<f64>()) - (score_avg * score_avg);
        line.push_str(&(ref_id.to_string() + "," + &novel.to_string() + "," + &annot.to_string() + "," + &var.to_string() + "\n"));

        if var >= 0.3{
            for read in sample_matches.iter(){
                if read.sample_id.contains("_novel"){
                    novel_30 += 1;
                }
                else{
                    annot_30 += 1;
                }
            }
        }
        if var >= 0.5{
            for read in sample_matches.iter(){
                if read.sample_id.contains("_novel"){
                    novel_50 += 1;
                }
                else{
                    annot_50 += 1;
                }
            }
        }
        if var >= 0.6{
            for read in sample_matches.iter(){
                if read.sample_id.contains("_novel"){
                    novel_60 += 1;
                }
                else{
                    annot_60 += 1;
                }
            }
        }
        if var >= 0.7{
            for read in sample_matches.iter(){
                if read.sample_id.contains("_novel"){
                    novel_70 += 1;
                }
                else{
                    annot_70 += 1;
                }
            }
        }

        if var >= v{
            for read in sample_matches.iter(){
                if read.sample_id.contains("_novel"){
                    novel_40 += 1;
                    //s_40_line.push_str(&("novel,".to_string() + &(read.similarity.to_string()) + "," + &(read.coverage.to_string()) + "\n"));
                }
                else{
                    annot_40 += 1;
                    //s_40_line.push_str(&("annot,".to_string() + &(read.similarity.to_string()) + "," + &(read.coverage.to_string()) + "\n"));
                }

                let d = read.coverage - (read.similarity as f32 / 10000.0);
                if d <= 0.15{
                    isoform_reads.push(read.clone());
                    sim_scores.push(read.similarity);
                    subsampled_reads_40.push(read.clone());

                    if read.sample_id.contains("_novel"){
                        novel_cov += 1;
                    }
                    else{
                        annot_cov += 1;
                    }       
                }  
            }
        }
    }

    let cand_to_select = ((s * total_reads as f64) - num_exon_gaps as f64) as usize;

    //println!("Isoforms with at least 5 read queries: {}", iso_ge_5);
    //println!("Query reads: {}", read_ge_5);

    sim_scores.sort();
    //let cand_to_select = (s * total_reads as f64) - num_exon_gaps as f64;

    let mut threshold_40 = 0;
    if sim_scores.len() < cand_to_select{
        threshold_40 = sim_scores.pop().expect("msg");
    }
    else{
        threshold_40 = sim_scores[cand_to_select as usize];
    }

    for read in subsampled_reads_40.iter(){
        //let r: u32 = rng.random_range(1..subsampled_reads_40.len() as u32);
        
        if read.similarity <= threshold_40{
            cov_40.entry(read.sample_id.clone()).or_insert(read.chunk);  

            if read.sample_id.contains("_novel"){
                novel_chosen += 1;
            }
            else{
                annot_chosen += 1;
            }
        }
    }
    

    /*
    sim_scores.sort();
    sim_scores_10.sort();
    sim_scores_15.sort();
    sim_scores_20.sort();
    sim_scores_25.sort();
    sim_scores_30.sort();
    sim_scores_35.sort();
    all_sim_scores.sort();
    //let mut threshold_collapse = -1.0;
    //let threshold_no_collapse = all_sim_scores[cand_to_select as usize] as f32 / 10000.0;

    let mut threshold_collapse = -1;
    let mut threshold_collapse_10 = -1;
    let mut threshold_collapse_15 = -1;
    let mut threshold_collapse_20 = -1;
    let mut threshold_collapse_25 = -1;
    let mut threshold_collapse_30 = -1;
    let mut threshold_collapse_35 = -1;
   // println!("{}", sim_scores.len());
    //println!("{}", all_sim_scores.len());

    println!("subsampled scores, : {}", sim_scores.len());
    println!("subsampled scores, max diff = .1: {}", sim_scores_10.len());
    println!("subsampled scores, max diff = .15: {}", sim_scores_15.len());
    println!("subsampled scores, max diff = .2:  {}", sim_scores_20.len());
    println!("subsampled scores, max diff = .25: {}", sim_scores_25.len());
    println!("subsampled scores, max diff = .3:  {}", sim_scores_30.len());
    println!("subsampled scores, max diff = .35: {}", sim_scores_35.len());

    println!("{}", threshold_collapse);
    println!("{}", threshold_collapse_10);
    println!("{}", threshold_collapse_15);
    println!("{}", threshold_collapse_20);
    println!("{}", threshold_collapse_25);
    println!("{}", threshold_collapse_30);
    println!("{}", threshold_collapse_35);
    */
    //println!("{}", num_exon_gaps);              //1972
    //println!("Cand to select: {}", cand_to_select);             //46292.5
    //println!("Collapse threshold is {}", threshold_collapse  as f32 / 10000.0);
    //println!("No Collapse threshold is {}", threshold_no_collapse  as f32 / 10000.0);

    let mut crit_2_novel = 0;
    let mut crit_2_annot = 0;
    let mut crit_2_novel_10 = 0;
    let mut crit_2_annot_10 = 0;
    let mut crit_2_novel_15 = 0;
    let mut crit_2_annot_15 = 0;
    let mut crit_2_novel_20 = 0;
    let mut crit_2_annot_20 = 0;
    let mut crit_2_novel_25 = 0;
    let mut crit_2_annot_25 = 0;
    let mut crit_2_novel_30 = 0;
    let mut crit_2_annot_30 = 0;
    let mut crit_2_novel_35 = 0;
    let mut crit_2_annot_35 = 0;

    /*
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
    for read in subsampled_reads_10{
        if read.similarity <= threshold_collapse_10{
            if read.sample_id.contains("novel"){
                crit_2_novel_10 += 1;
            }
            else{
                crit_2_annot_10 += 1;
            }
            cov_10.entry(read.sample_id).or_insert(read.chunk);  
        }
    }
    for read in subsampled_reads_20{
        if read.similarity <= threshold_collapse_15{
            if read.sample_id.contains("novel"){
                crit_2_novel_15 += 1;
            }
            else{
                crit_2_annot_15 += 1;
            }
            cov_20.entry(read.sample_id).or_insert(read.chunk);  
        }
    }
    for read in subsampled_reads_30{
        if read.similarity <= threshold_collapse_20{
            if read.sample_id.contains("novel"){
                crit_2_novel_20 += 1;
            }
            else{
                crit_2_annot_20 += 1;
            }
            cov_30.entry(read.sample_id).or_insert(read.chunk);  
        }
    }
    for read in subsampled_reads_40{
        if read.similarity <= threshold_collapse_25{
            if read.sample_id.contains("novel"){
                crit_2_novel_25 += 1;
            }
            else{
                crit_2_annot_25 += 1;
            }
            cov_40.entry(read.sample_id).or_insert(read.chunk);  
        }
    }
    for read in subsampled_reads_50{
        if read.similarity <= threshold_collapse_30{
            if read.sample_id.contains("novel"){
                crit_2_novel_30 += 1;
            }
            else{
                crit_2_annot_30 += 1;
            }
            cov_50.entry(read.sample_id).or_insert(read.chunk);  
        }
    }
    for read in subsampled_reads_60{
        if read.similarity <= threshold_collapse_35{
            if read.sample_id.contains("novel"){
                crit_2_novel_35 += 1;
            }
            else{
                crit_2_annot_35 += 1;
            }
            cov_60.entry(read.sample_id).or_insert(read.chunk);  
        }
    }
    */
    /*
        println!("Novel w/ 3'/5' AS: {}", crit_2_novel);
    println!("Annot w/ 3'/5' AS: {}", crit_2_annot);
    println!("Max diff = .1");
    println!("Novel w/ 3'/5' AS: {}", crit_2_novel_10);
    println!("Annot w/ 3'/5' AS: {}", crit_2_annot_10);
    println!("Max diff = .15");
    println!("Novel w/ 3'/5' AS: {}", crit_2_novel_15);
    println!("Annot w/ 3'/5' AS: {}", crit_2_annot_15);
    println!("Max diff = .2");
    println!("Novel w/ 3'/5' AS: {}", crit_2_novel_20);
    println!("Annot w/ 3'/5' AS: {}", crit_2_annot_20);
    println!("Max diff = .25");
    println!("Novel w/ 3'/5' AS: {}", crit_2_novel_25);
    println!("Annot w/ 3'/5' AS: {}", crit_2_annot_25);
    println!("Max diff = .3");
    println!("Novel w/ 3'/5' AS: {}", crit_2_novel_30);
    println!("Annot w/ 3'/5' AS: {}", crit_2_annot_30);
    println!("Max diff = .35");
    println!("Novel w/ 3'/5' AS: {}", crit_2_novel_35);
    println!("Annot w/ 3'/5' AS: {}", crit_2_annot_35);
     */

    /*
    println!("Novel reads ge 5: {}", novel_le_all);
    println!("Annotated reads ge 5: {}", annot_le_all);

    println!("Novel reads le: {}", novel_le);
    println!("Annotated reads le:  {}", annot_le);

    println!("Novel reads chosen: {}", novel_chosen);
    println!("Annotated reads chosen: {}", annot_chosen);
     */

    println!("Novel reads within .15 diff: {}", novel_cov);
    println!("Annot reads within .15 diff: {}", annot_cov);

    println!("Novel reads w/ var > .3: {}", novel_30);
    println!("Annot reads w/ var > .3: {}", annot_30);
    println!("Novel reads w/ var > {} : {}", v, novel_40);
    println!("Annot reads w/ var > {}: {}", v, annot_40);
    println!("Novel reads w/ var > .5: {}", novel_50);
    println!("Annot reads w/ var > .5: {}", annot_50);
    println!("Novel reads w/ var > .6: {}", novel_60);
    println!("Annot reads w/ var > .6: {}", annot_60);
    println!("Novel reads w/ var > .7: {}", novel_70);
    println!("Annot reads w/ var > .7: {}", annot_70);

    println!("Novel reads chosen: {}", novel_chosen);       //from var > .4
    println!("Annot reads chosen: {}", annot_chosen);
    

    //println!("{}", f);
    //println!("{}", line);

    //fs::write(f_path, line);
    //fs::write(s_40_path, s_40_line);
    //fs::write(diff_path, diff_line);
    //fs::write(s_path, s_line);
    //fs::write(f_path_10, line_10);
    //fs::write(f_path_15, line_15);
    //fs::write(f_path_20, line_20);
    //fs::write(f_path_25, line_25);
    //fs::write(f_path_30, line_30);
    //fs::write(f_60, line_35);
    */
    (cov_40)
    
}
*/