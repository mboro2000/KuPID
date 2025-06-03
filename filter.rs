use std::collections::HashMap;

use crate::types;
use std::path::Path;
use std::fs;

pub fn filter(sample_data:&Vec<HashMap<String, String>>, ref_map:HashMap<String, Vec<types::Match>>, output:String, nb:i32, m:i32, bin_s:f32, s:f64){

    let candidates_path = Path::new(&output);
    let mut candidates_outline = "".to_string();

    let mut num_exon_gaps = 0;
    let total_reads = sample_data.len();
    let mut sim_scores:Vec<i32> = Vec::new();
    let mut all_sim_scores:Vec<i32> = Vec::new();
    let mut subsampled_reads:Vec<types::Match> = Vec::new();

    println!("Num. of Isoforms in Map: {}", ref_map.len());
    let mut i_count = 0;

    for (ref_id, sample_matches) in ref_map.iter(){
        i_count += 1;
        println!("{}", i_count);
        let mut bins:Vec<types::Bin> = Vec::new();
        for i in 1 .. nb+1{     
            bins.push(types::build_Bin(((i as f32) / nb as f32), Vec::new(), false, 2.0, 0, 0, Vec::new()))
        }

        for samp_match in sample_matches.iter(){
            if samp_match.gap > m{
                num_exon_gaps += 1;
                candidates_outline.push_str(&(samp_match.sample_id.clone() + "\n" + &sample_data[samp_match.chunk].get(&samp_match.sample_id).expect("msg") + "\n"));
            }
            else{
                for read in sample_matches{
                    for i in 1 .. (nb+1) as usize{
                        if read.similarity >= ((i as f32 - 1.0) / nb as f32) && read.similarity < ((i as f32) / nb as f32){
                            /*
                            if read.read_type == "novel"{
                                bins[i-1].contain_novel = true;
                                bins[i-1].num_novel += 1;
                            }
                             */
                            
                            bins[i-1].num_total += 1;
                            if read.similarity < bins[i-1].min_val{
                                bins[i-1].min_val = read.similarity;
                            }       
                            bins[i-1].scores.push((10000.0 * read.similarity) as i32);
                            bins[i-1].matches.push(read.clone());
                        break;  //
                        }   
                    }
                }
            }
        }
        for bin in bins{
            if bin.num_total > 0{
                let mut bin_scores = bin.scores;
                bin_scores.sort();
                let subsample = bin_scores[((bin.num_total as f32 * bin_s).ceil() as usize - 1)] as f32 / 10000.0;

                for read in bin.matches{
                    all_sim_scores.push((10000.0 * read.similarity) as i32);
                    if read.similarity <= subsample{
                        sim_scores.push((10000.0 * read.similarity) as i32);
                        subsampled_reads.push(read.clone());
                    }
                }
            }
        }
    }

    let cand_to_select = (s * total_reads as f64) - num_exon_gaps as f64;

    sim_scores.sort();
    all_sim_scores.sort();
    let mut threshold_collapse = 0.0;
    let threshold_no_collapse = all_sim_scores[cand_to_select as usize] as f32 / 10000.0;

    if cand_to_select >= sim_scores.len() as f64{
        threshold_collapse = *sim_scores.last().expect("msg") as f32 / 10000.0;
    }
    else{
        threshold_collapse = sim_scores[cand_to_select as usize] as f32 / 10000.0;
    }
    
    println!("{}", num_exon_gaps);              //1972
    println!("Cand to select: {}", cand_to_select);             //46292.5
    println!("Collapse threshold is {}", threshold_collapse);
    println!("No Collapse threshold is {}", threshold_no_collapse);


    for read in subsampled_reads{
        if read.similarity <= threshold_collapse{
            candidates_outline.push_str(&(read.sample_id.clone() + "\n" + &sample_data[read.chunk].get(&read.sample_id).expect("msg") + "\n"));
        }
    }

    fs::write(candidates_path, candidates_outline);

}
