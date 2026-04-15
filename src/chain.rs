
use crate::types::*;
use std::collections::HashMap;

//Finds the largest set of colinear matches between a query and reference sketch
pub fn kmer_chain_unique(q_sketch:&HashMap<i64, Vec<i32>>, r_sketch:&HashMap<i64, Vec<i32>>, m:i32) -> (i32, i32, (i32, i32), (i32, i32), usize){    
    let mut g:Vec<DpCell> = Vec::new(); 
    let mut anchors:Vec<(i32, i32)> = Vec::new();
    let mut num_repeat = 0;
    let default:Vec<i32> = Vec::new();

    for (kmer, pos) in q_sketch{
        let ref_pos = r_sketch.get(kmer).unwrap_or_else(|| &default);
        if (pos.len() == 1) && (ref_pos.len() == 1){
            anchors.push((pos[0], ref_pos[0]));
        }
        else if ref_pos.len() > 0{
            num_repeat += ref_pos.len() + 1;
        }
    }     

    if anchors.len() == 0{
        num_repeat = 0;
        for (kmer, pos) in q_sketch{
            let ref_pos = r_sketch.get(kmer).unwrap_or_else(|| &default);
            if ref_pos.len() > 0{
                anchors.push((pos[0], ref_pos[0]));
            }
        }
    }

    anchors.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));   
    for i in 0.. anchors.len() as i32{
        let mut cell:DpCell = build_cell(1, 0, i, i);
        if i == 0{            
            g.push(cell);
        }
        else{
            for l in 0..m{
                let j = i - l - 1;
                if j >= 0{
                    let j = j as usize;             
                    //Determines optimal chain where anchor i is included
                    let prev = g[j].last_anchor;
                    if prev != -1{
                        let prev = prev as usize;
                        if (anchors[prev].1 < anchors[i as usize].1) && (anchors[prev].0 < anchors[i as usize].0){
                            let gap = ((anchors[prev].1 - anchors[i as usize].1) - (anchors[prev].0 - anchors[i as usize].0)).abs();
                            if g[j].OptChain_score + 1 >= cell.OptChain_score{
                                cell = build_cell(g[j].OptChain_score + 1, std::cmp::max(gap, g[j].opt_gap as i32), g[j].first_anchor, i);                        
                            }                   
                        }
                    }
                }
            }
            g.push(cell);
        }
    }    
    let ret = g.pop().expect("msg");    
    (ret.OptChain_score, ret.opt_gap, anchors[ret.first_anchor as usize], anchors[ret.last_anchor as usize], num_repeat)
}
