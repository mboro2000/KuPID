use crate::types::*;
use std::collections::HashMap;
use std::cmp;

//Finds the largest set of colinear matches between a query and reference sketch
pub fn kmer_chain(q_sketch:&HashMap<i64, Vec<i32>>, r_sketch:&HashMap<i64, Vec<i32>>, m:i32) -> (i32, i32, (i32, i32), (i32, i32)){    
    let mut g_add:Vec<DpCell> = Vec::new();  
    let mut g_not:Vec<DpCell> = Vec::new();  
    let mut anchors:Vec<(i32, i32)> = Vec::new();
    for (kmer, pos) in q_sketch{
        if r_sketch.contains_key(kmer){
            for i in pos{
                for j in r_sketch.get(kmer).expect("msg"){
                    anchors.push((*i, *j));
                }
            }            
        }
    }       
    anchors.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));      
    for i in 0.. anchors.len() as i32{
        let mut add_cell:DpCell = build_cell(1, 0, i, i);
        let mut not_cell:DpCell = build_cell(0, 0, -1, -1);
        if i == 0{            
            g_add.push(add_cell);
            g_not.push(not_cell);
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
                    //Determines optimal chain where anchor i is not included
                    if g_not[j].OptChain_score >= not_cell.OptChain_score{
                        not_cell = copy_cell(&g_not[j]);                 
                    }  
                    if g_add[j].OptChain_score >= not_cell.OptChain_score{
                        not_cell = copy_cell(&g_add[j]);                      
                    }                
                    //Determines optimal chain where anchor i is included
                    let prev = g_add[j].last_anchor;
                    if prev != -1{
                        let prev = prev as usize;
                        if (anchors[prev].1 < anchors[i as usize].1) && (anchors[prev].0 < anchors[i as usize].0){
                            let gap = ((anchors[prev].1 - anchors[i as usize].1) - (anchors[prev].0 - anchors[i as usize].0)).abs();
                            if g_add[j].OptChain_score + 1 >= add_cell.OptChain_score{
                                add_cell = build_cell(g_add[j].OptChain_score + 1, std::cmp::max(gap, g_add[j].opt_gap as i32), g_add[j].first_anchor, i);                        
                            }                       
                        }
                    }
                    let prev = g_not[j].last_anchor;
                    if prev != -1{
                        let prev = prev as usize;
                        if (anchors[prev].1 < anchors[i as usize].1) && (anchors[prev].0 < anchors[i as usize].0){
                            let gap = ((anchors[prev].1 - anchors[i as usize].1) - (anchors[prev].0 - anchors[i as usize].0)).abs();                 
                            if g_not[j].OptChain_score + 1 >= add_cell.OptChain_score{
                                add_cell = build_cell(g_not[j].OptChain_score + 1, std::cmp::max(gap, g_not[j].opt_gap as i32), g_not[j].first_anchor, i);                        
                            }                       
                        }
                    }
                }
                else{
                    unique_idx = m;
                }
            }
            
            g_add.push(add_cell);
            g_not.push(not_cell);
        }
    }    
    let ret_add = g_add.pop().expect("msg");
    let ret_not = g_not.pop().expect("msg");
    if ret_add.OptChain_score > ret_not.OptChain_score{
        let first_anchor = anchors[ret_add.first_anchor as usize];
        let last_anchor = anchors[ret_add.last_anchor as usize];
        (ret_add.OptChain_score, ret_add.opt_gap, first_anchor, last_anchor)
    }
    else{
        let first_anchor = anchors[ret_not.first_anchor as usize];
        let last_anchor = anchors[ret_not.last_anchor as usize];
        (ret_not.OptChain_score, ret_not.opt_gap, first_anchor, last_anchor)
    }
}