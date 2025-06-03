use std::collections::HashMap;

#[derive(Clone)]
pub struct Sketch{
    pub id:String,
    pub size:i32,
    pub kmers:HashMap<i64, Vec<i32>>
}

#[derive(Clone)]
pub struct Match{
    pub sample_id:String,
    pub similarity:f32,
    pub gap:i32,
    pub chunk:usize
}

pub fn build_Sketch(id:String, size:i32, kmers:HashMap<i64, Vec<i32>>) -> Sketch{
    Sketch{
        id:id,
        size:size,
        kmers:kmers
    }
}

pub fn build_Match(sample_id:String, similarity:f32, gap:i32, chunk:usize) -> Match{
    Match {
        sample_id: sample_id,
        similarity: similarity, 
        gap: gap,
        chunk:chunk}
}

pub struct Bin{
    pub threshold: f32,
    pub matches: Vec<Match>,
    pub contain_novel: bool,
    pub min_val: f32,
    pub num_novel: i32,
    pub num_total: i32,
    pub scores: Vec<i32>,
}

pub fn build_Bin(threshold:f32, matches:Vec<Match>, contain_novel:bool, min_val:f32, num_novel:i32, num_total:i32, scores:Vec<i32>) -> Bin{
    Bin { 
        threshold: threshold,
        matches: matches,
        contain_novel: contain_novel,
        min_val: min_val,
        num_novel: num_novel,
        num_total: num_total,
        scores:scores
    }
}
