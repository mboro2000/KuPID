use std::collections::HashMap;

#[derive(Clone)]
pub struct Sketch{
    pub size:i32,
    pub kmers:HashMap<i64, Vec<i32>>,
    pub seq_len:usize
}

#[derive(Clone)]
pub struct Match{
    pub sample_id:String,
    pub similarity:f32,
    pub chunk:usize,
    pub num_matches:f32,
    pub query_len:usize
}

pub fn build_sketch(size:i32, kmers:HashMap<i64, Vec<i32>>, seq_len:usize) -> Sketch{
    Sketch{
        size:size,
        kmers:kmers,
        seq_len:seq_len
    }
}

pub fn build_match(sample_id:String, similarity:f32, chunk:usize, num_matches:f32, query_len:usize) -> Match{
    Match {
        sample_id: sample_id,
        similarity: similarity,
        chunk:chunk,
    num_matches:num_matches,
query_len:query_len}
}

