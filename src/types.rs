use std::collections::HashMap;

#[derive(Clone)]
pub struct Sketch{
    pub size:i32,
    pub kmers:HashMap<i64, Vec<i32>>,
    pub seq_len:usize
}

pub fn build_sketch(size:i32, kmers:HashMap<i64, Vec<i32>>, seq_len:usize) -> Sketch{
    Sketch{
        size:size,
        kmers:kmers,
        seq_len:seq_len
    }
}

#[derive(Clone)]
pub struct Match{
    pub sample_id:String,
    pub similarity:f32,
    pub chunk:usize
}

pub fn build_match(sample_id:String, similarity:f32, chunk:usize) -> Match{
    Match {
        sample_id: sample_id,
        similarity: similarity,
        chunk:chunk}
}

#[derive(Clone)]
pub struct DpCell{
    pub OptChain_score:i32,
    pub opt_gap:i32,
    pub first_anchor: i32,
    pub last_anchor: i32
}

pub fn build_cell(OptChain_score:i32, opt_gap:i32, first_anchor:i32, last_anchor:i32) -> DpCell{
    DpCell { OptChain_score: OptChain_score, opt_gap: opt_gap, first_anchor: first_anchor, last_anchor: last_anchor }
}

pub fn copy_cell(cell: &DpCell) -> DpCell{
    DpCell { OptChain_score: cell.OptChain_score, opt_gap: cell.opt_gap, first_anchor: cell.first_anchor, last_anchor: cell.last_anchor }
}

pub struct OptChain{
    pub score:i32,
    pub ref_match:String,
    pub gap:i32,
    pub query_gap:i32,
    pub ref_gap:i32,
    pub similarity:f32
}

pub fn init_chain() -> OptChain{
    OptChain { score: 0, ref_match: "".to_string(), gap: 0, query_gap: 0, ref_gap: 0, similarity: -1.0 * f32::INFINITY}
}

pub fn build_chain(score:i32, ref_match:String, gap:i32, query_gap:i32, ref_gap:i32, similarity:f32) -> OptChain{
    OptChain { score: score, ref_match: ref_match, gap: gap, query_gap: query_gap, ref_gap: ref_gap, similarity: similarity}
}