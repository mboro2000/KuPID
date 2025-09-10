use std::collections::HashMap;
use bio::io::fasta;
use std::str;

pub fn read_input(t:i32, file:&str) -> (Vec<HashMap<String, String>>){

    let mut iso_chunks:Vec<HashMap<String, String>> = Vec::new();
    let mut iso_data:HashMap<String, String> = HashMap::new();
    let mut i_count = 0;
    
    for i in 0 .. t{
        let mut iso_data:HashMap<String, String>= HashMap::new();
        iso_chunks.push(iso_data);
    } 

    let p = "gencode.v47.chr1_22.transcripts.fa";
    let read_results = fasta::Reader::from_file(file);
    for reader in read_results{
        for result in reader.records(){
        
            let record = result.expect("Error during fasta record parsing");
            let seq = str::from_utf8(record.seq()).unwrap().to_string();

            let i = (i_count % t) as usize;
            iso_chunks[i].entry(record.id().to_string()).or_insert(seq);     
            //iso_data.entry(record.id().to_string()).or_insert(seq);          
            i_count += 1;
        }
    }  

    (iso_chunks)
}