use std::collections::HashMap;
use std::error::Error;
use csv::ReaderBuilder;

fn main(){
    let AS_params:HashMap<String, f32> = HashMap::from([("ES", 0.4),("IR", 0.1),("5p", 0.25),("3p", 0.25)]);

    let file_path = "gencode.vM25.annotation.gtf";
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t') // GTF files are tab-delimited
        .has_headers(false) // GTF files typically don't have headers
        .from_path(file_path)?;

    for result in reader.records() {
        let record = result?;
        println!("{:?}", record);
    }
}
