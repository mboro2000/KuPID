use std::io;
use std::collections::HashMap;
use std::collections::HashSet;

use std::fs::read_to_string;
use std::fs;
use std::path::Path;
use std::io::Write;
use bio::io::fasta;

//got hint about bam from chatGPT
use std::str;
use std::env;
use std::cmp;

struct Exon{
    name: String,
    start:usize,
    end:usize,
    length:usize
}
   
fn main() {   

    let high_k = 12;
    let low_k = 12;
    let w:f64 = 0.1 as f64;
    let t:i32 = 10;
    //let e = 0.014;
    let e = 0.05 as f64;

    let a = fastrand::i32(..);

    env::set_var("RUST_BACKTRACE", "1");

    iso_sets_new(high_k, low_k, w, e, a)

    //(same 0 bits btw old and maj) | (same 1 bits btw old and maj)

    //max advantage = # of neighbors in the set
    //max 7 neighbors => 3 bits

    //let advantage be 64 bit
    //each neighbor can have max length of 21 bits => 10 nucleotides

    /*
    PLAN
    ~ keep array for advantage bits

    ~ get the bits where old row = new row
    ~ apply to bits where old != new:
        if old != majority bits:
            adv += 2
        if old == majority bits:
            if adv = 1, flip the majority bit
            else: adv -= 2

    add 2 to adv for bit i:
        (1 << 1) << (i * (size for storing advantage bit))

    subtract 2 to adv for bit i:
        get mask value for the bit
       s = get opposite of mask val, add 2
       adv - orig val + (opposite of s)
    
     */
    
    /*
    let now = Instant::now();

    env::set_var("RUST_BACKTRACE", "1");

    iso_sets(high_k, low_k, w, e);

    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    */
}

fn read_annotation() -> (HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, String>, HashSet<String>){
    let mut exons_per_gene:HashMap<String, HashSet<String>> = HashMap::new();
    let mut isos_per_gene:HashMap<String, HashSet<String>> = HashMap::new();
    let mut exons_per_iso:HashMap<String, HashSet<String>> = HashMap::new();
    let mut exon_seqs:HashMap<String, String> = HashMap::new();
    let mut known_isoforms:HashSet<String> = HashSet::new();


    let mut iso_seqs:HashMap<String, String> = HashMap::new();

    let mut e_count = 0;


    let p = "gencode.v47.chr1.transcripts.fa";
    let read_results = fasta::Reader::from_file(p);
    for reader in read_results{
        for result in reader.records(){
        
            let record = result.expect("Error during fasta record parsing");
            let seq = str::from_utf8(record.seq()).unwrap();

            //let iso = record.id().split("|").collect::<Vec<&str>>()[0];
            let iso = record.id();
            iso_seqs.entry(iso.to_string()).or_insert(seq.to_string());
      
        }
    }  

    println!("Read exon seqs");

    let mut unique_exons:HashSet<String> = HashSet::new();

    for i in 0..1{
        let p = "human_p14_annotation.csv";
        let s = read_to_string(p).unwrap();
        let lines = s.lines();
        for line in lines{
            let cols = line.split(",").collect::<Vec<&str>>();
            
            let gene = cols[0].to_string();
            let iso = cols[1].to_string();
            let exon = cols[2].to_string();

            if gene != "gene"{                
                known_isoforms.insert(cols[1].to_string());

                exons_per_gene.entry(gene.clone()).or_insert(HashSet::new()).insert(exon.clone());
                exons_per_iso.entry(iso.clone()).or_insert(HashSet::new()).insert(exon.clone());
                isos_per_gene.entry(gene.clone()).or_insert(HashSet::new()).insert(iso.clone());

            }         
        }
    }

    println!("read annotation data");

    //(exons_per_gene, exons_per_iso, isos_per_gene, exon_seqs, known_isoforms)
    (exons_per_gene, exons_per_iso, isos_per_gene, iso_seqs, known_isoforms)

}


//Returns a unique number to represent the string *s
fn dna_encode(first:&str, last:&str, prev:i32, p0:i32, pw:i32) -> (i32, i32){

    //prev = score of previous kmer
    //pw = k-1
    //p0 = char from previous kmer that's being left behind
        
    let mut label:i32= 0;
    let mut p = -1;

    if "A" == first || "a" == first{
        p = 0;
    }
    else if "C" == first || "c" == first {
        p = 1;
    }
    else if "G" == first || "g" == first {
        p = 2;
    }
    else if "T" == first || "t" == first {
        p = 3;
    }


    label = prev - (p0 << (pw << 1));
    label = label << 2;

        
    if "A" == last || "a" == last{
        label += 0;
    }
    else if "C" == last || "c" == last {
        label +=  1;
    }
    else if "G" == last || "g" == last {
        label += 2;
    }
    else if "T" == last || "t" == last {
        label +=  3;
    }

    (label, p)
        
}



//Returns a table of minimizers selected using the FracMinHash method
//table = <score of minimizer m, set of positions where m is located in the sequence>
//fn FracMinHash(seq:&String, seq_len:i32, k:u32, s:f64, add_on:i32) -> HashMap<i32, HashSet<i32>> {
fn FracMinHash(seq:&str, k:i32, s:f64, a:i32) -> (HashSet<i32>) {
    let mut H = (1 << (2 * k-1));

    //let mut H =0x7FFFFFFF;
    let max = 0x7F;

    let Hs: f64 = H as f64 * s;

    let rounds = (seq.len() as i32 - k + 1);

    let mut table:HashSet<i32> = HashSet::new();


    let mut label:i32 = 0;

    let mut q : VecDeque<i8> = VecDeque::new();

    if seq.len() as i32 <= k{
        for item in seq.chars() {
            label  = label << 2;
    
            if 'A' == item || 'a' == item{
                label += 0;
            }
            else if 'C' == item || 'c' == item {
                label += 1;
            }
            else if 'G' == item || 'g' == item {
                label += 2;
            }
            else if 'T' == item || 't' == item {
                label += 3;
            }
        }
    
        let ax = a as i64 * label as i64;
        let ax_label = ax as i32 % max;
        if label as f64 <= Hs{
            table.insert(label);

            q.push_back(ax_label as i8);
        }
    }

    else{
        let s = &seq[0 .. k as usize];

        for item in s.chars() {
            label  = label << 2;

            if 'A' == item || 'a' == item{
                label += 0;
            }
            else if 'C' == item || 'c' == item {
                label += 1;
            }
            else if 'G' == item || 'g' == item {
                label += 2;
            }
            else if 'T' == item || 't' == item {
                label += 3;
            }
        }

        let ax = a as i64 * label as i64;
        let ax_label = ax as i32 % max;
        if label as f64 <= Hs{
            table.insert(label);

            q.push_back(ax_label as i8);
        }
        

        let first = &seq[0..1];
        let mut prev = label;
        let mut p0:i32 = 0;

        if "A" == first || "a" == first{
            p0 = 0;
        }
        else if "C" == first || "c" == first {
            p0 = 1;
        }
        else if "G" == first || "g" == first {
            p0 = 2;
        }
        else if "T" == first || "t" == first {
            p0 = 3;
        }

        for i in 1..rounds{
            let (score, p) = dna_encode(&seq[i as usize ..(i+1) as usize], &seq[(i+k-1) as usize .. (i+k) as usize], prev, p0, k-1);
            prev = score;
            p0 = p;

            let ax = a as i64 * score as i64;
            let ax_score = ax as i32 % max;
            if score as f64 <= Hs{
                table.insert(score);

                q.push_back(ax_score as i8);
            }
        }
    }

    (table)
    
}

fn create_iso_min(high_k:i32, low_k:i32, w:f64, a:i32) -> (HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>,
     HashSet<String>, HashMap<String, String>, HashMap<String, HashMap<String, HashSet<i64>>>, HashMap<i32, HashSet<String>>, HashMap<String, usize>){

    let n:usize = 32;

    let mut exon_min_data:HashMap<String, HashMap<String, HashSet<i64>>> = HashMap::new();

    //let (exons_per_gene, exons_per_iso, isos_per_gene, exon_seqs, known_isoforms) = read_annotation();
    let (exons_per_gene, exons_per_iso, isos_per_gene, iso_seqs, known_isoforms) = read_annotation();
    
    let mut h_kmer_map:HashMap<i32, HashSet<String>> = HashMap::new();

    let mut i_count = 0;

    let mut iso_info:HashMap<String, usize> = HashMap::new();

    let m_p = "ax_mod_values".to_string() + "_w" + &w.to_string() + ".csv";    
    let m = Path::new(&m_p);
    let mut m_line = "ax_mod\n".to_string();

    
    for iso in iso_seqs.keys(){
        let seq = iso_seqs.get(iso).expect("invalid iso key");

        i_count += 1;

        if i_count % 100 == 0{
            println!("{}", i_count);
        }

        let (high_min) = FracMinHash(seq, high_k, w, a);
        for kmer in high_min.iter(){
            h_kmer_map.entry(*kmer).or_insert(HashSet::new()).insert(iso.clone());

            let l = (*kmer).to_string();
            m_line.push_str(&l);

        }

        iso_info.insert(iso.to_string(), high_min.len());
    }

    
    fs::write(m, m_line);

    println!("read isoform minimizers");
    (exons_per_gene, exons_per_iso, isos_per_gene, known_isoforms, iso_seqs, exon_min_data, h_kmer_map, iso_info)

}

fn iso_sets_new(high_k:i32, low_k:i32, w:f64, e:f64, a:i32){

    let (exons_per_gene, exons_per_iso, isos_per_gene, known_isoforms, iso_seqs, exon_min_data, h_kmer_map, iso_info) = create_iso_min(high_k, low_k, w, a);

    let FP = 0.11;

    let mut num_novel = 0;

    let mut i_count = 0;

    let mut iso_found:HashSet<String> = HashSet::new();

    let mut total_reads = 0;

    let read_results = fasta::Reader::from_file("simulated_human_p14_as_reads.fa");

    
    //let c_p = "KMER_DETECT_annotated_weighted_PacSequel_high".to_string() + &high_k.to_string() + "_modax__w" + &w.to_string() + ".csv";
    //let d_p = "KMER_DETECT_novel_weighted_PacSequel_high".to_string() + &high_k.to_string() + "_modax_w" + &w.to_string() + ".csv";
    //let novel_p = "novel_candidates_PacSequel".to_string() + "_modax_w" + &w.to_string() + ".csv";
    
    //let c = Path::new(&c_p);
    //let d = Path::new(&d_p);
    //let novel = Path::new(&novel_p);

    /*
    let mut c_line = "jaccard\n".to_string();
    let mut d_line = "jaccard\n".to_string();
    let mut n_line = "".to_string();
     */
    
    let mut db_line = "distance\n".to_string();

    for reader in read_results{
        for result in reader.records(){
    
            let record = result.expect("Error during fasta record parsing") ;
            let seq = str::from_utf8(record.seq()).unwrap();            
            let id = record.id();

            total_reads += 1;

            let iso_orig = id.split("=").collect::<Vec<&str>>()[1];
        
            if i_count % 100 == 0{
                println!("Reads {}", i_count);
            }

            i_count += 1;
                
            let (high_k_table) = FracMinHash(seq, high_k, w, a);
            let read_len = high_k_table.len() as f32;

            iso_found.insert(iso_orig.to_string());

            if !known_isoforms.contains(iso_orig){
               num_novel += 1;
            }

            let mut iso_occur:HashMap<String, i32> = HashMap::new();
            
            for min in high_k_table.iter(){
                if h_kmer_map.contains_key(min){

                    for iso in h_kmer_map.get(min).expect("msg"){
                        *iso_occur.entry(iso.clone()).or_insert(0) += 1;
                    }
                }
            }

            let mut best_weighted_jaccard:f32 = 0.0;
            
            for (iso, num_match) in iso_occur.iter(){
                let iso_len = *iso_info.get(iso).expect("msg") as f32;
                let iso_inter = *num_match as f32;        // = num_match

                let mut cov = 0.0;
                if read_len < iso_len{
                    cov = read_len / iso_len;
                }
                else{
                    cov = iso_len / read_len;
                }

                let wj = cov * iso_inter / (read_len + iso_len - iso_inter);

                if wj > best_weighted_jaccard{
                    best_weighted_jaccard = wj;
                }
            }

            /* 
            let l = best_weighted_jaccard.to_string() + "\n"; 
                if known_isoforms.contains(iso_orig){      
                    c_line.push_str(&l);
                }
                else{
                    d_line.push_str(&l);
                }

            if best_weighted_jaccard <= 0.33{
                let l = ">".to_string() + id + "\n" + seq + "\n";
                n_line.push_str(&l);
            }
            */
            
        }
    }

    println!("{:#?}", known_isoforms.len());

    //println!("Reads Processed: {}", r_count);
    println!("Total Reads: {}", total_reads);
    
    //fs::write(novel, novel_line);

    /*
    fs::write(c, c_line);
    fs::write(d, d_line);
    fs::write(novel, n_line);
     */
    

    println!("Number of Novel Transcripts: {}", num_novel);
    println!("Number of isoforms in sample reads: {}", iso_found.len());

}
