use std::collections::VecDeque;
use std::hash::Hash;
use std::io;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufRead;
use std::thread;
use std::fs::read_to_string;
use std::fs;
use std::path::Path;
use std::io::Write;
use bio::io::fasta;

//got hint about bam from chatGPT
use std::str;
use std::env;

struct IsoGroup {
    name: i32,
    gene: String,
    exons: HashSet<String>,
    reads: Vec<Read>,
    length: i32,
   }

struct Read{
    name: String,
    length: i32,
    num_matches: HashMap<String, i32>,
    num_low_k: i32,
    jaccard:f32,
    weighted_jaccard:f32
}

struct Exon{
    name: String,
    start:usize,
    end:usize,
    length:usize
}

struct Overlap{
    max_start:i32,
    max_end:i32,
    exons:Vec<Exon>
}
    

fn build_iso_group(name: i32, gene:String, exons: HashSet<String>, reads: Vec<Read>, length: i32) -> IsoGroup {
    IsoGroup {
        name: name,
        gene:gene,
        exons: exons,
        reads: reads,
        length: length,
    }
}

fn build_read(name: String, length: i32, num_matches: HashMap<String, i32>, num_low_k: i32, jaccard:f32, weighted_jaccard:f32) -> Read{
    Read { name: name, length: length, num_matches: num_matches, num_low_k: num_low_k, jaccard, weighted_jaccard }
}

fn main() {
    let high_k = 20;
    let low_k = 12;
    let w:f64 = 0.2 as f64;
    let t:i32 = 10;
    let e = 0.014;
            
    env::set_var("RUST_BACKTRACE", "1");

    iso_sets(high_k, low_k, w, e);

}

fn read_annotation() -> (HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, String>, HashSet<String>){
    let mut exons_per_gene:HashMap<String, HashSet<String>> = HashMap::new();
    let mut isos_per_gene:HashMap<String, HashSet<String>> = HashMap::new();
    let mut exons_per_iso:HashMap<String, HashSet<String>> = HashMap::new();
    let mut exon_seqs:HashMap<String, String> = HashMap::new();
    let mut known_isoforms:HashSet<String> = HashSet::new();

    let mut e_per_g:HashMap<String, Vec<Exon>> = HashMap::new();

    let mut e_count = 0;

    for i in 0..1{
    //for i in 0..10{
        //let p =  "mouse_reduced_exons_".to_string() + &i.to_string() + ".fasta";
        let p = "sub_vM25_exons.fasta";
        let read_results = fasta::Reader::from_file(p);
        for reader in read_results{
            for result in reader.records(){

                println!("{}", e_count);
                e_count += 1;
            
                let record = result.expect("Error during fasta record parsing");
                let seq = str::from_utf8(record.seq()).unwrap();
                exon_seqs.entry(record.id().to_string()).or_insert(seq.to_string());
          
            }
        }   
    }

    println!("Read exon seqs");

    let mut unique_exons:HashSet<String> = HashSet::new();

    for i in 0..1{
    //for i in 0..10{
        //let p =  "mus_annot_".to_string() + &i.to_string() + ".csv";
        let p = "sub_known_vM25_annotation.csv";
        let s = read_to_string(p).unwrap();
        let lines = s.lines();
        for line in lines{
            let cols = line.split(",").collect::<Vec<&str>>();
            
            let gene = cols[0].to_string();
            let iso = cols[1].to_string();
            let exon = cols[2].to_string();

            if gene != "gene"{
                let start:i32 = cols[3].parse().expect("invalid");
                let end:i32 = cols[4].parse().expect("invalid");
                
                known_isoforms.insert(cols[1].to_string());

                exons_per_gene.entry(gene.clone()).or_insert(HashSet::new()).insert(exon.clone());
                exons_per_iso.entry(iso.clone()).or_insert(HashSet::new()).insert(exon.clone());
                isos_per_gene.entry(gene.clone()).or_insert(HashSet::new()).insert(iso.clone());

                /*
                if ! unique_exons.contains(&exon){
                    e_per_g.entry(gene).or_insert(Vec::new()).push(Exon{name:exon.to_string(), start:start, end:end, length:end-start});
                    unique_exons.insert(exon);
                    exons_per_gene.entry(gene).or_insert(HashSet::new()).insert(exon);
                }
                */
            }         
        }
    } 

    let num_genes = e_per_g.keys().len();

    /*
    let mut avg_exons = 0;
    for key in e_per_g.keys(){
        let c = e_per_g.get(key).expect("msg");
        println!("{}", e_per_g.get(key).expect("msg").len());
        avg_exons += e_per_g.get(key).expect("msg").len();
    }

    let mut exon_len = 0;

    for exon in unique_exons.iter(){
        exon_len += exon_seqs.get(exon).expect("msg").len();
    }

    println!("Sum of exons per gene: {}", avg_exons as f64);
    println!("Average exons per gene: {}", avg_exons as f64 / num_genes as f64);
    println!("Avg exon length: {}", exon_len as f64 / unique_exons.len() as f64);
    */
    //

    println!("read annotation data");

    /*
    for (gene, exon_data) in e_per_g.iter_mut(){
        let mut overlaps:Vec<Overlap> = Vec::new();
        
        exon_data.sort_by(|a, b| b.length.cmp(&a.length));      

        for exon in exon_data{

            let mut overlapped = 0;

            if overlaps.len() == 0{
                let mut exons:Vec<Exon> = Vec::new();
                exons.push(Exon{name:exon.name.to_string(), start:exon.start, end:exon.end, length:exon.length});
                overlaps.push(Overlap{max_start:exon.start, max_end:exon.end, exons:exons});
            }

            else{
                for i in 0..overlaps.len(){
                    let o_s = overlaps[i].max_start;
                    let o_e = overlaps[i].max_end;

                    if exon.start == o_s && exon.end == o_e{
                        overlapped += 1;
                        overlaps[i].exons.push(Exon{name:exon.name.to_string(), start:exon.start, end:exon.end, length:exon.length});
                    }
                    
                    else if exon.start == o_s && exon.end != o_e{
                        overlapped += 1;
                        if exon.end > o_e{
                            overlaps[i].max_end = exon.end;
                        }
                        overlaps[i].exons.push(Exon{name:exon.name.to_string(), start:exon.start, end:exon.end, length:exon.length});
                    }

                    else if exon.end == o_e && exon.start != o_s{
                        overlapped += 1;
                        if exon.start < o_e{
                            overlaps[i].max_end = exon.end;
                        }
                        overlaps[i].exons.push(Exon{name:exon.name.to_string(), start:exon.start, end:exon.end, length:exon.length});                  
                    }

                    else if exon.start > o_s && exon.start < o_e {
                        overlapped += 1;
                        overlaps[i].exons.push(Exon{name:exon.name.to_string(), start:exon.start, end:exon.end, length:exon.length});
                        if exon.end > o_e{
                            overlaps[i].max_end = exon.end;
                        }
                    }

                    else if exon.end > o_s && exon.end < o_e{
                        overlapped += 1;
                        overlaps[i].exons.push(Exon{name:exon.name.to_string(), start:exon.start, end:exon.end, length:exon.length});
                        if exon.start < o_e{
                            overlaps[i].max_end = exon.end;
                        }  
                    }
                }
                if overlapped == 0{
                    let mut exons:Vec<Exon> = Vec::new();
                    exons.push(Exon{name:exon.name.to_string(), start:exon.start, end:exon.end, length:exon.length});
                    overlaps.push(Overlap{max_start:exon.start, max_end:exon.end, exons:exons});
                }            
            }
        }

        for overlap in overlaps{
            if overlap.exons.len() == 1{
                for exon in overlap.exons{
                    if exon.length >= low_k{
                        exons_per_gene.entry(gene.to_string()).or_insert(HashSet::new()).insert(exon.name);
                    }
                }
            }
            else{
                let mut max_len = 0;
                let mut max_exon = "".to_string();
                for exon in overlap.exons{
                    if exon.length > max_len{
                        max_len = exon.length;
                        max_exon = exon.name.to_string();
                    }
                }
                if max_len >= low_k{
                    exons_per_gene.entry(gene.to_string()).or_insert(HashSet::new()).insert(max_exon.to_string());
                }                
            }
        }
    }
    println!("Removed overlaps");
    */

    (exons_per_gene, exons_per_iso, isos_per_gene, exon_seqs, known_isoforms)

}


//Returns a unique number to represent the string *s
fn dna_encode(s:&str, k:i32, prev:i64, p0:i64, pw:i64) -> (i64, i64){
        
    let mut label:i64 = 0;
    let mut p = -1 as i64;

    let first_char = s.chars().nth(0).expect("no char");

    if 'A' == first_char || 'a' == first_char{
        p = 0;
    }
    else if 'C' == first_char || 'c' == first_char {
        p = 1;
    }
    else if 'G' == first_char || 'g' == first_char {
        p = 2;
    }
    else if 'T' == first_char || 't' == first_char {
        p = 3;
    }

    if prev == -1{
        let mut enc: Vec<i64> = Vec::new();

        for item in s.chars() {
            if 'A' == item || 'a' == item{
                enc.push(0);
            }
            else if 'C' == item || 'c' == item {
                enc.push(1);
            }
            else if 'G' == item || 'g' == item {
                enc.push(2);
            }
            else if 'T' == item || 't' == item {
                enc.push(3);
            }
        }

        for i in 0..k{
            let c_val = &enc[i as usize];
            label *= 4;
            label += *c_val;
        }
    }
    
    else{
        let last_char = s.chars().nth(s.len() - 1).expect("no char");
        label = prev - (pw * p0);
        label *= 4;
        
        if 'A' == last_char || 'a' == last_char{
            label += 0;
        }
        else if 'C' == last_char || 'c' == last_char {
            label +=  1;
        }
        else if 'G' == last_char || 'g' == last_char {
            label += 2;
        }
        else if 'T' == last_char || 't' == last_char {
            label +=  3;
        }
    }

    (label, p)
        
}

//Returns a table of minimizers selected using the FracMinHash method
//table = <score of minimizer m, set of positions where m is located in the sequence>
//fn FracMinHash(seq:&String, seq_len:i32, k:u32, s:f64, add_on:i32) -> HashMap<i32, HashSet<i32>> {
fn FracMinHash(seq:&str, k:i32, s:f64) -> HashSet<i64> {
    let mut Hs = 1 as f64;
    for i in 0..k{
        Hs *= 4 as f64;
    }
    Hs *= s;
    Hs -= 1 as f64;

    let rounds = (seq.len() as i32 + 1 - k);

    let mut table:HashSet<i64> = HashSet::new();

    let mut prev = -1 as i64;
    let mut p0 = -1 as i64;

    let mut pw = 1 as i64;
    for i in 0..k-1{
        pw *= 4;
    }


    for i in 0..rounds{
        let a = i as usize;
        let b = a+k as usize;
        let sub = &seq[a..b];

        let (score, p) = dna_encode(sub, k as i32, prev, p0, pw);
        prev = score;
        p0 = p;

        if score as f64 <= Hs{
            table.insert(score);
        }
       
    }
    table
}

fn create_iso_min(high_k:i32, low_k:i32, w:f64) -> (HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashSet<String>, HashMap<String, String>, HashMap<String, HashMap<String, HashSet<i64>>>, HashMap<i64, HashSet<String>>){

    let n:usize = 32;

    let mut exon_min_data:HashMap<String, HashMap<String, HashSet<i64>>> = HashMap::new();

    let (exons_per_gene, exons_per_iso, isos_per_gene, exon_seqs, known_isoforms) = read_annotation();
    
    let mut h_kmer_map:HashMap<i64, HashSet<String>> = HashMap::new();
    let mut e_count = 0;

    for gene in exons_per_gene.keys(){

        let exons = exons_per_gene.get(gene).expect("Incorrect key").iter();
        
        for exon in exons{
            if exon != "exon"{

                if e_count % 1000 == 0{
                    println!("{}", e_count);
                }

                e_count += 1;

                let seq = exon_seqs.get(exon).expect("invalid exon key");
                let mut min_data:HashMap<String, HashSet<i64>> = HashMap::new();
                let high_min = FracMinHash(seq, high_k, w);
                
                for kmer in high_min.iter(){
                    h_kmer_map.entry(*kmer).or_insert(HashSet::new()).insert(gene.clone());
                }

                min_data.insert("high k".to_string(), high_min);
                min_data.insert("low k".to_string(), FracMinHash(seq, low_k, w));              
        
                exon_min_data.insert(exon.clone(), min_data);
            }
        }

    }

    println!("read exon minimizers");
    (exons_per_gene, exons_per_iso, isos_per_gene, known_isoforms, exon_seqs, exon_min_data, h_kmer_map)

}

fn iso_sets(high_k:i32, low_k:i32, w:f64, e:f64){

    let (exons_per_gene, exons_per_iso, isos_per_gene, known_isoforms, exon_seqs, exon_min_data, h_kmer_map) = create_iso_min(high_k, low_k, w);

    let FP = 0.11;

    let mut num_novel = 0;
    let mut correct_contain_ref:Vec<f64> = Vec::new();
    let mut correct_contain_query:Vec<f64> = Vec::new();
    let mut novel_contain_ref:Vec<f64> = Vec::new();
    let mut novel_contain_query:Vec<f64> = Vec::new();

    let mut correct_weighted_contain_ref:Vec<f64> = Vec::new();
    let mut correct_weighted_contain_query:Vec<f64> = Vec::new();
    let mut novel_weighted_contain_ref:Vec<f64> = Vec::new();
    let mut novel_weighted_contain_query:Vec<f64> = Vec::new();

    let mut orig_gene_matches:Vec<i32> = Vec::new();
    let mut false_gene_matches:Vec<i32> = Vec::new();

    let mut correct_jaccard:Vec<f32> = Vec::new();
    let mut novel_jaccard:Vec<f32> = Vec::new();

    let mut correct_weighted_jaccard:Vec<f32> = Vec::new();
    let mut novel_weighted_jaccard:Vec<f32> = Vec::new();

    let mut iso_groups:Vec<IsoGroup> = Vec::new();
    let mut group_count = 0;

    let mut i_count = 0;

    let mut iso_found:HashSet<String> = HashSet::new();

    let mut total_reads = 0;
    let mut no_gene_match = 0;
    let mut novel_wo_gene = 0;

    

    let novel_p = "novel_candidates_PacSequel.csv";
    let novel = Path::new(&novel_p);

    let mut novel_line = "".to_string();


    let mut novel_found:HashSet<String> = HashSet::new();
    /*
    
    let new = orig_gene_matches[i].to_string() + "\n";
        l += &new;        
        
        a_line.push_str(&l);

     */



    let read_results = fasta::Reader::from_file("sim_vM25_reads.fa");
    //let read_results = fasta::Reader::from_file("src/sim_reads.fasta");
    for reader in read_results{
        for result in reader.records(){
    
            let record = result.expect("Error during fasta record parsing") ;
            let seq = str::from_utf8(record.seq()).unwrap();            
            let id = record.id();

            total_reads += 1;

            //let iso = id.split("_").collect::<Vec<&str>>()[0];
            let iso = id.split("_").collect::<Vec<&str>>()[1];
        
            if i_count % 100 == 0{
                println!("{}", i_count);
            }

            i_count += 1;
            let mut min_match = 0;
            let mut match_gene = "";
                
            let high_k_table = FracMinHash(seq, high_k, w);
            let low_k_table = FracMinHash(seq, low_k, w);

            iso_found.insert(iso.to_string());

            if !known_isoforms.contains(iso){
               num_novel += 1;
            }

            let mut pot_genes:HashSet<String> = HashSet::new();
            
            for min in high_k_table.iter(){
                if h_kmer_map.contains_key(min){
                    for gene in h_kmer_map.get(min).expect("msg"){
                        pot_genes.insert(gene.clone());
                    }
                }
            }

            let mut g = 0;

            for gene in pot_genes.iter(){
                
                g += 1;

                let mut kmer_match = 0;
                let mut gene_min:HashSet<i64> = HashSet::new();

                let exons = exons_per_gene.get(gene).expect("Incorrect key").iter();
                
                for exon in exons{
                    if exon != "exon"{
                        let exon_keys = exon_min_data.get(exon).expect("incorrect key").get("high k").expect("incorrect key").iter();
                        for key in exon_keys{
                            gene_min.insert(*key);
                        }
                    }
                }

                let inter = gene_min.intersection(&high_k_table);
                let inter_len = inter.count() as i32;
                kmer_match += inter_len;

                if isos_per_gene.get(gene).expect("msg").contains(iso){
                    orig_gene_matches.push(kmer_match.clone());
                }
                else{
                    false_gene_matches.push(kmer_match.clone());
                }

                if inter_len >= min_match{
                    match_gene = gene;
                    min_match = kmer_match;
                }
            }

            if match_gene == ""{
                no_gene_match += 1;

                if !known_isoforms.contains(iso){
                    novel_wo_gene += 1;
                    novel_found.insert(iso.to_string());
                }            

                


                let line = ">".to_string() + id + "\n" + seq + "\n";

                    novel_line.push_str(&line);

            }

            if match_gene != ""{
                let mut e_matches:HashMap<String, i32> = HashMap::new();
                let mut exon_set:HashSet<String> = HashSet::new();

                let mut match_length = 0;
                let mut group_length = 0;
                let mut group_min = 0;
            
                for exon in exons_per_gene.get(match_gene).expect("Incorrect key").iter(){

                    if exon != "exon"{
                        let exon_mins = exon_min_data.get(exon).expect("incorrect key").get("low k").expect("incorrect key");
                        let inter = exon_mins.intersection(&low_k_table);
                        let inter_len = inter.count() as i32;

                        let FP_match = exon_mins.len() as f32 * FP;

                        if inter_len as f32 > FP_match{
                            e_matches.entry(exon.to_string()).or_insert(inter_len);
                            exon_set.insert(exon.to_string());
                            match_length += inter_len;
                            group_length += exon_seqs.get(exon).expect("incorrect key").len();
                            group_min += exon_mins.len();
                        }                      
                    }
                }
            
                let l_len = low_k_table.len() as i32;
                let mut group_match = 0;

                let mut best_jaccard:f32 = 0.0;
                let mut best_weighted_jaccard:f32 = 0.0;


                for iso in isos_per_gene.get(match_gene).expect("invalid key"){
                    let mut iso_kmers:HashSet<i64> = HashSet::new();

                    for exon in exons_per_iso.get(iso).expect("msg"){
                        for kmer in exon_min_data.get(exon).expect("msg").get("low k").expect("msg"){
                            iso_kmers.insert(kmer.clone());
                        }
                    }

                    let iso_len = iso_kmers.len() as i32 + low_k;
                    let cov = seq.len() as f32 / iso_len as f32;

                    let jaccard = iso_kmers.intersection(&low_k_table).count() as f32 / iso_kmers.union(&low_k_table).count() as f32;
                    let weighted_jaccard = jaccard * cov;

                    if jaccard > best_jaccard{
                        best_jaccard = jaccard;
                    }
                    if weighted_jaccard > best_weighted_jaccard{
                        best_weighted_jaccard = weighted_jaccard;
                    }
                }

                if best_weighted_jaccard < 1.5{
                    let line = ">".to_string() + id + "\n" + seq + "\n";

                    novel_line.push_str(&line);

                    if !known_isoforms.contains(iso){
                        novel_found.insert(iso.to_string());
                    }

                    
                }

                if iso_groups.len() == 0{
                    let mut reads:Vec<Read> = Vec::new();
                    reads.push(build_read(id.to_string(), seq.len() as i32, e_matches.clone(), l_len, best_jaccard, best_weighted_jaccard));
                    iso_groups.push(build_iso_group(group_count, match_gene.to_string(), exon_set.clone(), reads, group_length as i32));
                    group_count += 1;
                }

                else{
                    for i in 0..iso_groups.len(){

                        if exon_set.clone().difference(&iso_groups[i].exons).count() == 0 && (iso_groups[i].exons).difference(&exon_set).count() == 0{
                            iso_groups[i].reads.push(build_read(id.to_string(), seq.len() as i32, e_matches.clone(), l_len, best_jaccard, best_weighted_jaccard));
                            group_match += 1;
                        }
                    }

                    if group_match == 0{
                        let mut reads:Vec<Read> = Vec::new();
                        reads.push(build_read(id.to_string(), seq.len() as i32, e_matches.clone(), l_len, best_jaccard, best_weighted_jaccard));
                        iso_groups.push(build_iso_group(group_count, match_gene.to_string(), exon_set.clone(), reads, group_length as i32));
                        group_count += 1;
                    }  
                }
            }
        }
    }

    println!("{:#?}", known_isoforms);

    let mut r_count = 0;

    for i in 0..iso_groups.len(){

        //if iso_groups[i].reads.len() > 2 && iso_groups[i].exons.len() > 0{
        if iso_groups[i].reads.len() > 0 {

            let mut best_ref_iso = "";
            let mut best_ref_con = 0.0;
            let mut num_best_ref:usize = 0;
            let mut best_query_iso = "";
            let mut best_query_con = 0.0;

            let isos = isos_per_gene.get(&iso_groups[i].gene).expect("msg");
            for iso in isos{
                let ref_exons = exons_per_iso.get(iso).expect("msg");

                let inter = iso_groups[i].exons.intersection(ref_exons).count();

                let ref_con = inter as f32 / ref_exons.len() as f32;
                let query_con = inter as f32 / iso_groups[i].exons.len() as f32;


                if ref_con > best_ref_con{
                    best_ref_con = ref_con;
                    best_ref_iso = iso;
                    num_best_ref = ref_exons.len();
                }

                if ref_con == best_ref_con{
                    if ref_exons.len() > num_best_ref{
                        best_ref_con = ref_con;
                        best_ref_iso = iso;
                        num_best_ref = ref_exons.len();
                    }
                }

                if query_con > best_query_con{
                    best_query_con = query_con;
                    best_query_iso = iso;
                }
            }

            let mut num_q_iso_min:usize = 0;
            let mut num_r_iso_min:usize = 0;

            //set of exons from the isoform with the highest proportion of the read's exons

            if best_query_iso != ""{
                let best_query_exons = exons_per_iso.get(best_query_iso).expect("msg");

                for exon in best_query_exons{
                    num_q_iso_min += exon_min_data.get(exon).expect("incorrect key").get("low k").expect("incorrect key").len();
                }
            }
            
            if best_ref_iso != ""{
                let best_ref_exons = exons_per_iso.get(best_ref_iso).expect("msg");

                for exon in best_ref_exons{
                    num_r_iso_min += exon_min_data.get(exon).expect("incorrect key").get("low k").expect("incorrect key").len();
                }
            }
            
            let q_iso_len = num_q_iso_min as f64 + low_k as f64;
            let r_iso_len = num_r_iso_min as f64 + low_k as f64;

            for read in iso_groups[i].reads.iter(){

                r_count += 1;

                let iso = (*read).name.split("_").collect::<Vec<&str>>()[1];
                //let iso = (*read).name.split("_").collect::<Vec<&str>>()[0];
                                
                let mut con_ref = 0;
                let mut con_query = 0;

                if best_query_iso != ""{
                    for exon in exons_per_iso.get(best_query_iso).expect("msg"){
                        
                        if read.num_matches.contains_key(exon){
                            con_ref += read.num_matches.get(exon).expect("msg");            
                        }
                                        
                    }  
                }    

                if best_ref_iso != ""{

                    for exon in exons_per_iso.get(best_ref_iso).expect("msg"){
                        if read.num_matches.contains_key(exon){
                            con_query += read.num_matches.get(exon).expect("msg");
                        }                    
                    }
                }
                
                if !known_isoforms.contains(iso){
                    novel_contain_ref.push(con_ref as f64 / num_q_iso_min as f64);
                    novel_contain_query.push(con_query as f64 / read.num_low_k as f64);

                    novel_weighted_contain_ref.push((con_ref as f64 * q_iso_len) / num_q_iso_min as f64);
                    novel_weighted_contain_query.push((con_query as f64 * r_iso_len) / read.num_low_k as f64);


                    novel_jaccard.push(read.jaccard);
                    novel_weighted_jaccard.push(read.weighted_jaccard);
                }
                
                if known_isoforms.contains(iso){
                    correct_contain_ref.push(con_ref as f64 / num_q_iso_min as f64);
                    correct_contain_query.push(con_query as f64 / read.num_low_k as f64);

                    correct_weighted_contain_ref.push((con_ref as f64 * q_iso_len) / num_q_iso_min as f64);
                    correct_weighted_contain_query.push((con_query as f64 * r_iso_len) / read.num_low_k as f64);

                    correct_jaccard.push(read.jaccard);
                    correct_weighted_jaccard.push(read.weighted_jaccard);
                }
            }
        }
    }

    println!("Reads Processed: {}", r_count);
    println!("Total Reads: {}", total_reads);
    println!("Reads with no gene match: {}", no_gene_match);
    println!("Novel reads with no gene match: {}", novel_wo_gene);

    let a_p = "true_gene_matches_P4C2_high".to_string() + &high_k.to_string() + "_w" + &w.to_string() + ".csv";
    let b_p = "false_gene_matches_P4C2_high".to_string() + &high_k.to_string() + "_w" + &w.to_string() + ".csv";
    let c_p = "KMER_DETECT_annotated_weighted_P5C3_containment_filtered_high".to_string() + &high_k.to_string() + "_w" + &w.to_string() + ".csv";
    let d_p = "KMER_DETECT_novel_containment_P5C3_weighted_filtered_high".to_string() + &high_k.to_string() + "_w" + &w.to_string() + ".csv";

    let a = Path::new(&a_p);
    let b = Path::new(&b_p);
    let c = Path::new(&c_p);
    let d = Path::new(&d_p);

    let mut a_line = "true gene\n".to_string();
    let mut b_line = "false gene\n".to_string();
    let mut c_line = "jaccard,two way contain\n".to_string();
    let mut d_line = "jaccard,two way contain\n".to_string();

    for i in 0..orig_gene_matches.len(){
        let mut l = "".to_string();
        
        let new = orig_gene_matches[i].to_string() + "\n";
        l += &new;        
        
        a_line.push_str(&l);
    }

    for i in 0..false_gene_matches.len(){
        let mut l = "".to_string();
        
        let new = false_gene_matches[i].to_string() + "\n";
        l += &new;        
        
        b_line.push_str(&l);
    }

    for i in 0..correct_weighted_contain_ref.len(){
        let mut l = "".to_string();
        if correct_weighted_contain_ref[i] < correct_weighted_contain_query[i]{
            let new = correct_weighted_jaccard[i].to_string() + "," + &correct_weighted_contain_ref[i].to_string() + "\n";
            l += &new;
        }
        else{
            let new = correct_weighted_jaccard[i].to_string() + "," + &correct_weighted_contain_query[i].to_string() + "\n";
            l += &new;
        }
        c_line.push_str(&l);
    }

    for i in 0..novel_weighted_contain_ref.len(){
        let mut l = "".to_string();
        if novel_weighted_contain_ref[i] < novel_weighted_contain_query[i]{
            let new = novel_weighted_jaccard[i].to_string() + "," + &novel_weighted_contain_ref[i].to_string() + "\n";
            l += &new;
        }
        else{
            let new = novel_weighted_jaccard[i].to_string() + "," + &novel_weighted_contain_query[i].to_string() + "\n";
            l += &new;
        }

        d_line.push_str(&l);
    }

    
    fs::write(novel, novel_line);

    /*
    fs::write(a, a_line);
    fs::write(b, b_line);
    fs::write(c, c_line);
    fs::write(d, d_line);
    */    

    println!("Number of Novel Transcripts: {}", num_novel);
    println!("Number of isoforms in sample reads: {}", iso_found.len());
    println!("Number of isoforms found in selected candidates: {}", novel_found.len())

}
