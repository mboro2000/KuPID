//Returns a unique number to represent the string *s
fn dna_encode(s:&str, k:i32) -> i32{
    let mut enc: Vec<i32> = Vec::new();;
    let mut label = 0;
    
    for item in s.chars() {
        if 'A' == item {
            enc.push(0);
        }
        else if 'a' == item {
            enc.push(0);
        }
        else if 'C' == item {
            enc.push(1);
        }
        else if 'c' == item {
            enc.push(1);
        }
        else if 'G' == item {
            enc.push(2);
        }
        else if 'g' == item {
            enc.push(2);
        }
        else if 'T' == item {
            enc.push(3);
        }
        else if 't' == item {
            enc.push(3);
        }
    }

    let FOUR = 4 as i32;
    for i in 0..k{
        let c_val = &enc[i as usize];
        let power = FOUR.pow(i as u32);
        label += *c_val * power;
    }

    label
        
}

//Returns a table of minimizers selected using the FracMinHash method
//table = <score of minimizer m, set of positions where m is located in the sequence>
fn FracMinHash(seq:&String, seq_len:i32, k:u32, s:f64, add_on:i32) -> HashMap<i32, HashSet<i32>> {
    let FOUR = 4 as i32;

    let H = FOUR.pow(k) - 1;
    let Hs = H as f64 * s;
    let mut table:HashMap<i32, HashSet<i32>> = HashMap::new();

    for i in 0..seq_len + 1 -k as i32{
        let a = i as usize;
        let b = a+k as usize;
        let sub = &seq[a..b];
        println!("{sub}");

        let score = dna_encode(sub, k as i32);
        let mut pos = HashSet::new();
        pos.insert(i+add_on);

        let val = table.entry(score).or_insert(HashSet::new());
        val.insert(i+add_on);
    }
    table
}
