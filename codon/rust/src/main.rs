// Exact Motif-Free Codon Design - Rust Port
// Implements Aho-Corasick automaton, Exact DP, GC bounds, and validation.

use std::collections::{HashMap, VecDeque};
use std::env;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

// =========================================================================
// 1. Genetic Code Data Model
// =========================================================================
fn codon_gc(codon: &str) -> i32 {
    codon.chars().filter(|&c| c == 'G' || c == 'C').count() as i32
}

fn amino_acid_codons(aa: char) -> &'static [&'static str] {
    match aa {
        'A' => &["GCA", "GCC", "GCG", "GCT"],
        'C' => &["TGC", "TGT"],
        'D' => &["GAC", "GAT"],
        'E' => &["GAA", "GAG"],
        'F' => &["TTC", "TTT"],
        'G' => &["GGA", "GGC", "GGG", "GGT"],
        'H' => &["CAC", "CAT"],
        'I' => &["ATA", "ATC", "ATT"],
        'K' => &["AAA", "AAG"],
        'L' => &["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
        'M' => &["ATG"],
        'N' => &["AAC", "AAT"],
        'P' => &["CCA", "CCC", "CCG", "CCT"],
        'Q' => &["CAA", "CAG"],
        'R' => &["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
        'S' => &["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
        'T' => &["ACA", "ACC", "ACG", "ACT"],
        'V' => &["GTA", "GTC", "GTG", "GTT"],
        'W' => &["TGG"],
        'Y' => &["TAC", "TAT"],
        _ => &[],
    }
}

// =========================================================================
// 2. Aho-Corasick Automaton
// =========================================================================
#[derive(Clone, Default)]
struct AcNode {
    next: [usize; 4],
    fail: usize,
    out: bool,
}

fn base_index(b: char) -> usize {
    match b {
        'T' | 'U' => 0,
        'C' => 1,
        'A' => 2,
        'G' => 3,
        _ => 0,
    }
}

struct AhoCorasick {
    nodes: Vec<AcNode>,
}

impl AhoCorasick {
    fn new(patterns: &[&str]) -> Self {
        let mut nodes = vec![AcNode::default()];
        for &pat in patterns {
            let mut state = 0;
            for ch in pat.chars() {
                let b = base_index(ch);
                if nodes[state].next[b] == 0 {
                    nodes.push(AcNode::default());
                    nodes[state].next[b] = nodes.len() - 1;
                }
                state = nodes[state].next[b];
            }
            nodes[state].out = true;
        }

        // BFS failure links and completed transition table
        let mut queue = VecDeque::new();
        for b in 0..4 {
            let v = nodes[0].next[b];
            if v != 0 {
                nodes[v].fail = 0;
                queue.push_back(v);
            }
        }

        while let Some(u) = queue.pop_front() {
            for b in 0..4 {
                let v = nodes[u].next[b];
                if v != 0 {
                    let fail_target = nodes[nodes[u].fail].next[b];
                    nodes[v].fail = fail_target;
                    nodes[v].out = nodes[v].out || nodes[fail_target].out;
                    queue.push_back(v);
                } else {
                    nodes[u].next[b] = nodes[nodes[u].fail].next[b];
                }
            }
        }

        Self { nodes }
    }

    fn step(&self, state: usize, base: char) -> (usize, bool) {
        let b = base_index(base);
        let nxt = self.nodes[state].next[b];
        (nxt, self.nodes[nxt].out)
    }

    fn search(&self, text: &str) -> bool {
        let mut state = 0;
        for ch in text.chars() {
            let (nxt, out) = self.step(state, ch);
            state = nxt;
            if out {
                return true;
            }
        }
        false
    }
}

// =========================================================================
// 3. Exact Dynamic Programming
// =========================================================================
struct DpResult {
    feasible: bool,
    score: i64,
    total_states: usize,
}

fn solve_dp(protein: &str, motifs: &[&str], lo: i32, hi: i32) -> DpResult {
    let ac = AhoCorasick::new(motifs);
    // Key: (automaton state, GC count)
    // Value: best score
    let mut cur: HashMap<(usize, i32), i64> = HashMap::new();
    cur.insert((0, 0), 0);

    let mut total_states = 1;

    for aa in protein.chars() {
        let mut nxt: HashMap<(usize, i32), i64> = HashMap::new();
        let allowed = amino_acid_codons(aa);

        for (&(state, gc), &score) in &cur {
            for &codon in allowed {
                let mut qq = state;
                let mut bad = false;
                for ch in codon.chars() {
                    let (ns, out) = ac.step(qq, ch);
                    qq = ns;
                    if out {
                        bad = true;
                        break;
                    }
                }
                if bad {
                    continue;
                }
                let new_gc = gc + codon_gc(codon);
                if new_gc > hi {
                    continue;
                }

                let codon_score = 100 + codon_gc(codon) as i64 * 10;
                let new_score = score + codon_score;
                let entry = nxt.entry((qq, new_gc)).or_insert(i64::MIN);
                if new_score > *entry {
                    *entry = new_score;
                }
            }
        }
        total_states += nxt.len();
        cur = nxt;
        if cur.is_empty() {
            break;
        }
    }

    let mut best_score = i64::MIN;
    for (&(_, gc), &score) in &cur {
        if gc >= lo && gc <= hi && score > best_score {
            best_score = score;
        }
    }

    if best_score > i64::MIN {
        DpResult { feasible: true, score: best_score, total_states }
    } else {
        DpResult { feasible: false, score: 0, total_states }
    }
}

fn run_selftest() -> bool {
    let ac = AhoCorasick::new(&["GAATTC", "GGATCC"]);
    let found = ac.search("ACTGAATTCAAG");
    let not_found = ac.search("ACTGAGTTCAAG");
    let mut ok = found && !not_found;
    if ok {
        println!("PASS: Aho-Corasick motif search");
    } else {
        println!("FAIL: Aho-Corasick motif search");
    }

    let res = solve_dp("MKG", &["GAATTC"], 2, 7);
    if res.feasible {
        println!("PASS: Codon exact DP small instance (score={})", res.score);
    } else {
        println!("FAIL: Codon exact DP small instance infeasible");
        ok = false;
    }
    ok
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() > 1 && args[1] == "--test" {
        let ok = run_selftest();
        std::process::exit(if ok { 0 } else { 1 });
    }

    let out_csv = if args.len() > 2 && args[1] == "--bench" {
        &args[2]
    } else {
        "../results/codon_rust_results.csv"
    };

    println!("Codon Design Rust Port running...");
    run_selftest();

    let mut file = match File::create(out_csv) {
        Ok(f) => f,
        Err(e) => {
            println!("Could not open {}: {}", out_csv, e);
            return;
        }
    };

    writeln!(file, "language,experiment,method,n,motif_set,gc_mode,seed,rep,metric,value,unit").unwrap();

    let motifs = ["GAATTC", "GGATCC", "AAGCTT"];
    for &n in &[10, 20, 50, 100] {
        let protein: String = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"
            .chars()
            .cycle()
            .take(n)
            .collect();

        let start = Instant::now();
        let res = solve_dp(&protein, &motifs, 0, (3 * n) as i32);
        let dur_us = start.elapsed().as_micros();

        writeln!(file, "rust,E2,exact_dp,{},motifs_3,standard,42,0,runtime_us,{},us", n, dur_us).unwrap();
        writeln!(file, "rust,E2,exact_dp,{},motifs_3,standard,42,0,dp_states,{},count", n, res.total_states).unwrap();
        writeln!(file, "rust,E2,exact_dp,{},motifs_3,standard,42,0,feasible,{},bool", n, if res.feasible { 1 } else { 0 }).unwrap();
    }

    println!("Rust benchmark finished. Output saved to: {}", out_csv);
}

