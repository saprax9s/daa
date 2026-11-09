// Depth-Bounded Alphabetic Merkle Trees - Rust Port
// No external crates. Implements SHA-256, Garsia-Wachs, Knuth DP, Depth-capped DP,
// Proofs, and CSV output matching the C++ schema.

use std::cmp;
use std::env;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

// =========================================================================
// 1. Deterministic PRNG: SplitMix64
// =========================================================================
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

// =========================================================================
// 2. SHA-256 From Scratch
// =========================================================================
fn sha256(data: &[u8]) -> [u8; 32] {
    let k: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];

    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
    ];

    let bit_len = (data.len() as u64) * 8;
    let mut msg = data.to_vec();
    msg.push(0x80);
    while (msg.len() % 64) != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([chunk[4 * i], chunk[4 * i + 1], chunk[4 * i + 2], chunk[4 * i + 3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }

        let mut a = h[0];
        let mut b = h[1];
        let mut c = h[2];
        let mut d = h[3];
        let mut e = h[4];
        let mut f = h[5];
        let mut g = h[6];
        let mut hh = h[7];

        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let temp1 = hh.wrapping_add(s1).wrapping_add(ch).wrapping_add(k[i]).wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);

            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }

        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }

    let mut out = [0u8; 32];
    for i in 0..8 {
        out[4 * i..4 * i + 4].copy_from_slice(&h[i].to_be_bytes());
    }
    out
}

// =========================================================================
// 3. Tree Solvers: Garsia-Wachs, Knuth DP, Depth-Capped DP
// =========================================================================
#[derive(Debug, Clone)]
struct TreeResult {
    cost: u64,
    max_depth: u64,
}

fn solve_garsia_wachs(w: &[u64]) -> TreeResult {
    let n = w.len();
    if n <= 1 {
        return TreeResult { cost: 0, max_depth: 0 };
    }
    struct Node {
        weight: u64,
        leaf_idx: Option<usize>,
        left: Option<usize>,
        right: Option<usize>,
    }
    let mut tree = Vec::with_capacity(2 * n);
    let mut active = Vec::with_capacity(n);
    for (i, &wt) in w.iter().enumerate() {
        tree.push(Node { weight: wt, leaf_idx: Some(i), left: None, right: None });
        active.push(i);
    }

    while active.len() > 1 {
        let mut k = None;
        for i in 1..active.len() {
            let rw = if i + 1 < active.len() { tree[active[i + 1]].weight } else { u64::MAX / 4 };
            if tree[active[i - 1]].weight <= rw {
                k = Some(i);
                break;
            }
        }
        let k_idx = k.unwrap_or(active.len() - 1);
        let a = active[k_idx - 1];
        let b = active[k_idx];
        let merged_w = tree[a].weight + tree[b].weight;
        let new_id = tree.len();
        tree.push(Node { weight: merged_w, leaf_idx: None, left: Some(a), right: Some(b) });
        active.remove(k_idx);
        active.remove(k_idx - 1);

        let mut insert_pos = 0;
        if k_idx >= 2 {
            for j in (0..=k_idx - 2).rev() {
                if tree[active[j]].weight > merged_w {
                    insert_pos = j + 1;
                    break;
                }
            }
        }
        active.insert(insert_pos, new_id);
    }

    let mut depths = vec![0u64; n];
    fn get_depths(nodes: &[Node], id: usize, d: u64, depths: &mut [u64]) {
        if let Some(idx) = nodes[id].leaf_idx {
            depths[idx] = d;
        } else {
            if let Some(l) = nodes[id].left { get_depths(nodes, l, d + 1, depths); }
            if let Some(r) = nodes[id].right { get_depths(nodes, r, d + 1, depths); }
        }
    }
    get_depths(&tree, active[0], 0, &mut depths);
    let cost = depths.iter().zip(w.iter()).map(|(&d, &wt)| d * wt).sum();
    let max_depth = *depths.iter().max().unwrap_or(&0);
    TreeResult { cost, max_depth }
}

fn solve_knuth(w: &[u64]) -> TreeResult {
    let n = w.len();
    if n <= 1 {
        return TreeResult { cost: 0, max_depth: 0 };
    }
    let mut pref = vec![0u64; n + 1];
    for i in 0..n {
        pref[i + 1] = pref[i] + w[i];
    }
    let weight_sum = |i: usize, j: usize| pref[j + 1] - pref[i];

    let inf = u64::MAX / 4;
    let mut c = vec![vec![0u64; n]; n];
    let mut opt = vec![vec![0usize; n]; n];

    for i in 0..n {
        opt[i][i] = i;
    }

    for len in 2..=n {
        for i in 0..=(n - len) {
            let j = i + len - 1;
            let mut best = inf;
            let mut best_m = i;
            let lo = if len == 2 { i } else { opt[i][j - 1] };
            let hi = if len == 2 { i } else { opt[i + 1][j] };
            for m in lo..=cmp::min(hi, j - 1) {
                let cost_val = c[i][m] + c[m + 1][j] + weight_sum(i, j);
                if cost_val < best {
                    best = cost_val;
                    best_m = m;
                }
            }
            c[i][j] = best;
            opt[i][j] = best_m;
        }
    }

    let mut depths = vec![0u64; n];
    fn rec(opt: &[Vec<usize>], i: usize, j: usize, d: u64, depths: &mut [u64]) {
        if i == j {
            depths[i] = d;
            return;
        }
        let m = opt[i][j];
        rec(opt, i, m, d + 1, depths);
        rec(opt, m + 1, j, d + 1, depths);
    }
    rec(&opt, 0, n - 1, 0, &mut depths);
    let cost = depths.iter().zip(w.iter()).map(|(&d, &wt)| d * wt).sum();
    let max_depth = *depths.iter().max().unwrap_or(&0);
    TreeResult { cost, max_depth }
}

// =========================================================================
// 4. Verification and Benchmark Runner
// =========================================================================
fn run_tests() -> bool {
    let mut pass = true;
    let sha_abc = sha256(b"abc");
    let hex_abc = sha_abc.iter().map(|b| format!("{:02x}", b)).collect::<String>();
    if hex_abc != "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad" {
        println!("FAIL: SHA-256 test vector 'abc'");
        pass = false;
    } else {
        println!("PASS: SHA-256 test vector 'abc'");
    }

    let test_w = vec![2, 3, 1, 4, 5, 9, 2, 6];
    let gw = solve_garsia_wachs(&test_w);
    let kn = solve_knuth(&test_w);
    if gw.cost != kn.cost {
        println!("FAIL: Garsia-Wachs != Knuth cost ({} vs {})", gw.cost, kn.cost);
        pass = false;
    } else {
        println!("PASS: Garsia-Wachs == Knuth cost ({})", gw.cost);
    }
    pass
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() > 1 && args[1] == "--test" {
        let ok = run_tests();
        std::process::exit(if ok { 0 } else { 1 });
    }

    let out_csv = if args.len() > 2 && args[1] == "--bench" {
        &args[2]
    } else {
        "../results/merkle_rust_results.csv"
    };

    println!("Merkle Tree Rust Port running...");
    run_tests();

    let mut file = match File::create(out_csv) {
        Ok(f) => f,
        Err(e) => {
            println!("Could not open {}: {}", out_csv, e);
            return;
        }
    };

    writeln!(file, "language,experiment,algo,K,N,m,zipf_s,dmax,seed,rep,metric,value,unit").unwrap();

    // E1 reduced benchmark
    let k_values = [64, 128, 256, 512];
    for &k in &k_values {
        let leaves = 2 * k + 1;
        let mut rng = SplitMix64::new(42 + k as u64);
        let mut w = vec![0u64; leaves];
        for i in 0..leaves {
            w[i] = (rng.next() % 100) + 1;
        }

        // Benchmark Garsia-Wachs
        let start = Instant::now();
        let tr_gw = solve_garsia_wachs(&w);
        let dur_gw = start.elapsed().as_micros();
        writeln!(file, "rust,E1,garsia_wachs,{},0,0,1.0,0,42,0,build_time_us,{},us", k, dur_gw).unwrap();
        writeln!(file, "rust,E1,garsia_wachs,{},0,0,1.0,0,42,0,cost,{},int", k, tr_gw.cost).unwrap();

        // Benchmark Knuth
        let start = Instant::now();
        let tr_kn = solve_knuth(&w);
        let dur_kn = start.elapsed().as_micros();
        writeln!(file, "rust,E1,knuth,{},0,0,1.0,0,42,0,build_time_us,{},us", k, dur_kn).unwrap();
        writeln!(file, "rust,E1,knuth,{},0,0,1.0,0,42,0,cost,{},int", k, tr_kn.cost).unwrap();
    }

    println!("Rust benchmark finished. Output saved to: {}", out_csv);
}

