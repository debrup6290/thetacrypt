//! Threshold ML-DSA parameters per PQShield (CCS '25 / USENIX Security '26).
//!
//! This module defines the parameters specific to the threshold variant of ML-DSA,
//! drawn from the paper's Section 3.5 (Figures 8, 9a, 9b) and Algorithm 6 for
//! the hardcoded RSSRecover partitions.
//!
//! ## Parameter Sources
//!
//! - **Figure 8**: ML-DSA-44 parameters (ν = 3) for all (T, N) with 2 ≤ T ≤ N ≤ 6
//! - **Figure 9a**: ML-DSA-65 parameters (ν = 6) for 2 ≤ T ≤ N ≤ 6
//! - **Figure 9b**: ML-DSA-87 parameters (ν = 7) for 2 ≤ T ≤ N ≤ 6
//! - **Algorithm 6**: Optimal RSS partitions for each (T, N) combination
//!
//! ## How Parameters Are Used
//!
//! | Parameter | Used in        | Purpose                                         |
//! |-----------|----------------|------------------------------------------------|
//! | N, T      | rss, protocol  | Party count and signing threshold               |
//! | K         | protocol       | Run K parallel instances to boost success prob   |
//! | ν (nu)    | hyperball      | Expansion factor: scale z^(1) coords by ν        |
//! | r         | hyperball      | Target hyperball radius (HRej acceptance check)  |
//! | r'        | hyperball      | Randomness hyperball radius (masking vector)     |
//! | M         | hyperball      | Per-party rejection bound: M = (r'/r)^{n(k+l)}  |
//!
//! ## Key Design Decisions
//!
//! - **ν is constant per security level**, not per (T, N). This simplifies
//!   implementation and matches the paper exactly.
//! - **r and r' are stored as integers** (the paper tabulates them this way).
//!   The HRej acceptance check is ||z||₂ ≤ r, i.e. ||z||₂² ≤ r².
//! - **M is derived**, not stored: M = (r'/r)^{n(k+l)}.

use crate::params::{Params, N as POLY_N};

//  Threshold Parameter Set 

/// Complete parameter set for one threshold ML-DSA configuration.
///
/// Each configuration is uniquely identified by (N, T, security_level).
#[derive(Debug, Clone, Copy)]
pub struct ThresholdParams {
    /// Reference to the base ML-DSA parameter set.
    pub base: &'static Params,

    // ── Party configuration ──

    /// Total number of parties (2 ≤ N ≤ 6).
    pub n_parties: usize,

    /// Signing threshold: any T parties can cooperate to sign (2 ≤ T ≤ N).
    pub threshold: usize,

    // ── Parallel instances ──

    /// Number of parallel signing instances (K).
    ///
    /// The protocol runs K independent signing attempts in parallel.
    /// K = ceil(1/p) where p is the per-instance success probability,
    /// targeting overall success probability ≥ 1/2.
    pub num_instances: usize,

    // ── Hyperball rejection sampling ──

    /// Expansion factor ν (nu) for imbalanced rejection sampling.
    ///
    /// In HRej (paper Figure 4), the first l coordinates of the candidate
    /// vector are divided by ν before the norm check, then multiplied by ν
    /// after acceptance.  This creates an asymmetric acceptance region.
    ///
    /// Constant per security level:
    ///   ML-DSA-44: ν = 3
    ///   ML-DSA-65: ν = 6
    ///   ML-DSA-87: ν = 7
    pub nu: u32,

    /// Target hyperball radius r.
    ///
    /// A candidate vector z is accepted by HRej only if ||z||₂ ≤ r.
    /// Stored as an integer (paper tabulates these directly).
    pub r: u64,

    /// Randomness hyperball radius r'.
    ///
    /// The masking randomness is sampled uniformly from the hyperball of
    /// radius r' in dimension n(k+l).  Must satisfy r' > r.
    pub r_prime: u64,

    // ── Derived / cached ──

    /// Number of RSS shares: C(N, N-T+1).
    pub num_shares: usize,
}

impl ThresholdParams {
    /// Dimension of the hyperball: n · (k + l), where n = 256.
    pub fn hyperball_dim(&self) -> usize {
        POLY_N * (self.base.k + self.base.l)
    }

    /// Size of each RSS subset: N - T + 1.
    pub fn subset_size(&self) -> usize {
        self.n_parties - self.threshold + 1
    }

    /// Number of shares held by each individual party: C(N-1, N-T).
    pub fn shares_per_party(&self) -> usize {
        binomial(self.n_parties - 1, self.n_parties - self.threshold)
    }

    /// Maximum number of shares any single party signs with in
    /// RSSRecover: ceil(C(N, N-T+1) / T).
    ///
    /// This directly determines the partial secret norm bound B
    /// in the paper's parameter selection (Section 3.4).
    pub fn max_shares_per_signer(&self) -> usize {
        (self.num_shares + self.threshold - 1) / self.threshold
    }
}

//  Combinatorial Helpers 

/// Binomial coefficient C(n, k) = n! / (k! · (n-k)!).
///
/// For the small values in threshold ML-DSA (n ≤ 6), this never overflows.
/// Returns 0 if k > n.
pub fn binomial(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }
    let k = k.min(n - k);
    let mut result = 1usize;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

/// Enumerate all subsets of {0, 1, ..., n-1} with exactly `size` elements.
///
/// Returns subsets in lexicographic order.  Each subset is a sorted `Vec<usize>`.
pub fn enumerate_subsets(n: usize, size: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::with_capacity(binomial(n, size));
    let mut current = Vec::with_capacity(size);
    enumerate_subsets_rec(n, size, 0, &mut current, &mut result);
    result
}

fn enumerate_subsets_rec(
    n: usize,
    size: usize,
    start: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == size {
        result.push(current.clone());
        return;
    }
    let remaining_needed = size - current.len();
    let remaining_available = n - start;
    if remaining_available < remaining_needed {
        return;
    }
    for i in start..n {
        current.push(i);
        enumerate_subsets_rec(n, size, i + 1, current, result);
        current.pop();
    }
}

/// Return the indices of all shares held by a given party.
///
/// Party `party_id` holds share `s_I` for every subset I containing `party_id`.
pub fn shares_for_party(n_parties: usize, threshold: usize, party_id: usize) -> Vec<usize> {
    let subset_size = n_parties - threshold + 1;
    let subsets = enumerate_subsets(n_parties, subset_size);
    subsets.iter()
        .enumerate()
        .filter(|(_, subset)| subset.contains(&party_id))
        .map(|(idx, _)| idx)
        .collect()
}

//  RSSRecover 
//
// Paper Algorithm 6: hardcoded optimal partitions for each (T, N).
//
// The partition assigns each share (subset) to exactly one active signer,
// minimizing the maximum number of shares any single party handles.
//
// The paper gives partitions for a canonical active set {0, 1, ..., T-1}.
// For an arbitrary active set, we permute indices accordingly.

/// Compute the optimal RSS partition for a given active set.
///
/// Returns a `Vec` of length `n_parties`, where entry `i` (for i ∈ act) is
/// the list of share subset indices assigned to party i, and entry `i`
/// (for i ∉ act) is empty.
///
/// The share subsets are indexed in the same order as `enumerate_subsets(n, n-t+1)`.
pub fn rss_recover(n_parties: usize, threshold: usize, act: &[usize]) -> Vec<Vec<usize>> {
    assert_eq!(act.len(), threshold, "Active set size must equal threshold");

    let subset_size = n_parties - threshold + 1;
    let all_subsets = enumerate_subsets(n_parties, subset_size);

    // T = N: trivial case, each party holds exactly one share (its singleton)
    if threshold == n_parties {
        let mut partition = vec![vec![]; n_parties];
        for (idx, subset) in all_subsets.iter().enumerate() {
            let party = subset[0];
            if act.contains(&party) {
                partition[party].push(idx);
            }
        }
        return partition;
    }

    // Get the canonical partition (for act = {0, 1, ..., T-1}) from Algorithm 6.
    let canonical = canonical_partition(n_parties, threshold);

    // Build permutation φ: map canonical indices → actual party indices.
    let mut active_sorted: Vec<usize> = act.to_vec();
    active_sorted.sort();
    let inactive_sorted: Vec<usize> = (0..n_parties)
        .filter(|x| !active_sorted.contains(x))
        .collect();

    let mut phi = vec![0usize; n_parties];
    for (i, &a) in active_sorted.iter().enumerate() {
        phi[i] = a;
    }
    for (i, &a) in inactive_sorted.iter().enumerate() {
        phi[threshold + i] = a;
    }

    // Translate canonical partition through the permutation.
    let mut partition = vec![vec![]; n_parties];
    for (canonical_idx, canonical_subsets) in canonical.iter().enumerate() {
        let actual_party = phi[canonical_idx];
        for canonical_subset in canonical_subsets {
            let mut actual_subset: Vec<usize> = canonical_subset.iter()
                .map(|&x| phi[x])
                .collect();
            actual_subset.sort();

            if let Some(share_idx) = all_subsets.iter().position(|s| *s == actual_subset) {
                partition[actual_party].push(share_idx);
            }
        }
    }

    partition
}

/// Canonical partitions from Paper Algorithm 6.
///
/// Returns partitions for the canonical active set {0, 1, ..., T-1},
/// where each inner Vec contains the subset tuples assigned to that signer.
/// Inactive parties (indices T..N) receive empty assignments.
fn canonical_partition(n: usize, t: usize) -> Vec<Vec<Vec<usize>>> {
    let mut result = vec![vec![]; n];

    match (n, t) {
        (3, 2) => {
            result[0] = vec![vec![0, 1], vec![0, 2]];
            result[1] = vec![vec![1, 2]];
        }
        (4, 2) => {
            result[0] = vec![vec![0, 1, 3], vec![0, 2, 3]];
            result[1] = vec![vec![0, 1, 2], vec![1, 2, 3]];
        }
        (4, 3) => {
            result[0] = vec![vec![0, 1], vec![0, 3]];
            result[1] = vec![vec![1, 2], vec![1, 3]];
            result[2] = vec![vec![2, 3], vec![0, 2]];
        }
        (5, 2) => {
            result[0] = vec![vec![0, 1, 3, 4], vec![0, 2, 3, 4], vec![0, 1, 2, 4]];
            result[1] = vec![vec![1, 2, 3, 4], vec![0, 1, 2, 3]];
        }
        (5, 3) => {
            result[0] = vec![vec![0, 3, 4], vec![0, 1, 3], vec![0, 1, 4], vec![0, 2, 3]];
            result[1] = vec![vec![0, 1, 2], vec![1, 2, 3], vec![1, 2, 4], vec![1, 3, 4]];
            result[2] = vec![vec![2, 3, 4], vec![0, 2, 4]];
        }
        (5, 4) => {
            result[0] = vec![vec![0, 1], vec![0, 3], vec![0, 4]];
            result[1] = vec![vec![1, 2], vec![1, 3], vec![1, 4]];
            result[2] = vec![vec![2, 3], vec![0, 2], vec![2, 4]];
            result[3] = vec![vec![3, 4]];
        }
        (6, 2) => {
            result[0] = vec![vec![0, 2, 3, 4, 5], vec![0, 1, 2, 3, 5], vec![0, 1, 2, 4, 5]];
            result[1] = vec![vec![1, 2, 3, 4, 5], vec![0, 1, 2, 3, 4], vec![0, 1, 3, 4, 5]];
        }
        (6, 3) => {
            result[0] = vec![
                vec![0, 1, 3, 4], vec![0, 1, 2, 4], vec![0, 1, 3, 5],
                vec![0, 3, 4, 5], vec![0, 1, 2, 5],
            ];
            result[1] = vec![
                vec![0, 1, 4, 5], vec![1, 3, 4, 5], vec![1, 2, 3, 5],
                vec![1, 2, 3, 4], vec![1, 2, 4, 5],
            ];
            result[2] = vec![
                vec![0, 2, 3, 5], vec![0, 2, 4, 5], vec![0, 2, 3, 4],
                vec![0, 1, 2, 3], vec![2, 3, 4, 5],
            ];
        }
        (6, 4) => {
            result[0] = vec![
                vec![0, 1, 4], vec![0, 2, 3], vec![0, 1, 5],
                vec![0, 1, 2], vec![0, 4, 5],
            ];
            result[1] = vec![
                vec![1, 3, 5], vec![1, 3, 4], vec![1, 2, 5],
                vec![1, 4, 5], vec![1, 2, 4],
            ];
            result[2] = vec![
                vec![2, 4, 5], vec![0, 2, 4], vec![2, 3, 5],
                vec![2, 3, 4], vec![0, 2, 5],
            ];
            result[3] = vec![
                vec![0, 3, 4], vec![0, 1, 3], vec![1, 2, 3],
                vec![3, 4, 5], vec![0, 3, 5],
            ];
        }
        (6, 5) => {
            result[0] = vec![vec![0, 1], vec![0, 2], vec![0, 5]];
            result[1] = vec![vec![1, 2], vec![1, 3], vec![1, 5]];
            result[2] = vec![vec![2, 3], vec![2, 4], vec![2, 5]];
            result[3] = vec![vec![0, 3], vec![3, 4], vec![3, 5]];
            result[4] = vec![vec![4, 5], vec![0, 4], vec![1, 4]];
        }
        _ => panic!("Unsupported (N, T) = ({}, {}) for canonical partition", n, t),
    }

    result
}

//  Parameter Tables 
//
// Transcribed from the paper's Section 3.5:
//   - Figure 8:  ML-DSA-44 (ν = 3)
//   - Figure 9a: ML-DSA-65 (ν = 6)
//   - Figure 9b: ML-DSA-87 (ν = 7)

/// Look up threshold parameters for a given configuration.
///
/// Returns `None` if the (N, T, base) combination is not supported.
pub fn lookup(n_parties: usize, threshold: usize, base: &'static Params) -> Option<ThresholdParams> {
    if n_parties < 2 || n_parties > 6 || threshold < 2 || threshold > n_parties {
        return None;
    }

    let num_shares = binomial(n_parties, n_parties - threshold + 1);
    let (num_instances, r, r_prime, nu) = lookup_table_entry(n_parties, threshold, base)?;

    Some(ThresholdParams {
        base,
        n_parties,
        threshold,
        num_instances,
        nu,
        r,
        r_prime,
        num_shares,
    })
}

/// Internal: look up (K, r, r', ν) from the paper's parameter tables.
fn lookup_table_entry(
    n_parties: usize,
    threshold: usize,
    base: &Params,
) -> Option<(usize, u64, u64, u32)> {
    //         (K,   r,       r',      ν)
    match (base.name, n_parties, threshold) {
        //  ML-DSA-44 (Figure 8, ν = 3) 
        ("ML-DSA-44", 2, 2) => Some((  2, 252778, 252833, 3)),
        ("ML-DSA-44", 3, 2) => Some((  3, 310060, 310138, 3)),
        ("ML-DSA-44", 3, 3) => Some((  4, 246490, 246546, 3)),
        ("ML-DSA-44", 4, 2) => Some((  3, 310060, 310138, 3)),
        ("ML-DSA-44", 4, 3) => Some((  7, 279235, 279314, 3)),
        ("ML-DSA-44", 4, 4) => Some((  8, 243463, 243519, 3)),
        ("ML-DSA-44", 5, 2) => Some((  3, 285363, 285459, 3)),
        ("ML-DSA-44", 5, 3) => Some(( 14, 282800, 282912, 3)),
        ("ML-DSA-44", 5, 4) => Some(( 30, 259427, 259526, 3)),
        ("ML-DSA-44", 5, 5) => Some(( 16, 239924, 239981, 3)),
        ("ML-DSA-44", 6, 2) => Some((  4, 300265, 300362, 3)),
        ("ML-DSA-44", 6, 3) => Some(( 19, 277014, 277139, 3)),
        ("ML-DSA-44", 6, 4) => Some(( 74, 268705, 268831, 3)),
        ("ML-DSA-44", 6, 5) => Some((100, 250590, 250686, 3)),
        ("ML-DSA-44", 6, 6) => Some(( 37, 219245, 219301, 3)),

        //  ML-DSA-65 (Figure 9a, ν = 6) 
        ("ML-DSA-65", 2, 2) => Some((   3, 501495, 501613, 6)),
        ("ML-DSA-65", 3, 2) => Some((   5, 540212, 540378, 6)),
        ("ML-DSA-65", 3, 3) => Some((   9, 510387, 510504, 6)),
        ("ML-DSA-65", 4, 2) => Some((   6, 540212, 540378, 6)),
        ("ML-DSA-65", 4, 3) => Some((  20, 506761, 506928, 6)),
        ("ML-DSA-65", 4, 4) => Some((  26, 433594, 433711, 6)),
        ("ML-DSA-65", 5, 2) => Some((   8, 552371, 552575, 6)),
        ("ML-DSA-65", 5, 3) => Some((  62, 552909, 553145, 6)),
        ("ML-DSA-65", 5, 4) => Some(( 205, 474331, 474535, 6)),
        ("ML-DSA-65", 5, 5) => Some((  78, 425914, 426032, 6)),
        ("ML-DSA-65", 6, 2) => Some((   8, 571208, 571412, 6)),
        ("ML-DSA-65", 6, 3) => Some((  95, 536793, 537058, 6)),
        ("ML-DSA-65", 6, 4) => Some(( 804, 488704, 488969, 6)),
        ("ML-DSA-65", 6, 5) => Some((1200, 461324, 461529, 6)),
        ("ML-DSA-65", 6, 6) => Some(( 250, 414896, 415013, 6)),

        //  ML-DSA-87 (Figure 9b, ν = 7) 
        ("ML-DSA-87", 2, 2) => Some((   3, 503119, 503192, 7)),
        ("ML-DSA-87", 3, 2) => Some((   4, 631601, 631703, 7)),
        ("ML-DSA-87", 3, 3) => Some((   6, 483107, 483180, 7)),
        ("ML-DSA-87", 4, 2) => Some((   4, 632903, 633006, 7)),
        ("ML-DSA-87", 4, 3) => Some((  11, 551752, 551854, 7)),
        ("ML-DSA-87", 4, 4) => Some((  14, 487958, 488031, 7)),
        ("ML-DSA-87", 5, 2) => Some((   5, 607694, 607820, 7)),
        ("ML-DSA-87", 5, 3) => Some((  26, 577400, 577546, 7)),
        ("ML-DSA-87", 5, 4) => Some((  70, 518384, 518510, 7)),
        ("ML-DSA-87", 5, 5) => Some((  35, 468214, 468287, 7)),
        ("ML-DSA-87", 6, 2) => Some((   5, 665106, 665232, 7)),
        ("ML-DSA-87", 6, 3) => Some((  39, 577541, 577704, 7)),
        ("ML-DSA-87", 6, 4) => Some(( 208, 517689, 517853, 7)),
        ("ML-DSA-87", 6, 5) => Some(( 295, 479692, 479819, 7)),
        ("ML-DSA-87", 6, 6) => Some((  87, 424124, 424197, 7)),

        _ => None,
    }
}

//  Tests 

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87};

    // ── Binomial coefficient ──

    #[test]
    fn test_binomial_base_cases() {
        assert_eq!(binomial(0, 0), 1);
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 5), 1);
        assert_eq!(binomial(3, 4), 0);
    }

    #[test]
    fn test_binomial_known_values() {
        assert_eq!(binomial(4, 2), 6);
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(6, 3), 20);
        assert_eq!(binomial(6, 4), 15);
    }

    #[test]
    fn test_binomial_symmetry() {
        for n in 0..=6 {
            for k in 0..=n {
                assert_eq!(binomial(n, k), binomial(n, n - k));
            }
        }
    }

    #[test]
    fn test_binomial_pascals_triangle() {
        for n in 1..=6 {
            for k in 1..n {
                assert_eq!(binomial(n, k), binomial(n - 1, k - 1) + binomial(n - 1, k));
            }
        }
    }

    // ── Subset enumeration ──

    #[test]
    fn test_enumerate_subsets_count() {
        for n in 1..=6 {
            for k in 1..=n {
                assert_eq!(enumerate_subsets(n, k).len(), binomial(n, k));
            }
        }
    }

    #[test]
    fn test_enumerate_subsets_3_choose_2() {
        assert_eq!(enumerate_subsets(3, 2), vec![vec![0, 1], vec![0, 2], vec![1, 2]]);
    }

    // ── Share assignment ──

    #[test]
    fn test_shares_for_party_3_2() {
        assert_eq!(shares_for_party(3, 2, 0), vec![0, 1]);
        assert_eq!(shares_for_party(3, 2, 1), vec![0, 2]);
        assert_eq!(shares_for_party(3, 2, 2), vec![1, 2]);
    }

    #[test]
    fn test_shares_cover_all_subsets() {
        for n in 2..=5 {
            for t in 2..=n {
                let num_shares = binomial(n, n - t + 1);
                let active: Vec<usize> = (0..t).collect();
                let mut covered = vec![false; num_shares];
                for &party in &active {
                    for idx in shares_for_party(n, t, party) {
                        covered[idx] = true;
                    }
                }
                assert!(covered.iter().all(|&c| c));
            }
        }
    }

    // ── ThresholdParams lookup ──

    #[test]
    fn test_lookup_basic_44() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        assert_eq!(tp.n_parties, 3);
        assert_eq!(tp.threshold, 2);
        assert_eq!(tp.num_shares, 3);
        assert_eq!(tp.num_instances, 3);
        assert_eq!(tp.nu, 3);
        assert_eq!(tp.r, 310060);
        assert_eq!(tp.r_prime, 310138);
    }

    #[test]
    fn test_lookup_rejects_invalid() {
        assert!(lookup(1, 1, &ML_DSA_44).is_none());
        assert!(lookup(7, 2, &ML_DSA_44).is_none());
        assert!(lookup(3, 4, &ML_DSA_44).is_none());
        assert!(lookup(3, 1, &ML_DSA_44).is_none());
    }

    #[test]
    fn test_lookup_all_configs() {
        for base in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            for n in 2..=6 {
                for t in 2..=n {
                    assert!(lookup(n, t, base).is_some(),
                        "{} ({},{}) should be supported", base.name, n, t);
                }
            }
        }
    }

    #[test]
    fn test_nu_constant_per_level() {
        for n in 2..=6 {
            for t in 2..=n {
                assert_eq!(lookup(n, t, &ML_DSA_44).unwrap().nu, 3);
                assert_eq!(lookup(n, t, &ML_DSA_65).unwrap().nu, 6);
                assert_eq!(lookup(n, t, &ML_DSA_87).unwrap().nu, 7);
            }
        }
    }

    #[test]
    fn test_r_prime_greater_than_r() {
        for base in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            for n in 2..=6 {
                for t in 2..=n {
                    let tp = lookup(n, t, base).unwrap();
                    assert!(tp.r_prime > tp.r,
                        "{} ({},{}): r' must exceed r", base.name, n, t);
                }
            }
        }
    }

    // ── Derived values ──

    #[test]
    fn test_threshold_params_derived() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        assert_eq!(tp.subset_size(), 2);
        assert_eq!(tp.shares_per_party(), 2);
        assert_eq!(tp.hyperball_dim(), 256 * 8);
        assert_eq!(tp.max_shares_per_signer(), 2);
    }

    // ── RSSRecover ──

    #[test]
    fn test_rss_recover_3_2_canonical() {
        let partition = rss_recover(3, 2, &[0, 1]);
        assert_eq!(partition[0], vec![0, 1]); // {01}, {02}
        assert_eq!(partition[1], vec![2]);     // {12}
        assert!(partition[2].is_empty());
    }

    #[test]
    fn test_rss_recover_covers_all_shares() {
        for n in 2..=6 {
            for t in 2..=n {
                let act: Vec<usize> = (0..t).collect();
                let partition = rss_recover(n, t, &act);
                let num_shares = binomial(n, n - t + 1);
                let mut assigned = vec![false; num_shares];
                for i in 0..n {
                    for &idx in &partition[i] {
                        assert!(!assigned[idx], "({},{}) share {} assigned twice", n, t, idx);
                        assigned[idx] = true;
                    }
                }
                assert!(assigned.iter().all(|&a| a), "({},{}) not all shares assigned", n, t);
            }
        }
    }

    #[test]
    fn test_rss_recover_party_owns_shares() {
        for n in 3..=6 {
            for t in 2..=n {
                let act: Vec<usize> = (0..t).collect();
                let partition = rss_recover(n, t, &act);
                let all_subsets = enumerate_subsets(n, n - t + 1);
                for party in 0..n {
                    for &share_idx in &partition[party] {
                        assert!(all_subsets[share_idx].contains(&party),
                            "({},{}) party {} assigned share {:?} it doesn't own",
                            n, t, party, all_subsets[share_idx]);
                    }
                }
            }
        }
    }

    #[test]
    fn test_rss_recover_balanced() {
        for n in 2..=6 {
            for t in 2..=n {
                let act: Vec<usize> = (0..t).collect();
                let partition = rss_recover(n, t, &act);
                let num_shares = binomial(n, n - t + 1);
                let expected_max = (num_shares + t - 1) / t;
                let actual_max = partition.iter().map(|p| p.len()).max().unwrap();
                assert!(actual_max <= expected_max,
                    "({},{}) max {} > expected {}", n, t, actual_max, expected_max);
            }
        }
    }

    #[test]
    fn test_rss_recover_permuted_active_set() {
        let partition = rss_recover(4, 3, &[1, 2, 3]);
        let num_shares = binomial(4, 2);
        let mut assigned = vec![false; num_shares];
        for i in 0..4 {
            for &idx in &partition[i] {
                assigned[idx] = true;
            }
        }
        assert!(assigned.iter().all(|&a| a));
        assert!(partition[0].is_empty());
    }
}