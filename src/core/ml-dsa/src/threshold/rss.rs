//! Replicated Short Secret Sharing (RSS) for Threshold ML-DSA.
//!
//! Implements the key generation procedure from Paper Figure 5:
//!
//! 1. **RSS**: Sample C(N, N-T+1) independent short secret vectors s_I from χ_s,
//!    one for each subset I ⊆ [N] of size N-T+1. The full secret is s = Σ s_I.
//!    Each party i receives all s_I where i ∈ I.
//!
//! 2. **Keygen**: Expand ρ → A, compute t = [A|I]·s, apply Power2Round,
//!    produce the standard ML-DSA verification key vk = (ρ, t₁).
//!
//! 3. **Partial secret reconstruction**: Given an RSSRecover partition,
//!    each party sums its assigned shares to form s_part_i for signing.
//!
//! ## Security Property
//!
//! Since there are at most T-1 corrupted parties, at least one subset I of
//! size N-T+1 is fully non-corrupted. That subset's secret s_I remains
//! hidden from the adversary, so the public key retains pseudorandomness
//! under the same MLWE assumption as standard ML-DSA.
//!
//! ## Example (N=3, T=2)
//!
//! ```text
//! Subsets of size 2: I₀={0,1}, I₁={0,2}, I₂={1,2}
//!
//! Sample: s_{01}, s_{02}, s_{12}  (each with coefficients in [-η, η])
//! Full secret: s = s_{01} + s_{02} + s_{12}
//!
//! Party 0 holds: {s_{01}, s_{02}}     (subsets containing 0)
//! Party 1 holds: {s_{01}, s_{12}}     (subsets containing 1)
//! Party 2 holds: {s_{02}, s_{12}}     (subsets containing 2)
//!
//! If act = {0, 2} (parties 0 and 2 sign):
//!   RSSRecover → party 0 uses {s_{01}, s_{02}}, party 2 uses {s_{12}}
//!   s_part_0 = s_{01} + s_{02}
//!   s_part_2 = s_{12}
//!   s_part_0 + s_part_2 = s  ✓
//! ```

use crate::arithmetic;
use crate::encoding;
use crate::params::N as POLY_N;
use crate::poly::{self, Poly, PolyVec};
use crate::sampling;
use crate::shake;

use super::params::{self as threshold_params, ThresholdParams};

// ---

/// ---
///
/// ---
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RssShare {
    /// ---
    /// ---
    pub subset: Vec<usize>,

    /// ---
    pub s1: PolyVec,

    /// ---
    pub s2: PolyVec,
}

/// ---
///
/// ---
/// ---
/// ---
#[derive(serde::Serialize, serde::Deserialize)]
pub struct PartyKey {
    pub party_id: usize,
    #[serde(with = "serde_arrays")]
    pub tr: [u8; 64],
    pub shares: Vec<RssShare>,
}

/// ---
///
/// ---
/// ---
/// ---
#[derive(serde::Serialize, serde::Deserialize)]
pub struct VerificationKey {
    pub rho: [u8; 32],       // ---
    pub t1: PolyVec,
    pub packed: Vec<u8>,
    #[serde(with = "serde_arrays")]
    pub tr: [u8; 64],
}

/// ---
pub struct ThresholdKeyBundle {
    /// ---
    pub vk: VerificationKey,

    /// ---
    pub party_keys: Vec<PartyKey>,

    /// ---
    pub params: ThresholdParams,
}

// ---
// ---
//
// ---
// ---
// ---
// ---
// ---
// ---
// ---

/// ---
///
/// ---
/// ---
/// ---
///
/// ---
/// ---
pub fn keygen(xi: &[u8; 32], tp: &ThresholdParams) -> ThresholdKeyBundle {
    let base = tp.base;

    // ---
    //
    // ---
    // ---
    let total_bytes = 32 + tp.num_shares * 64;
    let expanded = shake::shake256(xi, total_bytes);

    let mut rho = [0u8; 32];
    rho.copy_from_slice(&expanded[0..32]);

    // ---
    let a_hat = sampling::expand_a(&rho, base);

    // ---
    let all_subsets = threshold_params::enumerate_subsets(tp.n_parties, tp.subset_size());
    debug_assert_eq!(all_subsets.len(), tp.num_shares);

    let mut all_shares: Vec<RssShare> = Vec::with_capacity(tp.num_shares);

    // ---
    let mut s1_full = poly::polyvec_zero(base.l);
    let mut s2_full = poly::polyvec_zero(base.k);

    for (share_idx, subset) in all_subsets.iter().enumerate() {
        // ---
        let seed_offset = 32 + share_idx * 64;
        let mut share_seed = [0u8; 64];
        share_seed.copy_from_slice(&expanded[seed_offset..seed_offset + 64]);

        // ---
        let (s1_i, s2_i) = sampling::expand_s(&share_seed, base);

        // ---
        s1_full = poly::polyvec_add(&s1_full, &s1_i);
        s2_full = poly::polyvec_add(&s2_full, &s2_i);

        all_shares.push(RssShare {
            subset: subset.clone(),
            s1: s1_i,
            s2: s2_i,
        });
    }

    // ---
    // ---
    let s1_hat = poly::polyvec_ntt(&s1_full);
    let t_hat = poly::mat_vec_mul_ntt(&a_hat, &s1_hat);
    let t_normal = poly::polyvec_inv_ntt(&t_hat);
    let t = poly::polyvec_add(&t_normal, &s2_full);

    // ---
    let mut t1 = Vec::with_capacity(base.k);
    for i in 0..base.k {
        let mut t1_poly = Poly::zero();
        for j in 0..POLY_N {
            let (hi, _lo) = arithmetic::power2round(t[i].coeffs[j]);
            t1_poly.coeffs[j] = hi;
        }
        t1.push(t1_poly);
    }

    // ---
    let pk_bytes = encoding::pk_encode(&rho, &t1, base);
    let tr_vec = shake::shake256(&pk_bytes, 64);
    let mut tr = [0u8; 64];
    tr.copy_from_slice(&tr_vec);

    // ---
    let mut party_keys = Vec::with_capacity(tp.n_parties);
    for party_id in 0..tp.n_parties {
        let party_share_indices = threshold_params::shares_for_party(
            tp.n_parties, tp.threshold, party_id,
        );
        let shares: Vec<RssShare> = party_share_indices.iter()
            .map(|&idx| all_shares[idx].clone())
            .collect();

        party_keys.push(PartyKey {
            party_id,
            tr,
            shares,
        });
    }

    ThresholdKeyBundle {
        vk: VerificationKey {
            rho,
            t1,
            packed: pk_bytes,
            tr,
        },
        party_keys,
        params: *tp,
    }
}

// ---
//
// ---
// ---
// ---

/// ---
///
/// ---
/// ---
/// ---
///
/// ---
pub fn partial_secret(
    party_key: &PartyKey,
    assigned_share_indices: &[usize],
    tp: &ThresholdParams,
) -> (PolyVec, PolyVec) {
    let base = tp.base;
    let all_subsets = threshold_params::enumerate_subsets(tp.n_parties, tp.subset_size());

    let mut s1_part = poly::polyvec_zero(base.l);
    let mut s2_part = poly::polyvec_zero(base.k);

    for &share_idx in assigned_share_indices {
        let target_subset = &all_subsets[share_idx];

        // ---
        let share = party_key.shares.iter()
            .find(|s| s.subset == *target_subset)
            .unwrap_or_else(|| panic!(
                "Party {} does not hold share for subset {:?}",
                party_key.party_id, target_subset
            ));

        s1_part = poly::polyvec_add(&s1_part, &share.s1);
        s2_part = poly::polyvec_add(&s2_part, &share.s2);
    }

    (s1_part, s2_part)
}

/// ---
///
/// ---
/// ---
pub fn to_mldsa_public_key(vk: &VerificationKey) -> crate::signature::PublicKey {
    crate::signature::PublicKey {
        rho: vk.rho,
        t1: vk.t1.clone(),
        packed: vk.packed.clone(),
    }
}

// ---

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87};
    use super::threshold_params::{lookup, rss_recover};

    // ---

    #[test]
    fn test_keygen_produces_valid_pk_size() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x42; 32], &tp);

        assert_eq!(bundle.vk.packed.len(), ML_DSA_44.pk_bytes);
        assert_eq!(bundle.vk.t1.len(), ML_DSA_44.k);
        assert_eq!(bundle.party_keys.len(), 3);
    }

    #[test]
    fn test_keygen_deterministic() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let b1 = keygen(&[0x01; 32], &tp);
        let b2 = keygen(&[0x01; 32], &tp);

        assert_eq!(b1.vk.packed, b2.vk.packed);
        assert_eq!(b1.vk.tr, b2.vk.tr);
    }

    #[test]
    fn test_keygen_different_seeds_differ() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let b1 = keygen(&[0x01; 32], &tp);
        let b2 = keygen(&[0x02; 32], &tp);

        assert_ne!(b1.vk.packed, b2.vk.packed);
    }

    #[test]
    fn test_keygen_t1_range() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0xAA; 32], &tp);

        for (k, poly) in bundle.vk.t1.iter().enumerate() {
            for i in 0..POLY_N {
                assert!(poly.coeffs[i] < 1024,
                    "t1[{}][{}] = {} doesn't fit in 10 bits", k, i, poly.coeffs[i]);
            }
        }
    }

    // ---

    #[test]
    fn test_party_share_counts() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x10; 32], &tp);

        // ---
        assert_eq!(bundle.party_keys[0].shares.len(), 2);
        assert_eq!(bundle.party_keys[1].shares.len(), 2);
        assert_eq!(bundle.party_keys[2].shares.len(), 2);
    }

    #[test]
    fn test_party_shares_have_correct_subsets() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x20; 32], &tp);

        // ---
        let subsets_0: Vec<&Vec<usize>> = bundle.party_keys[0].shares.iter()
            .map(|s| &s.subset)
            .collect();
        assert!(subsets_0.contains(&&vec![0, 1]));
        assert!(subsets_0.contains(&&vec![0, 2]));
    }

    #[test]
    fn test_shares_have_bounded_coefficients() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x30; 32], &tp);
        let eta = ML_DSA_44.eta;

        for party in &bundle.party_keys {
            for share in &party.shares {
                for p in &share.s1 {
                    assert!(p.check_norm(eta + 1),
                        "s1 coefficients exceed η={}", eta);
                }
                for p in &share.s2 {
                    assert!(p.check_norm(eta + 1),
                        "s2 coefficients exceed η={}", eta);
                }
            }
        }
    }

    // ---

    #[test]
    fn test_full_secret_sum_matches_pk() {
        // ---
        // ---
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x40; 32], &tp);
        let base = tp.base;

        // ---
        let all_subsets = threshold_params::enumerate_subsets(tp.n_parties, tp.subset_size());
        let mut s1_full = poly::polyvec_zero(base.l);
        let mut s2_full = poly::polyvec_zero(base.k);

        // ---
        // ---
        for subset in &all_subsets {
            // ---
            let party = subset[0]; // ---
            let share = bundle.party_keys[party].shares.iter()
                .find(|s| &s.subset == subset)
                .unwrap();
            s1_full = poly::polyvec_add(&s1_full, &share.s1);
            s2_full = poly::polyvec_add(&s2_full, &share.s2);
        }

        // ---
        let a_hat = sampling::expand_a(&bundle.vk.rho, base);
        let s1_hat = poly::polyvec_ntt(&s1_full);
        let t_hat = poly::mat_vec_mul_ntt(&a_hat, &s1_hat);
        let t_normal = poly::polyvec_inv_ntt(&t_hat);
        let t = poly::polyvec_add(&t_normal, &s2_full);

        // ---
        for i in 0..base.k {
            for j in 0..POLY_N {
                let (hi, _) = arithmetic::power2round(t[i].coeffs[j]);
                assert_eq!(hi, bundle.vk.t1[i].coeffs[j],
                    "t1 mismatch at [{i}][{j}]");
            }
        }
    }

    // ---

    #[test]
    fn test_partial_secrets_sum_to_full() {
        // ---
        // ---
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x50; 32], &tp);
        let base = tp.base;

        let act = vec![0, 2]; // ---
        let partition = rss_recover(tp.n_parties, tp.threshold, &act);

        // ---
        let (s1_part_0, s2_part_0) = partial_secret(
            &bundle.party_keys[0], &partition[0], &tp);
        let (s1_part_2, s2_part_2) = partial_secret(
            &bundle.party_keys[2], &partition[2], &tp);

        // ---
        let s1_sum = poly::polyvec_add(&s1_part_0, &s1_part_2);
        let s2_sum = poly::polyvec_add(&s2_part_0, &s2_part_2);

        // ---
        let all_subsets = threshold_params::enumerate_subsets(tp.n_parties, tp.subset_size());
        let mut s1_full = poly::polyvec_zero(base.l);
        let mut s2_full = poly::polyvec_zero(base.k);
        for subset in &all_subsets {
            let party = subset[0];
            let share = bundle.party_keys[party].shares.iter()
                .find(|s| &s.subset == subset).unwrap();
            s1_full = poly::polyvec_add(&s1_full, &share.s1);
            s2_full = poly::polyvec_add(&s2_full, &share.s2);
        }

        // ---
        for j in 0..base.l {
            for i in 0..POLY_N {
                assert_eq!(s1_sum[j].coeffs[i], s1_full[j].coeffs[i],
                    "s1 partial sum mismatch at [{j}][{i}]");
            }
        }
        for j in 0..base.k {
            for i in 0..POLY_N {
                assert_eq!(s2_sum[j].coeffs[i], s2_full[j].coeffs[i],
                    "s2 partial sum mismatch at [{j}][{i}]");
            }
        }
    }

    #[test]
    fn test_partial_secrets_all_active_sets() {
        // ---
        let tp = lookup(4, 3, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x60; 32], &tp);
        let base = tp.base;

        // ---
        let all_subsets = threshold_params::enumerate_subsets(tp.n_parties, tp.subset_size());
        let mut s1_full = poly::polyvec_zero(base.l);
        let mut s2_full = poly::polyvec_zero(base.k);
        for subset in &all_subsets {
            let party = subset[0];
            let share = bundle.party_keys[party].shares.iter()
                .find(|s| &s.subset == subset).unwrap();
            s1_full = poly::polyvec_add(&s1_full, &share.s1);
            s2_full = poly::polyvec_add(&s2_full, &share.s2);
        }

        // ---
        let active_sets = threshold_params::enumerate_subsets(4, 3);
        for act in &active_sets {
            let partition = rss_recover(tp.n_parties, tp.threshold, act);

            let mut s1_sum = poly::polyvec_zero(base.l);
            let mut s2_sum = poly::polyvec_zero(base.k);

            for &party in act {
                let (s1_p, s2_p) = partial_secret(
                    &bundle.party_keys[party], &partition[party], &tp);
                s1_sum = poly::polyvec_add(&s1_sum, &s1_p);
                s2_sum = poly::polyvec_add(&s2_sum, &s2_p);
            }

            for j in 0..base.l {
                for i in 0..POLY_N {
                    assert_eq!(s1_sum[j].coeffs[i], s1_full[j].coeffs[i],
                        "act={:?} s1 mismatch at [{j}][{i}]", act);
                }
            }
            for j in 0..base.k {
                for i in 0..POLY_N {
                    assert_eq!(s2_sum[j].coeffs[i], s2_full[j].coeffs[i],
                        "act={:?} s2 mismatch at [{j}][{i}]", act);
                }
            }
        }
    }

    // ---

    #[test]
    fn test_vk_converts_to_mldsa_pk() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0x70; 32], &tp);

        let pk = to_mldsa_public_key(&bundle.vk);
        assert_eq!(pk.packed.len(), ML_DSA_44.pk_bytes);
        assert_eq!(pk.rho, bundle.vk.rho);
        assert_eq!(pk.t1.len(), ML_DSA_44.k);
    }

    // ---

    #[test]
    fn test_keygen_ml_dsa_65() {
        let tp = lookup(2, 2, &ML_DSA_65).unwrap();
        let bundle = keygen(&[0x80; 32], &tp);

        assert_eq!(bundle.vk.packed.len(), ML_DSA_65.pk_bytes);
        assert_eq!(bundle.vk.t1.len(), ML_DSA_65.k);
        assert_eq!(bundle.party_keys.len(), 2);

        // ---
        for pk in &bundle.party_keys {
            assert_eq!(pk.shares.len(), 1);
        }
    }

    #[test]
    fn test_keygen_ml_dsa_87() {
        let tp = lookup(3, 2, &ML_DSA_87).unwrap();
        let bundle = keygen(&[0x90; 32], &tp);

        assert_eq!(bundle.vk.packed.len(), ML_DSA_87.pk_bytes);
        assert_eq!(bundle.vk.t1.len(), ML_DSA_87.k);
        assert_eq!(bundle.party_keys.len(), 3);
    }

    // ---

    #[test]
    fn test_keygen_5_3() {
        let tp = lookup(5, 3, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0xA0; 32], &tp);

        assert_eq!(bundle.party_keys.len(), 5);
        // ---
        // ---
        for pk in &bundle.party_keys {
            assert_eq!(pk.shares.len(), 6);
        }

        // ---
        let tr0 = bundle.party_keys[0].tr;
        for pk in &bundle.party_keys {
            assert_eq!(pk.tr, tr0);
        }
    }

    #[test]
    fn test_partial_secret_5_3() {
        let tp = lookup(5, 3, &ML_DSA_44).unwrap();
        let bundle = keygen(&[0xB0; 32], &tp);
        let base = tp.base;

        // ---
        let all_subsets = threshold_params::enumerate_subsets(tp.n_parties, tp.subset_size());
        let mut s1_full = poly::polyvec_zero(base.l);
        let mut s2_full = poly::polyvec_zero(base.k);
        for subset in &all_subsets {
            let party = subset[0];
            let share = bundle.party_keys[party].shares.iter()
                .find(|s| &s.subset == subset).unwrap();
            s1_full = poly::polyvec_add(&s1_full, &share.s1);
            s2_full = poly::polyvec_add(&s2_full, &share.s2);
        }

        // ---
        let act = vec![1, 3, 4];
        let partition = rss_recover(5, 3, &act);

        let mut s1_sum = poly::polyvec_zero(base.l);
        let mut s2_sum = poly::polyvec_zero(base.k);
        for &party in &act {
            let (s1_p, s2_p) = partial_secret(
                &bundle.party_keys[party], &partition[party], &tp);
            s1_sum = poly::polyvec_add(&s1_sum, &s1_p);
            s2_sum = poly::polyvec_add(&s2_sum, &s2_p);
        }

        for j in 0..base.l {
            for i in 0..POLY_N {
                assert_eq!(s1_sum[j].coeffs[i], s1_full[j].coeffs[i]);
            }
        }
    }
}