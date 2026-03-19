//! Combine partial signatures into a standard ML-DSA signature (Paper Figure 7).
//!
//! The combiner assembles partial signatures from T parties into a valid
//! ML-DSA signature (c̃, z^(1), h) that passes standard `signature::verify()`.
//!
//! ## Combine Algorithm (Paper Figure 7)
//!
//! For each parallel instance j where all T parties accepted:
//!   1. z^(1) = Σ_{i∈act} z^(1)_{i,j}
//!   2. δ = w_j - (A·z^(1) - 2^d·c_j·t₁)
//!   3. h = MakeHint(δ, A·z^(1) - 2^d·c_j·t₁, 2γ₂)
//!   4. If ||z^(1)||∞ ≥ γ₁-β or ||δ||∞ > γ₂ or ||h||₁ > ω: try next instance
//!   5. Return sig = (c̃_j, z^(1), h)
//!
//! The output is a **standard ML-DSA signature** — verified by unchanged
//! `signature::verify()`.

use crate::arithmetic;
use crate::encoding;
use crate::params::{N as POLY_N, D};
use crate::poly::{self, Poly};
use crate::sampling;
use crate::shake;
use crate::signature::Signature;

use super::params::ThresholdParams;
use super::protocol::{Round2Msg, Round3Msg};
use super::rss::VerificationKey;

pub fn combine(
    vk: &VerificationKey,
    _act: &[usize],
    msg: &[u8],
    round2_messages: &[Round2Msg],
    round3_messages: &[Round3Msg],
    tp: &ThresholdParams,
) -> Option<Signature> {
    let base = tp.base;
    let k_instances = tp.num_instances;


    let a_hat = sampling::expand_a(&vk.rho, base);


    let mut mu_xof = shake::Shake256::new();
    mu_xof.absorb(&vk.tr);
    mu_xof.absorb(msg);
    let mu = mu_xof.squeeze_vec(64);


    for j in 0..k_instances {
    
        let mut all_accepted = true;
        for r3 in round3_messages {
            if r3.z1_values[j].is_none() {
                all_accepted = false;
                break;
            }
        }
        if !all_accepted {
            continue;
        }

    
        let mut z1 = poly::polyvec_zero(base.l);
        for r3 in round3_messages {
            let z1_ij = r3.z1_values[j].as_ref().unwrap();
            z1 = poly::polyvec_add(&z1, z1_ij);
        }

    
        let mut w_j = poly::polyvec_zero(base.k);
        for r2 in round2_messages {
            w_j = poly::polyvec_add(&w_j, &r2.w_values[j]);
        }

    
        let mut w1_j = Vec::with_capacity(base.k);
        for i in 0..base.k {
            let mut w1_poly = Poly::zero();
            for c in 0..POLY_N {
                w1_poly.coeffs[c] = arithmetic::high_bits(w_j[i].coeffs[c], base.gamma2);
            }
            w1_j.push(w1_poly);
        }

    
        let mut c_hash = shake::Shake256::new();
        c_hash.absorb(&mu);
        let w1_bits = base.w1_bits();
        for poly in &w1_j {
            c_hash.absorb(&encoding::simple_bit_pack(poly, w1_bits));
        }
        let c_tilde_j = c_hash.squeeze_vec(base.c_tilde_bytes);

    
        let c_j = sampling::sample_in_ball(&c_tilde_j, base.tau);
        let c_j_hat = poly::ntt(&c_j);

    
        if !poly::polyvec_check_norm(&z1, base.gamma1 - base.beta) {
            continue;
        }

    
    
        let z1_hat = poly::polyvec_ntt(&z1);
        let az1_hat = poly::mat_vec_mul_ntt(&a_hat, &z1_hat);

        let two_d = 1u32 << D;
        let mut ct1_2d_hat = Vec::with_capacity(base.k);
        for i in 0..base.k {
            let t1_scaled = vk.t1[i].scalar_mul(two_d);
            let t1_hat = poly::ntt(&t1_scaled);
            ct1_2d_hat.push(c_j_hat.pointwise_mul(&t1_hat));
        }

    
        let approx_hat = poly::polyvec_sub(&az1_hat, &ct1_2d_hat);
        let approx = poly::polyvec_inv_ntt(&approx_hat);

    
        let delta = poly::polyvec_sub(&w_j, &approx);

    
        if !poly::polyvec_check_norm(&delta, base.gamma2 + 1) {
            continue;
        }

    
        let mut h = vec![vec![false; POLY_N]; base.k];
        let mut hint_count = 0usize;
        for i in 0..base.k {
            for c in 0..POLY_N {
                let delta_centered = arithmetic::to_centered(delta[i].coeffs[c]);
                h[i][c] = arithmetic::make_hint(
                    delta_centered,
                    approx[i].coeffs[c],
                    base.gamma2,
                );
                if h[i][c] {
                    hint_count += 1;
                }
            }
        }

    
        if hint_count > base.omega {
            continue;
        }

    
        return Some(Signature {
            c_tilde: c_tilde_j,
            z: z1,
            h,
        });
    }

    None // ---
}

// ---

/// ---
///
/// ---
/// ---
/// at most (1/2)^10 ≈ 0.1%, which is negligible.
const MAX_SIGNING_ATTEMPTS: usize = 10;

/// ---
/// ---
///
/// ---
/// ---
/// ---
///
/// ---
/// ---
pub fn threshold_sign(
    vk: &VerificationKey,
    party_keys: &[&super::rss::PartyKey],
    act: &[usize],
    msg: &[u8],
    tp: &ThresholdParams,
    base_seed: &[u8],
) -> Option<Signature> {
    for attempt in 0..MAX_SIGNING_ATTEMPTS {
    
        let mut attempt_xof = shake::Shake256::new();
        attempt_xof.absorb(base_seed);
        attempt_xof.absorb(&(attempt as u32).to_le_bytes());
        let attempt_seed = attempt_xof.squeeze_vec(32);

        if let Some(sig) = threshold_sign_single_attempt(
            vk, party_keys, act, msg, tp, &attempt_seed,
        ) {
            return Some(sig);
        }
    }
    None
}

/// ---
fn threshold_sign_single_attempt(
    vk: &VerificationKey,
    party_keys: &[&super::rss::PartyKey],
    act: &[usize],
    msg: &[u8],
    tp: &ThresholdParams,
    attempt_seed: &[u8],
) -> Option<Signature> {
    use super::protocol::{share_sign_1, share_sign_2, share_sign_3};

    let t = act.len();


    let mut states: Vec<super::protocol::SigningState> = Vec::with_capacity(t);
    let mut r1_msgs: Vec<super::protocol::Round1Msg> = Vec::with_capacity(t);

    for &party_id in act {
        let pk = party_keys.iter().find(|k| k.party_id == party_id).unwrap();
        let mut seed_xof = shake::Shake256::new();
        seed_xof.absorb(attempt_seed);
        seed_xof.absorb(&(party_id as u32).to_le_bytes());
        let party_seed = seed_xof.squeeze_vec(32);

        let (state, r1_msg) = share_sign_1(vk, pk, tp, &party_seed);
        states.push(state);
        r1_msgs.push(r1_msg);
    }


    let mut r2_msgs: Vec<Round2Msg> = Vec::with_capacity(t);
    for idx in 0..t {
        let r1_clones: Vec<_> = r1_msgs.iter().map(|m| m.clone_msg()).collect();
        let r2 = share_sign_2(
            &mut states[idx], vk, act, msg, r1_clones,
        )?;
        r2_msgs.push(r2);
    }


    let mut r3_msgs: Vec<Round3Msg> = Vec::with_capacity(t);
    let r2_refs: Vec<Round2Msg> = r2_msgs.iter().map(|m| m.clone_msg()).collect();
    for idx in 0..t {
        let party_id = act[idx];
        let pk = party_keys.iter().find(|k| k.party_id == party_id).unwrap();
        let r3 = share_sign_3(
            &states[idx], vk, pk, &r2_refs,
        )?;
        r3_msgs.push(r3);
    }


    combine(vk, act, msg, &r2_msgs, &r3_msgs, tp)
}

// ---

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::ML_DSA_44;
    use crate::signature;
    use super::super::params::lookup;
    use super::super::rss;

    /// ---
    /// ---
    #[test]
    fn test_threshold_sign_verify_2_2_ml_dsa_44() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x01; 32], &tp);

        let act = vec![0, 1];
        let msg = b"Hello, threshold ML-DSA!";
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle.vk, &party_keys, &act, msg, &tp, b"signing session 1",
        ).expect("Threshold signing should succeed");

    
        let pk = rss::to_mldsa_public_key(&bundle.vk);
        let valid = signature::verify(&pk, msg, &sig, &ML_DSA_44);
        assert!(valid, "Threshold signature must pass standard ML-DSA verification");
    }

    #[test]
    fn test_threshold_sign_verify_3_2_ml_dsa_44() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x02; 32], &tp);

    
        let act = vec![0, 2];
        let msg = b"3-party threshold test, 2-of-3 signing";
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle.vk, &party_keys, &act, msg, &tp, b"session 3-2",
        ).expect("2-of-3 signing should succeed");

        let pk = rss::to_mldsa_public_key(&bundle.vk);
        assert!(signature::verify(&pk, msg, &sig, &ML_DSA_44));
    }

    #[test]
    fn test_threshold_sign_verify_3_3_ml_dsa_44() {
        let tp = lookup(3, 3, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x03; 32], &tp);

        let act = vec![0, 1, 2];
        let msg = b"3-of-3 full threshold";
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle.vk, &party_keys, &act, msg, &tp, b"session 3-3",
        ).expect("3-of-3 signing should succeed");

        let pk = rss::to_mldsa_public_key(&bundle.vk);
        assert!(signature::verify(&pk, msg, &sig, &ML_DSA_44));
    }

    #[test]
    fn test_threshold_sig_rejected_by_wrong_message() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x04; 32], &tp);

        let act = vec![0, 1];
        let msg = b"correct message";
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle.vk, &party_keys, &act, msg, &tp, b"session wrong",
        ).expect("Signing should succeed");

        let pk = rss::to_mldsa_public_key(&bundle.vk);
        assert!(!signature::verify(&pk, b"wrong message", &sig, &ML_DSA_44),
            "Threshold sig should be rejected for wrong message");
    }

    #[test]
    fn test_threshold_sig_rejected_by_wrong_key() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle1 = rss::keygen(&[0x05; 32], &tp);
        let bundle2 = rss::keygen(&[0x06; 32], &tp);

        let act = vec![0, 1];
        let msg = b"wrong key test";
        let party_keys: Vec<&rss::PartyKey> = bundle1.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle1.vk, &party_keys, &act, msg, &tp, b"session wk",
        ).expect("Signing should succeed");

        let wrong_pk = rss::to_mldsa_public_key(&bundle2.vk);
        assert!(!signature::verify(&wrong_pk, msg, &sig, &ML_DSA_44),
            "Threshold sig should be rejected for wrong public key");
    }

    #[test]
    fn test_threshold_sig_z_norm_bound() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x07; 32], &tp);

        let act = vec![0, 1];
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle.vk, &party_keys, &act, b"norm test", &tp, b"session norm",
        ).expect("Signing should succeed");

        let bound = ML_DSA_44.gamma1 - ML_DSA_44.beta;
        assert!(poly::polyvec_check_norm(&sig.z, bound),
            "z should satisfy ||z||∞ < γ₁ - β");
    }

    #[test]
    fn test_threshold_sig_hint_weight() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x08; 32], &tp);

        let act = vec![0, 1];
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle.vk, &party_keys, &act, b"hint test", &tp, b"session hint",
        ).expect("Signing should succeed");

        let hint_count: usize = sig.h.iter()
            .flat_map(|v| v.iter())
            .filter(|&&b| b)
            .count();
        assert!(hint_count <= ML_DSA_44.omega,
            "Hint weight {} exceeds ω={}", hint_count, ML_DSA_44.omega);
    }

    #[test]
    fn test_threshold_sign_empty_message() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x09; 32], &tp);

        let act = vec![0, 1];
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();

        let sig = threshold_sign(
            &bundle.vk, &party_keys, &act, b"", &tp, b"session empty",
        ).expect("Empty message signing should succeed");

        let pk = rss::to_mldsa_public_key(&bundle.vk);
        assert!(signature::verify(&pk, b"", &sig, &ML_DSA_44));
    }

    #[test]
    fn test_threshold_sign_different_active_sets() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x0A; 32], &tp);
        let msg = b"different active sets";
        let party_keys: Vec<&rss::PartyKey> = bundle.party_keys.iter().collect();
        let pk = rss::to_mldsa_public_key(&bundle.vk);

    
        let active_sets = vec![vec![0, 1], vec![0, 2], vec![1, 2]];
        for (idx, act) in active_sets.iter().enumerate() {
            let mut seed = b"session diff ".to_vec();
            seed.push(idx as u8);

            let sig = threshold_sign(
                &bundle.vk, &party_keys, act, msg, &tp, &seed,
            ).expect(&format!("Signing with act={:?} should succeed", act));

            assert!(signature::verify(&pk, msg, &sig, &ML_DSA_44),
                "Signature from act={:?} should verify", act);
        }
    }
}