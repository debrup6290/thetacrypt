//! Three-round threshold signing protocol (Paper Figure 6).
//!
//! Implements ShareSign1/2/3 — the per-party signing state machine with K
//! parallel instances to boost success probability.
//!
//! ## Protocol Flow (per signing attempt, K parallel instances)
//!
//! ```text
//! Round 1 — Commit:
//!   Each party i, for each instance j ∈ [K]:
//!     1. Sample randomness r_{i,j} ← χ_r (hyperball)
//!     2. Compute commitment w_{i,j} = A·r1_{i,j} + r2_{i,j}
//!     3. Publish cmt_{i,j} = H_cmt(vk, i, j, w_{i,j})
//!
//! Round 2 — Reveal:
//!   Each party i reveals all K commitment values w_{i,j}.
//!   All parties verify: cmt_{i,j} == H_cmt(vk, i, j, w_{i,j}).
//!
//! Round 3 — Respond:
//!   For each instance j:
//!     1. Aggregate w_j = Σ_{i∈act} w_{i,j}
//!     2. w1_j = HighBits(w_j, 2γ₂)
//!     3. c̃_j = H(μ || w1_j), c_j = SampleInBall(c̃_j)
//!     4. m = RSSRecover(act)
//!     5. s_part_i = Σ_{I∈m_i} s_I
//!     6. z_{i,j} = HRej(c_j · s_part_i, r, r', ν, M; r_{i,j})
//!     7. If z_{i,j} = ⊥: mark instance j as failed for party i
//!   Output: z^(1)_{i,j} for each non-failed instance j
//! ```

use crate::arithmetic;
use crate::encoding;
use crate::params::N as POLY_N;
use crate::poly::{self, Poly, PolyVec};
use crate::sampling;
use crate::shake;

use super::hyperball::{self, HyperballSample};
use super::params::{self as threshold_params, ThresholdParams};
use super::rss::{self, PartyKey, VerificationKey};


pub struct Round1Msg {
    pub party_id: usize,
    
    pub commitments: Vec<Vec<u8>>,
}


pub struct Round2Msg {
    pub party_id: usize,
    
    pub w_values: Vec<PolyVec>,
}


pub struct Round3Msg {
    pub party_id: usize,
    
    pub z1_values: Vec<Option<PolyVec>>,
}




pub struct SigningState {
    pub party_id: usize,
    pub tp: ThresholdParams,

    
    
    samples: Vec<HyperballSample>,
    
    w_values: Vec<PolyVec>,
    
    #[allow(dead_code)]
    own_commitments: Vec<Vec<u8>>,

    
    
    all_round1: Vec<Round1Msg>,
    
    msg: Vec<u8>,
    
    act: Vec<usize>,
}

impl SigningState {
    /// Overwrite the active set.  Used by the network protocol layer to
    /// replace the placeholder written by `share_sign_2` with the real
    /// T-sized active set computed after round 2.
    pub fn set_act(&mut self, act: Vec<usize>) {
        self.act = act;
    }
}

///


///


pub fn share_sign_1(
    vk: &VerificationKey,
    party_key: &PartyKey,
    tp: &ThresholdParams,
    signing_seed: &[u8],
) -> (SigningState, Round1Msg) {
    let base = tp.base;
    let k_instances = tp.num_instances;

    
    let a_hat = sampling::expand_a(&vk.rho, base);

    let mut samples = Vec::with_capacity(k_instances);
    let mut w_values = Vec::with_capacity(k_instances);
    let mut commitments = Vec::with_capacity(k_instances);

    for j in 0..k_instances {
        
        let mut seed_xof = shake::Shake256::new();
        seed_xof.absorb(signing_seed);
        seed_xof.absorb(&(party_key.party_id as u32).to_le_bytes());
        seed_xof.absorb(&(j as u32).to_le_bytes());
        let instance_seed = seed_xof.squeeze_vec(64);

        
        let sample = hyperball::sample_randomness(&instance_seed, tp);

        
        
        let ri_1_hat = poly::polyvec_ntt(&sample.ri_1);
        let w_hat = poly::mat_vec_mul_ntt(&a_hat, &ri_1_hat);
        let w_normal = poly::polyvec_inv_ntt(&w_hat);
        let w_ij = poly::polyvec_add(&w_normal, &sample.ri_2);

        
        let cmt = commitment_hash(&vk.packed, party_key.party_id, j, &w_ij, base);

        samples.push(sample);
        w_values.push(w_ij);
        commitments.push(cmt);
    }

    let state = SigningState {
        party_id: party_key.party_id,
        tp: *tp,
        samples,
        w_values: w_values.clone(),
        own_commitments: commitments.clone(),
        all_round1: Vec::new(),
        msg: Vec::new(),
        act: Vec::new(),
    };

    let msg = Round1Msg {
        party_id: party_key.party_id,
        commitments,
    };

    (state, msg)
}

pub fn share_sign_2(
    state: &mut SigningState,
    _vk: &VerificationKey,
    act: &[usize],
    msg: &[u8],
    round1_messages: Vec<Round1Msg>,
) -> Option<Round2Msg> {
    state.all_round1 = round1_messages;
    state.msg = msg.to_vec();
    state.act = act.to_vec();

    Some(Round2Msg {
        party_id: state.party_id,
        w_values: state.w_values.clone(),
    })
}

pub fn share_sign_3(
    state: &SigningState,
    vk: &VerificationKey,
    party_key: &PartyKey,
    round2_messages: &[Round2Msg],
) -> Option<Round3Msg> {
    let tp = &state.tp;
    let base = tp.base;
    let k_instances = tp.num_instances;
    let act = &state.act;

    
    for r2 in round2_messages {
        
        let r1 = state.all_round1.iter()
            .find(|m| m.party_id == r2.party_id)?;

        for j in 0..k_instances {
            let expected = commitment_hash(
                &vk.packed, r2.party_id, j, &r2.w_values[j], base,
            );
            if expected != r1.commitments[j] {
                return None; 
            }
        }
    }

    
    let partition = threshold_params::rss_recover(tp.n_parties, tp.threshold, act);

    
    let (s1_part, s2_part) = rss::partial_secret(
        party_key, &partition[state.party_id], tp,
    );

    
    let s1_part_hat = poly::polyvec_ntt(&s1_part);
    let s2_part_hat = poly::polyvec_ntt(&s2_part);

    
    let mut mu_xof = shake::Shake256::new();
    mu_xof.absorb(&vk.tr);
    mu_xof.absorb(&state.msg);
    let mu = mu_xof.squeeze_vec(64);

    let mut z1_values = Vec::with_capacity(k_instances);

    for j in 0..k_instances {
        
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

        
        let mut v1 = Vec::with_capacity(base.l);
        for p in 0..base.l {
            v1.push(poly::inv_ntt(&c_j_hat.pointwise_mul(&s1_part_hat[p])));
        }
        let mut v2 = Vec::with_capacity(base.k);
        for p in 0..base.k {
            v2.push(poly::inv_ntt(&c_j_hat.pointwise_mul(&s2_part_hat[p])));
        }

        
        let z1 = hyperball::hrej(&v1, &v2, &state.samples[j], tp);
        z1_values.push(z1);
    }

    Some(Round3Msg {
        party_id: state.party_id,
        z1_values,
    })
}


fn commitment_hash(
    vk_packed: &[u8],
    party_id: usize,
    instance: usize,
    w: &PolyVec,
    base: &crate::params::Params,
) -> Vec<u8> {
    let mut xof = shake::Shake256::new();
    xof.absorb(vk_packed);
    xof.absorb(&(party_id as u32).to_le_bytes());
    xof.absorb(&(instance as u32).to_le_bytes());
    
    
    for poly in w {
        xof.absorb(&encoding::simple_bit_pack(poly, 23));
    }
    xof.squeeze_vec(base.c_tilde_bytes)
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::ML_DSA_44;
    use super::super::params::lookup;
    use super::super::rss;

    fn setup_2_2() -> (ThresholdParams, rss::ThresholdKeyBundle) {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let bundle = rss::keygen(&[0x42; 32], &tp);
        (tp, bundle)
    }

    #[test]
    fn test_round1_produces_k_commitments() {
        let (tp, bundle) = setup_2_2();

        let (state, r1_msg) = share_sign_1(
            &bundle.vk, &bundle.party_keys[0], &tp, b"signing seed 0",
        );

        assert_eq!(r1_msg.commitments.len(), tp.num_instances);
        assert_eq!(state.w_values.len(), tp.num_instances);
        assert_eq!(state.samples.len(), tp.num_instances);

        
        for cmt in &r1_msg.commitments {
            assert_eq!(cmt.len(), ML_DSA_44.c_tilde_bytes);
        }
    }

    #[test]
    fn test_round1_deterministic() {
        let (tp, bundle) = setup_2_2();

        let (_, r1a) = share_sign_1(
            &bundle.vk, &bundle.party_keys[0], &tp, b"same seed",
        );
        let (_, r1b) = share_sign_1(
            &bundle.vk, &bundle.party_keys[0], &tp, b"same seed",
        );

        assert_eq!(r1a.commitments, r1b.commitments);
    }

    #[test]
    fn test_round1_different_parties_differ() {
        let (tp, bundle) = setup_2_2();

        let (_, r1_0) = share_sign_1(
            &bundle.vk, &bundle.party_keys[0], &tp, b"shared seed",
        );
        let (_, r1_1) = share_sign_1(
            &bundle.vk, &bundle.party_keys[1], &tp, b"shared seed",
        );

        
        assert_ne!(r1_0.commitments[0], r1_1.commitments[0]);
    }

    #[test]
    fn test_round2_reveals_w_values() {
        let (tp, bundle) = setup_2_2();

        let (mut state0, r1_0) = share_sign_1(
            &bundle.vk, &bundle.party_keys[0], &tp, b"seed0",
        );
        let (_, r1_1) = share_sign_1(
            &bundle.vk, &bundle.party_keys[1], &tp, b"seed1",
        );

        let act = vec![0, 1];
        let r2 = share_sign_2(
            &mut state0, &bundle.vk, &act, b"test message",
            vec![r1_0, r1_1],
        ).unwrap();

        assert_eq!(r2.w_values.len(), tp.num_instances);
        assert_eq!(r2.party_id, 0);
    }

    #[test]
    fn test_full_3_round_protocol() {
        let (tp, bundle) = setup_2_2();
        let act = vec![0, 1];
        let msg = b"threshold signing test message";

        
        
        let mut succeeded = false;
        for attempt in 0u8..10 {
            let mut seed0 = b"seed_party_0_attempt_".to_vec();
            seed0.push(attempt);
            let mut seed1 = b"seed_party_1_attempt_".to_vec();
            seed1.push(attempt);

            
            let (mut state0, r1_0) = share_sign_1(
                &bundle.vk, &bundle.party_keys[0], &tp, &seed0,
            );
            let (mut state1, r1_1) = share_sign_1(
                &bundle.vk, &bundle.party_keys[1], &tp, &seed1,
            );

            
            let r2_0 = share_sign_2(
                &mut state0, &bundle.vk, &act, msg,
                vec![r1_0.clone_msg(), r1_1.clone_msg()],
            ).unwrap();
            let r2_1 = share_sign_2(
                &mut state1, &bundle.vk, &act, msg,
                vec![r1_0.clone_msg(), r1_1.clone_msg()],
            ).unwrap();

            
            let r3_0 = share_sign_3(
                &state0, &bundle.vk, &bundle.party_keys[0],
                &[r2_0.clone_msg(), r2_1.clone_msg()],
            ).unwrap();
            let r3_1 = share_sign_3(
                &state1, &bundle.vk, &bundle.party_keys[1],
                &[r2_0.clone_msg(), r2_1.clone_msg()],
            ).unwrap();

            
            for j in 0..tp.num_instances {
                if r3_0.z1_values[j].is_some() && r3_1.z1_values[j].is_some() {
                    succeeded = true;
                    break;
                }
            }
            if succeeded { break; }
        }
        assert!(succeeded, "Protocol should succeed within 10 attempts");
    }
}


// (In a real network protocol these would be serialized/deserialized,


impl Round1Msg {
    
    pub fn clone_msg(&self) -> Round1Msg {
        Round1Msg {
            party_id: self.party_id,
            commitments: self.commitments.clone(),
        }
    }
}

impl Round2Msg {
    
    pub fn clone_msg(&self) -> Round2Msg {
        Round2Msg {
            party_id: self.party_id,
            w_values: self.w_values.iter().map(|v|
                v.iter().map(|p| p.clone()).collect()
            ).collect(),
        }
    }
}

impl Round3Msg {
    
    pub fn clone_msg(&self) -> Round3Msg {
        Round3Msg {
            party_id: self.party_id,
            z1_values: self.z1_values.iter().map(|opt|
                opt.as_ref().map(|v| v.iter().map(|p| p.clone()).collect())
            ).collect(),
        }
    }
}
