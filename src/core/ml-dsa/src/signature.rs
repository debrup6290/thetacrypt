//! ML-DSA KeyGen, Sign, and Verify — FIPS 204 Algorithms 1–3.
//!
//! This module implements the three top-level operations of ML-DSA:
//!
//! - **KeyGen** (Algorithm 1): Generate a key pair (pk, sk) from a 32-byte seed ξ.
//! - **Sign** (Algorithm 2): Sign a message, producing a signature (deterministic variant).
//! - **Verify** (Algorithm 3): Verify a signature against a public key and message.
//!
//! ## Deterministic vs. Hedged Signing
//!
//! FIPS 204 defines two signing modes:
//!   - **Deterministic** (ML-DSA.Sign_internal with rnd = 0…0):
//!     ρ'' = H(K || rnd || μ) where rnd is 32 zero bytes.
//!   - **Hedged** (randomized): rnd is 32 random bytes.
//!
//! This implementation provides both: `sign` (deterministic) and
//! `sign_randomized` (hedged, if you supply your own randomness).

use crate::arithmetic;
use crate::encoding;
use crate::params::{N, D, Params};
use crate::poly::{self, Poly, PolyVec};
use crate::sampling;
use crate::shake;

//  Data Types 

/// ML-DSA key pair.
pub struct KeyPair {
    pub pk: PublicKey,
    pub sk: SecretKey,
}

/// Public key.
///
/// Contains the seed ρ, the high bits t1, and the packed byte encoding
/// (cached for hashing during signing/verification).
pub struct PublicKey {
    pub rho: [u8; 32],
    pub t1: PolyVec,
    /// The packed public key bytes (ρ || t1_packed).
    pub packed: Vec<u8>,
}

/// Secret key.
///
/// Contains everything needed for signing:
///   - ρ: seed for expanding A
///   - K: secret seed for deterministic signing
///   - tr: hash of the public key (tr = H(pk))
///   - s1, s2: short secret vectors
///   - t0: low bits from Power2Round
pub struct SecretKey {
    pub rho: [u8; 32],
    pub k_seed: [u8; 32],
    pub tr: [u8; 64],
    pub s1: PolyVec,
    pub s2: PolyVec,
    pub t0: PolyVec,
}

/// ML-DSA signature.
pub struct Signature {
    pub c_tilde: Vec<u8>,
    pub z: PolyVec,
    pub h: Vec<Vec<bool>>,
}


impl PublicKey {
    /// Reconstruct a PublicKey from its FIPS 204 packed byte encoding.
    pub fn from_packed(packed: Vec<u8>, params: &Params) -> Option<Self> {
        if packed.len() != params.pk_bytes { return None; }
        let mut rho = [0u8; 32];
        rho.copy_from_slice(&packed[0..32]);
        let bytes_per_t1_poly = params.t1_packed_bytes(); // 320
        let mut t1 = Vec::with_capacity(params.k);
        for i in 0..params.k {
            let start = 32 + i * bytes_per_t1_poly;
            t1.push(encoding::unpack_t1(&packed[start..start + bytes_per_t1_poly]));
        }
        Some(Self { rho, t1, packed })
    }
}

impl Signature {
    /// Encode this signature to FIPS 204 bytes.
    pub fn to_bytes(&self, params: &Params) -> Vec<u8> {
        let mut out = Vec::with_capacity(params.sig_bytes);
        // c̃
        out.extend_from_slice(&self.c_tilde);
        // z: l polynomials
        let bits = params.gamma1_bits();
        for poly in &self.z {
            out.extend(encoding::encode_z(poly, params.gamma1, bits));
        }
        // h
        out.extend(encoding::encode_hint(&self.h, params));
        out
    }

    /// Decode a signature from FIPS 204 bytes.
    pub fn from_bytes(bytes: &[u8], params: &Params) -> Option<Self> {
        if bytes.len() != params.sig_bytes { return None; }
        let mut pos = 0;

        let c_tilde = bytes[pos..pos + params.c_tilde_bytes].to_vec();
        pos += params.c_tilde_bytes;

        let bits = params.gamma1_bits();
        let bytes_per_z_poly = 256 * bits as usize / 8;
        let mut z = Vec::with_capacity(params.l);
        for _ in 0..params.l {
            z.push(encoding::decode_z(&bytes[pos..pos + bytes_per_z_poly], params.gamma1, bits));
            pos += bytes_per_z_poly;
        }

        let h = encoding::decode_hint(&bytes[pos..], params)?;
        Some(Self { c_tilde, z, h })
    }
}

//  KeyGen 
// FIPS 204, Algorithm 1: ML-DSA.KeyGen_internal(ξ)
//
// Input:  ξ — 32-byte seed
// Output: (pk, sk)
//
// Steps:
//   1. (ρ, ρ', K) ← H(ξ, 128)             [expand seed via SHAKE-256]
//   2. Â ← ExpandA(ρ)                       [public matrix in NTT domain]
//   3. (s1, s2) ← ExpandS(ρ')              [short secret vectors]
//   4. ŝ1 ← NTT(s1)
//   5. t̂ ← Â · ŝ1                          [matrix-vector multiply in NTT domain]
//   6. t ← NTT⁻¹(t̂) + s2
//   7. (t1, t0) ← Power2Round(t)
//   8. pk ← pkEncode(ρ, t1)
//   9. tr ← H(pk, 64)
//  10. sk ← skEncode(ρ, K, tr, s1, s2, t0)
//  11. return (pk, sk)

/// Generate an ML-DSA key pair from a 32-byte seed ξ.
pub fn keygen(xi: &[u8; 32], params: &Params) -> KeyPair {
    // Step 1: Expand ξ into (ρ, ρ', K)
    let expanded = shake::shake256(xi, 128);
    let mut rho = [0u8; 32];
    let mut rho_prime = [0u8; 64];
    let mut k_seed = [0u8; 32];
    rho.copy_from_slice(&expanded[0..32]);
    rho_prime.copy_from_slice(&expanded[32..96]);
    k_seed.copy_from_slice(&expanded[96..128]);

    // Step 2: Â = ExpandA(ρ)
    let a_hat = sampling::expand_a(&rho, params);

    // Step 3: (s1, s2) = ExpandS(ρ')
    let (s1, s2) = sampling::expand_s(&rho_prime, params);

    // Step 4: ŝ1 = NTT(s1)
    let s1_hat = poly::polyvec_ntt(&s1);

    // Step 5: t̂ = Â · ŝ1
    let t_hat = poly::mat_vec_mul_ntt(&a_hat, &s1_hat);

    // Step 6: t = NTT⁻¹(t̂) + s2
    let t_normal = poly::polyvec_inv_ntt(&t_hat);
    let t = poly::polyvec_add(&t_normal, &s2);

    // Step 7: (t1, t0) = Power2Round(t)
    let mut t1 = Vec::with_capacity(params.k);
    let mut t0 = Vec::with_capacity(params.k);
    for i in 0..params.k {
        let mut t1_poly = Poly::zero();
        let mut t0_poly = Poly::zero();
        for j in 0..N {
            let (hi, lo) = arithmetic::power2round(t[i].coeffs[j]);
            t1_poly.coeffs[j] = hi;
            t0_poly.coeffs[j] = arithmetic::from_centered(lo);
        }
        t1.push(t1_poly);
        t0.push(t0_poly);
    }

    // Step 8: pk = pkEncode(ρ, t1)
    let pk_bytes = encoding::pk_encode(&rho, &t1, params);

    // Step 9: tr = H(pk, 64)
    let tr_vec = shake::shake256(&pk_bytes, 64);
    let mut tr = [0u8; 64];
    tr.copy_from_slice(&tr_vec);

    // Return key pair
    KeyPair {
        pk: PublicKey {
            rho,
            t1,
            packed: pk_bytes,
        },
        sk: SecretKey {
            rho,
            k_seed,
            tr,
            s1,
            s2,
            t0,
        },
    }
}

//  Sign 
// FIPS 204, Algorithm 2: ML-DSA.Sign_internal(sk, M', rnd)
//
// This implements the deterministic variant (rnd = 0^32).
//
// Input:  sk, message M
// Output: signature σ = (c̃, z, h)
//
// Steps:
//   1. Â ← ExpandA(ρ)
//   2. μ ← H(tr || M, 64)
//   3. ρ'' ← H(K || rnd || μ, 64)          [rnd = 0^32 for deterministic]
//   4. ŝ1 ← NTT(s1), ŝ2 ← NTT(s2), t̂0 ← NTT(t0)
//   5. κ ← 0
//   6. (z, h) ← ⊥                          [rejection loop]
//   7. while (z, h) = ⊥:
//      a. y ← ExpandMask(ρ'', κ)
//      b. ŷ ← NTT(y)
//      c. w ← NTT⁻¹(Â · ŷ)
//      d. w1 ← HighBits(w)
//      e. c̃ ← H(μ || w1Encode(w1), λ/4)
//      f. c ← SampleInBall(c̃)
//      g. ĉ ← NTT(c)
//      h. cs1 ← NTT⁻¹(ĉ · ŝ1)
//      i. cs2 ← NTT⁻¹(ĉ · ŝ2)
//      j. z ← y + cs1
//      k. r0 ← LowBits(w - cs2)
//      l. if ||z||∞ ≥ γ₁-β  or  ||r0||∞ ≥ γ₂-β: continue
//      m. ct0 ← NTT⁻¹(ĉ · t̂0)
//      n. if ||ct0||∞ ≥ γ₂: continue
//      o. h ← MakeHint(-ct0, w - cs2 + ct0)
//      p. if count(h) > ω: continue
//      q. κ ← κ + l
//   8. σ ← sigEncode(c̃, z, h)
//   9. return σ

/// Sign a message (deterministic variant, rnd = 0^32).
///
/// Returns None only if the rejection loop exceeds `max_attempts`
/// (should not happen in practice with correct parameters).
pub fn sign(sk: &SecretKey, msg: &[u8], params: &Params) -> Option<Signature> {
    sign_internal(sk, msg, &[0u8; 32], params)
}

/// Sign a message with explicit randomness (hedged variant).
pub fn sign_randomized(sk: &SecretKey, msg: &[u8], rnd: &[u8; 32], params: &Params) -> Option<Signature> {
    sign_internal(sk, msg, rnd, params)
}

fn sign_internal(sk: &SecretKey, msg: &[u8], rnd: &[u8; 32], params: &Params) -> Option<Signature> {
    const MAX_ATTEMPTS: usize = 1000;

    // Step 1: Â = ExpandA(ρ)
    let a_hat = sampling::expand_a(&sk.rho, params);

    // Step 2: μ = H(tr || M, 64)
    let mut mu_xof = shake::Shake256::new();
    mu_xof.absorb(&sk.tr);
    mu_xof.absorb(msg);
    let mu = mu_xof.squeeze_vec(64);

    // Step 3: ρ'' = H(K || rnd || μ, 64)
    let mut rho_pp_xof = shake::Shake256::new();
    rho_pp_xof.absorb(&sk.k_seed);
    rho_pp_xof.absorb(rnd);
    rho_pp_xof.absorb(&mu);
    let rho_pp_vec = rho_pp_xof.squeeze_vec(64);
    let mut rho_pp = [0u8; 64];
    rho_pp.copy_from_slice(&rho_pp_vec);

    // Step 4: Pre-compute NTT of secret vectors
    let s1_hat = poly::polyvec_ntt(&sk.s1);
    let s2_hat = poly::polyvec_ntt(&sk.s2);
    let t0_hat = poly::polyvec_ntt(&sk.t0);

    // Steps 5–7: Rejection loop
    let mut kappa: u16 = 0;
    for _attempt in 0..MAX_ATTEMPTS {
        // (a) y = ExpandMask(ρ'', κ)
        let y = sampling::expand_mask(&rho_pp, kappa, params);

        // (b) ŷ = NTT(y)
        let y_hat = poly::polyvec_ntt(&y);

        // (c) w = NTT⁻¹(Â · ŷ)
        let w_hat = poly::mat_vec_mul_ntt(&a_hat, &y_hat);
        let w = poly::polyvec_inv_ntt(&w_hat);

        // (d) w1 = HighBits(w)
        let mut w1 = Vec::with_capacity(params.k);
        for i in 0..params.k {
            let mut w1_poly = Poly::zero();
            for j in 0..N {
                w1_poly.coeffs[j] = arithmetic::high_bits(w[i].coeffs[j], params.gamma2);
            }
            w1.push(w1_poly);
        }

        // (e) c̃ = H(μ || w1Encode(w1), c_tilde_bytes)
        let mut c_hash = shake::Shake256::new();
        c_hash.absorb(&mu);
        let w1_bits = params.w1_bits();
        for poly in &w1 {
            c_hash.absorb(&encoding::simple_bit_pack(poly, w1_bits));
        }
        let c_tilde = c_hash.squeeze_vec(params.c_tilde_bytes);

        // (f) c = SampleInBall(c̃)
        let c = sampling::sample_in_ball(&c_tilde, params.tau);

        // (g) ĉ = NTT(c)
        let c_hat = poly::ntt(&c);

        // (h) cs1 = NTT⁻¹(ĉ · ŝ1)  [per component]
        // (i) cs2 = NTT⁻¹(ĉ · ŝ2)  [per component]
        let mut cs1 = Vec::with_capacity(params.l);
        for j in 0..params.l {
            cs1.push(poly::inv_ntt(&c_hat.pointwise_mul(&s1_hat[j])));
        }
        let mut cs2 = Vec::with_capacity(params.k);
        for j in 0..params.k {
            cs2.push(poly::inv_ntt(&c_hat.pointwise_mul(&s2_hat[j])));
        }

        // (j) z = y + cs1
        let z = poly::polyvec_add(&y, &cs1);

        // (k) w - cs2
        let w_minus_cs2 = poly::polyvec_sub(&w, &cs2);

        // (l) Rejection check 1: ||z||∞ < γ₁ - β
        if !poly::polyvec_check_norm(&z, params.gamma1 - params.beta) {
            kappa += params.l as u16;
            continue;
        }

        // (l) Rejection check 2: ||LowBits(w - cs2)||∞ < γ₂ - β
        let mut r0_ok = true;
        for i in 0..params.k {
            for j in 0..N {
                let r0 = arithmetic::low_bits(w_minus_cs2[i].coeffs[j], params.gamma2);
                if r0.unsigned_abs() >= params.gamma2 - params.beta {
                    r0_ok = false;
                    break;
                }
            }
            if !r0_ok { break; }
        }
        if !r0_ok {
            kappa += params.l as u16;
            continue;
        }

        // (m) ct0 = NTT⁻¹(ĉ · t̂0)
        let mut ct0 = Vec::with_capacity(params.k);
        for j in 0..params.k {
            ct0.push(poly::inv_ntt(&c_hat.pointwise_mul(&t0_hat[j])));
        }

        // (n) Rejection check 3: ||ct0||∞ < γ₂
        if !poly::polyvec_check_norm(&ct0, params.gamma2) {
            kappa += params.l as u16;
            continue;
        }

        // (o) h = MakeHint(-ct0, w - cs2 + ct0)
        //     Equivalently: for each coefficient,
        //       MakeHint(neg_ct0, w_minus_cs2 + ct0)
        //     which checks if HighBits changes when adding ct0.
        let mut h = vec![vec![false; N]; params.k];
        let mut hint_count = 0usize;
        for i in 0..params.k {
            for j in 0..N {
                let neg_ct0 = arithmetic::to_centered(arithmetic::neg(ct0[i].coeffs[j]));
                let w_cs2_ct0 = arithmetic::add(w_minus_cs2[i].coeffs[j], ct0[i].coeffs[j]);
                h[i][j] = arithmetic::make_hint(neg_ct0, w_cs2_ct0, params.gamma2);
                if h[i][j] {
                    hint_count += 1;
                }
            }
        }

        // (p) Rejection check 4: number of hints ≤ ω
        if hint_count > params.omega {
            kappa += params.l as u16;
            continue;
        }

        // Success!
        return Some(Signature { c_tilde, z, h });
    }

    None // Should not happen in practice
}

//  Verify 
// FIPS 204, Algorithm 3: ML-DSA.Verify_internal(pk, M', σ)
//
// Input:  pk, message M, signature σ = (c̃, z, h)
// Output: true/false
//
// Steps:
//   1. Â ← ExpandA(ρ)
//   2. μ ← H(tr || M, 64)            where tr = H(pk, 64)
//   3. c ← SampleInBall(c̃)
//   4. ĉ ← NTT(c)
//   5. Check ||z||∞ < γ₁ - β
//   6. ŵ'_approx ← Â·NTT(z) - ĉ·NTT(t1·2^d)
//   7. w'_approx ← NTT⁻¹(ŵ'_approx)
//   8. w1' ← UseHint(h, w'_approx)
//   9. c̃' ← H(μ || w1Encode(w1'), λ/4)
//  10. Check hint weight ≤ ω
//  11. return c̃ == c̃'

/// Verify a signature on a message.
pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature, params: &Params) -> bool {
    // ── Structural validation: reject malformed signatures without panicking ──
    if sig.c_tilde.len() != params.c_tilde_bytes { return false; }
    if sig.z.len() != params.l { return false; }
    if sig.h.len() != params.k { return false; }
    for v in &sig.h {
        if v.len() != N { return false; }
    }
    if pk.t1.len() != params.k { return false; }

    // Step 1: Â = ExpandA(ρ)
    let a_hat = sampling::expand_a(&pk.rho, params);

    // Step 2: μ = H(tr || M, 64)  where tr = H(pk, 64)
    let tr = shake::shake256(&pk.packed, 64);
    let mut mu_xof = shake::Shake256::new();
    mu_xof.absorb(&tr);
    mu_xof.absorb(msg);
    let mu = mu_xof.squeeze_vec(64);

    // Step 3: c = SampleInBall(c̃)
    let c = sampling::sample_in_ball(&sig.c_tilde, params.tau);

    // Step 4: ĉ = NTT(c)
    let c_hat = poly::ntt(&c);

    // Step 5: Check ||z||∞ < γ₁ - β
    if !poly::polyvec_check_norm(&sig.z, params.gamma1 - params.beta) {
        return false;
    }

    // Step 6: ŵ'_approx = Â·NTT(z) - ĉ·NTT(t1·2^d)
    let z_hat = poly::polyvec_ntt(&sig.z);
    let az_hat = poly::mat_vec_mul_ntt(&a_hat, &z_hat);

    let two_d = 1u32 << D;
    let mut ct1_2d_hat = Vec::with_capacity(params.k);
    for i in 0..params.k {
        let t1_scaled = pk.t1[i].scalar_mul(two_d);
        let t1_hat = poly::ntt(&t1_scaled);
        ct1_2d_hat.push(c_hat.pointwise_mul(&t1_hat));
    }

    let w_prime_hat = poly::polyvec_sub(&az_hat, &ct1_2d_hat);

    // Step 7: w'_approx = NTT⁻¹(ŵ'_approx)
    let w_prime = poly::polyvec_inv_ntt(&w_prime_hat);

    // Step 8: w1' = UseHint(h, w'_approx)
    let mut w1_prime = Vec::with_capacity(params.k);
    for i in 0..params.k {
        let mut w1p = Poly::zero();
        for j in 0..N {
            w1p.coeffs[j] = arithmetic::use_hint(
                sig.h[i][j],
                w_prime[i].coeffs[j],
                params.gamma2,
            );
        }
        w1_prime.push(w1p);
    }

    // Step 9: c̃' = H(μ || w1Encode(w1'), c_tilde_bytes)
    let mut c_hash = shake::Shake256::new();
    c_hash.absorb(&mu);
    let w1_bits = params.w1_bits();
    for poly in &w1_prime {
        c_hash.absorb(&encoding::simple_bit_pack(poly, w1_bits));
    }
    let c_tilde_prime = c_hash.squeeze_vec(params.c_tilde_bytes);

    // Step 10: Check hint weight ≤ ω
    let hint_count: usize = sig.h.iter()
        .flat_map(|v| v.iter())
        .filter(|&&b| b)
        .count();
    if hint_count > params.omega {
        return false;
    }

    // Step 11: c̃ == c̃'
    c_tilde_prime == sig.c_tilde
}

//  Tests 

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87};

    // ── KeyGen ──

    #[test]
    fn test_keygen_deterministic() {
        let seed = [0x42u8; 32];
        let params = &ML_DSA_44;

        let kp1 = keygen(&seed, params);
        let kp2 = keygen(&seed, params);

        assert_eq!(kp1.pk.packed, kp2.pk.packed, "Same seed → same pk");
        assert_eq!(kp1.sk.tr, kp2.sk.tr, "Same seed → same tr");
        assert_eq!(kp1.sk.k_seed, kp2.sk.k_seed, "Same seed → same K");
    }

    #[test]
    fn test_keygen_pk_size() {
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let kp = keygen(&[0xAA; 32], params);
            assert_eq!(kp.pk.packed.len(), params.pk_bytes,
                "{} pk size", params.name);
        }
    }

    #[test]
    fn test_keygen_different_seeds_differ() {
        let params = &ML_DSA_44;
        let kp1 = keygen(&[1; 32], params);
        let kp2 = keygen(&[2; 32], params);
        assert_ne!(kp1.pk.packed, kp2.pk.packed, "Different seeds → different pk");
    }

    #[test]
    fn test_keygen_t1_range() {
        // t1 coefficients should fit in 10 bits
        let kp = keygen(&[0xBB; 32], &ML_DSA_44);
        for (k, poly) in kp.pk.t1.iter().enumerate() {
            for i in 0..N {
                assert!(poly.coeffs[i] < 1024,
                    "t1[{}][{}] = {} doesn't fit in 10 bits", k, i, poly.coeffs[i]);
            }
        }
    }

    #[test]
    fn test_keygen_s1_s2_bound() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0xCC; 32], params);
        for p in &kp.sk.s1 {
            assert!(p.check_norm(params.eta + 1), "s1 exceeds η bound");
        }
        for p in &kp.sk.s2 {
            assert!(p.check_norm(params.eta + 1), "s2 exceeds η bound");
        }
    }

    // ── Sign + Verify roundtrip ──

    #[test]
    fn test_sign_verify_ml_dsa_44() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x01; 32], params);
        let msg = b"Hello, ML-DSA-44!";

        let sig = sign(&kp.sk, msg, params)
            .expect("Signing should succeed");

        assert!(verify(&kp.pk, msg, &sig, params),
            "Signature should verify");
    }

    #[test]
    fn test_sign_verify_ml_dsa_65() {
        let params = &ML_DSA_65;
        let kp = keygen(&[0x02; 32], params);
        let msg = b"Hello, ML-DSA-65!";

        let sig = sign(&kp.sk, msg, params)
            .expect("Signing should succeed");

        assert!(verify(&kp.pk, msg, &sig, params),
            "Signature should verify");
    }

    #[test]
    fn test_sign_verify_ml_dsa_87() {
        let params = &ML_DSA_87;
        let kp = keygen(&[0x03; 32], params);
        let msg = b"Hello, ML-DSA-87!";

        let sig = sign(&kp.sk, msg, params)
            .expect("Signing should succeed");

        assert!(verify(&kp.pk, msg, &sig, params),
            "Signature should verify");
    }

    // ── Signature properties ──

    #[test]
    fn test_sign_deterministic() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x10; 32], params);
        let msg = b"deterministic test";

        let sig1 = sign(&kp.sk, msg, params).unwrap();
        let sig2 = sign(&kp.sk, msg, params).unwrap();

        assert_eq!(sig1.c_tilde, sig2.c_tilde, "Deterministic sign → same c̃");
    }

    #[test]
    fn test_verify_rejects_wrong_message() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x20; 32], params);
        let msg = b"correct message";

        let sig = sign(&kp.sk, msg, params).unwrap();

        assert!(!verify(&kp.pk, b"wrong message", &sig, params),
            "Should reject wrong message");
    }

    #[test]
    fn test_verify_rejects_wrong_key() {
        let params = &ML_DSA_44;
        let kp1 = keygen(&[0x30; 32], params);
        let kp2 = keygen(&[0x31; 32], params);

        let msg = b"test message";
        let sig = sign(&kp1.sk, msg, params).unwrap();

        assert!(!verify(&kp2.pk, msg, &sig, params),
            "Should reject wrong public key");
    }

    // ── Various messages ──

    #[test]
    fn test_sign_verify_empty_message() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x40; 32], params);

        let sig = sign(&kp.sk, b"", params).unwrap();
        assert!(verify(&kp.pk, b"", &sig, params));
    }

    #[test]
    fn test_sign_verify_long_message() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x50; 32], params);
        let msg = vec![0xAB; 10000];

        let sig = sign(&kp.sk, &msg, params).unwrap();
        assert!(verify(&kp.pk, &msg, &sig, params));
    }

    // ── Signature component checks ──

    #[test]
    fn test_signature_z_bound() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x60; 32], params);
        let sig = sign(&kp.sk, b"z bound test", params).unwrap();

        // z should satisfy ||z||∞ < γ₁ - β
        let bound = params.gamma1 - params.beta;
        assert!(poly::polyvec_check_norm(&sig.z, bound),
            "z should satisfy the norm bound");
    }

    #[test]
    fn test_signature_hint_weight() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x70; 32], params);
        let sig = sign(&kp.sk, b"hint weight test", params).unwrap();

        let hint_count: usize = sig.h.iter()
            .flat_map(|v| v.iter())
            .filter(|&&b| b)
            .count();
        assert!(hint_count <= params.omega,
            "Hint weight {} exceeds ω={}", hint_count, params.omega);
    }

    // ── Robustness: malformed signatures must not panic ──

    #[test]
    fn test_verify_rejects_malformed_signature() {
        let params = &ML_DSA_44;
        let kp = keygen(&[0x80; 32], params);

        // Empty hint vector (wrong dimensions)
        let bad_sig = Signature {
            c_tilde: vec![0u8; params.c_tilde_bytes],
            z: poly::polyvec_zero(params.l),
            h: vec![], // Wrong: should be k vectors of length N
        };
        assert!(!verify(&kp.pk, b"test", &bad_sig, params),
            "Should reject malformed h dimensions");

        // Wrong z length
        let bad_sig2 = Signature {
            c_tilde: vec![0u8; params.c_tilde_bytes],
            z: poly::polyvec_zero(1), // Wrong: should be l
            h: vec![vec![false; N]; params.k],
        };
        assert!(!verify(&kp.pk, b"test", &bad_sig2, params),
            "Should reject wrong z length");

        // Wrong c_tilde length
        let bad_sig3 = Signature {
            c_tilde: vec![0u8; 5], // Wrong length
            z: poly::polyvec_zero(params.l),
            h: vec![vec![false; N]; params.k],
        };
        assert!(!verify(&kp.pk, b"test", &bad_sig3, params),
            "Should reject wrong c_tilde length");
    }
}
