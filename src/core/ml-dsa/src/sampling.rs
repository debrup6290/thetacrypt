//! Sampling routines for ML-DSA per FIPS 204.
//!
//! These functions expand seeds into the structured randomness that ML-DSA needs:
//!
//! **Uniform sampling** (for the public matrix A):
//!   - `RejNTTPoly` (Algorithm 30): rejection-sample a uniform polynomial mod q
//!     from a SHAKE-128 stream, using `CoeffFromThreeBytes` (Algorithm 14).
//!   - `ExpandA` (Algorithm 32): expand seed ρ into the k×l matrix A in NTT domain.
//!
//! **Short coefficient sampling** (for secret vectors s1, s2):
//!   - `RejBoundedPoly` (Algorithm 31): rejection-sample a polynomial with
//!     coefficients in [-η, η] from a SHAKE-256 stream, using
//!     `CoeffFromHalfByte` (Algorithm 15).
//!   - `ExpandS` (Algorithm 33): expand seed ρ' into secret vectors s1, s2.
//!
//! **Masking vector** (for signing):
//!   - `ExpandMask` (Algorithm 34): expand seed ρ'' into a masking vector y
//!     with coefficients in [-(γ₁-1), γ₁-1].
//!
//! **Challenge polynomial** (for signing/verification):
//!   - `SampleInBall` (Algorithm 29): sample a polynomial c with exactly τ
//!     non-zero coefficients, each ±1.

use crate::arithmetic;
use crate::encoding;
use crate::params::{N, Q, Params};
use crate::poly::{Poly, PolyVec, PolyMat};
use crate::shake::{Shake128, Shake256};

//  CoeffFromThreeBytes 
// FIPS 204, Algorithm 14.
//
// Given three bytes b0, b1, b2, compute a candidate coefficient:
//   d = b0 + 256·b1 + 65536·(b2 mod 128)
//
// Accept d if d < q; otherwise reject and try the next three bytes.
// This gives a uniform distribution over [0, q).

/// Attempt to extract a coefficient from three bytes.
/// Returns Some(d) if d < q, None otherwise.
#[inline]
fn coeff_from_three_bytes(b0: u8, b1: u8, b2: u8) -> Option<u32> {
    let d = (b0 as u32) | ((b1 as u32) << 8) | (((b2 & 0x7F) as u32) << 16);
    if d < Q { Some(d) } else { None }
}

//  CoeffFromHalfByte 
// FIPS 204, Algorithm 15.
//
// Given a half-byte (4-bit value) z and parameter η, return a
// coefficient in [-η, η] if z is in the valid range, or reject.
//
// For η = 2: accept z < 15, return z mod 5 mapped to [-2, 2]
//            (values 0→0, 1→1, 2→2, 3→-2, 4→-1)
// For η = 4: accept z < 9, return z mapped to [0, 8] → [-4, 4]
//            (value z → z if z ≤ 4, else z - 9)

/// Attempt to extract a coefficient in [-η, η] from a 4-bit value.
/// Returns Some(coeff) as a standard representative in [0, q), or None if rejected.
#[inline]
fn coeff_from_half_byte(z: u8, eta: u32) -> Option<u32> {
    match eta {
        2 => {
            if z < 15 {
                // FIPS 204 Algorithm 15: return 2 - (b mod 5)
                // Maps: 0→2, 1→1, 2→0, 3→-1, 4→-2, 5→2, ...
                let centered = 2 - (z % 5) as i32;
                Some(arithmetic::from_centered(centered))
            } else {
                None // z = 15: reject
            }
        }
        4 => {
            if z < 9 {
                // FIPS 204 Algorithm 15: return 4 - b
                // Maps: 0→4, 1→3, ..., 4→0, ..., 8→-4
                let centered = 4 - z as i32;
                Some(arithmetic::from_centered(centered))
            } else {
                None // z ≥ 9: reject
            }
        }
        _ => panic!("Unsupported eta: {}", eta),
    }
}

//  RejNTTPoly 
// FIPS 204, Algorithm 30.
//
// Sample a polynomial with coefficients uniform in [0, q) by
// rejection-sampling from a SHAKE-128 XOF stream.

/// Rejection-sample a uniform polynomial from a SHAKE-128 stream.
///
/// Reads 3 bytes at a time, applies CoeffFromThreeBytes, and collects
/// 256 accepted coefficients.
pub fn rej_ntt_poly(xof: &mut Shake128) -> Poly {
    let mut p = Poly::zero();
    let mut j = 0usize;

    while j < N {
        let mut buf = [0u8; 3];
        xof.squeeze(&mut buf);

        if let Some(d) = coeff_from_three_bytes(buf[0], buf[1], buf[2]) {
            p.coeffs[j] = d;
            j += 1;
        }
    }

    p
}

//  RejBoundedPoly 
// FIPS 204, Algorithm 31.
//
// Sample a polynomial with coefficients in [-η, η] by rejection-sampling
// from a SHAKE-256 XOF stream.  Each byte provides two candidates
// (low nibble and high nibble).

/// Rejection-sample a short polynomial from a SHAKE-256 stream.
///
/// Reads one byte at a time, extracts two 4-bit candidates (low nibble
/// first, then high nibble), and accepts each via CoeffFromHalfByte.
pub fn rej_bounded_poly(xof: &mut Shake256, eta: u32) -> Poly {
    let mut p = Poly::zero();
    let mut j = 0usize;

    while j < N {
        let mut buf = [0u8; 1];
        xof.squeeze(&mut buf);
        let byte = buf[0];

        // Low nibble
        let z0 = byte & 0x0F;
        if let Some(c) = coeff_from_half_byte(z0, eta) {
            p.coeffs[j] = c;
            j += 1;
            if j >= N { break; }
        }

        // High nibble
        let z1 = byte >> 4;
        if let Some(c) = coeff_from_half_byte(z1, eta) {
            p.coeffs[j] = c;
            j += 1;
        }
    }

    p
}

//  ExpandA 
// FIPS 204, Algorithm 32.
//
// Expand the 32-byte seed ρ into a k × l matrix of polynomials,
// each uniform in NTT domain.
//
// For each entry (i, j):
//   A_hat[i][j] = RejNTTPoly(SHAKE-128(ρ || IntegerToByte(j) || IntegerToByte(i)))
//
// Note: the FIPS 204 convention absorbs j first, then i (column-major indexing).

/// Expand seed ρ into the public matrix A in NTT domain.
///
/// Returns a k×l matrix where each entry is a polynomial with
/// coefficients uniformly distributed in [0, q).
pub fn expand_a(rho: &[u8; 32], params: &Params) -> PolyMat {
    let mut mat = Vec::with_capacity(params.k);

    for i in 0..params.k {
        let mut row = Vec::with_capacity(params.l);
        for j in 0..params.l {
            let mut xof = Shake128::new();
            xof.absorb(rho);
            // FIPS 204: absorb IntegerToBytes(s, 1) || IntegerToBytes(r, 1)
            // where s = j (column) and r = i (row)
            xof.absorb(&[j as u8, i as u8]);
            row.push(rej_ntt_poly(&mut xof));
        }
        mat.push(row);
    }

    mat
}

//  ExpandS 
// FIPS 204, Algorithm 33.
//
// Expand the 64-byte seed ρ' into secret vectors s1 (length l) and
// s2 (length k), with coefficients in [-η, η].
//
// For each vector entry r:
//   s1[r] = RejBoundedPoly(SHAKE-256(ρ' || IntegerToBytes(r, 2)), η)
//   s2[r] = RejBoundedPoly(SHAKE-256(ρ' || IntegerToBytes(l + r, 2)), η)
//
// The counter is encoded as a 2-byte little-endian integer.

/// Expand seed ρ' into secret vectors s1 and s2.
pub fn expand_s(rho_prime: &[u8; 64], params: &Params) -> (PolyVec, PolyVec) {
    let mut s1 = Vec::with_capacity(params.l);
    let mut s2 = Vec::with_capacity(params.k);

    for r in 0..params.l {
        let mut xof = Shake256::new();
        xof.absorb(rho_prime);
        xof.absorb(&(r as u16).to_le_bytes());
        s1.push(rej_bounded_poly(&mut xof, params.eta));
    }

    for r in 0..params.k {
        let mut xof = Shake256::new();
        xof.absorb(rho_prime);
        xof.absorb(&((params.l + r) as u16).to_le_bytes());
        s2.push(rej_bounded_poly(&mut xof, params.eta));
    }

    (s1, s2)
}

//  ExpandMask 
// FIPS 204, Algorithm 34.
//
// Expand the 64-byte seed ρ'' and counter κ into a masking vector y
// of length l, with coefficients in [-(γ₁-1), γ₁-1].
//
// For each entry r:
//   y[r] = BitUnpack(SHAKE-256(ρ'' || IntegerToBytes(κ + r, 2)), γ₁-1, γ₁)
//
// The coefficients are packed using γ₁_bits bits each:
//   γ₁ = 2^17 → 18 bits per coeff → 576 bytes per poly
//   γ₁ = 2^19 → 20 bits per coeff → 640 bytes per poly

/// Expand seed ρ'' and counter κ into a masking vector y.
pub fn expand_mask(rho_pp: &[u8; 64], kappa: u16, params: &Params) -> PolyVec {
    let gamma1 = params.gamma1;
    let gamma1_bits = params.gamma1_bits();
    let byte_count = N * gamma1_bits as usize / 8;

    let mut y = Vec::with_capacity(params.l);

    for r in 0..params.l {
        let mut xof = Shake256::new();
        xof.absorb(rho_pp);
        xof.absorb(&(kappa + r as u16).to_le_bytes());

        // Squeeze enough bytes for one polynomial
        let buf = xof.squeeze_vec(byte_count);

        // Unpack using BitUnpack(buf, γ₁-1, γ₁)
        // This maps unsigned values in [0, 2γ₁-1] to centered [-(γ₁-1), γ₁]
        let p = encoding::bit_unpack(&buf, gamma1 - 1, gamma1);
        y.push(p);
    }

    y
}

//  SampleInBall 
// FIPS 204, Algorithm 29.
//
// Sample a challenge polynomial c with exactly τ non-zero coefficients,
// each equal to ±1, and the remaining coefficients zero.
//
// Algorithm:
//   1. Absorb the seed c̃ into SHAKE-256.
//   2. Squeeze 8 bytes → 64-bit sign mask.
//   3. For i from 256-τ to 255 (i.e., the last τ positions):
//      a. Sample j uniformly from [0, i] by rejection:
//         squeeze one byte; reject if > i.
//      b. Swap: c[i] ← c[j]; c[j] ← ±1 (sign from the sign mask).
//
// This is a Fisher-Yates shuffle restricted to the last τ positions.

/// Sample a challenge polynomial with exactly τ non-zero ±1 coefficients.
pub fn sample_in_ball(seed: &[u8], tau: usize) -> Poly {
    let mut c = Poly::zero();
    let mut xof = Shake256::from_data(seed);

    // Step 1: squeeze 8 bytes for the sign bits
    let mut sign_bytes = [0u8; 8];
    xof.squeeze(&mut sign_bytes);
    let signs = u64::from_le_bytes(sign_bytes);

    // Step 2: Fisher-Yates for positions (256-τ) through 255
    for i in (N - tau)..N {
        // Sample j uniformly from [0, i] by rejection
        let j = loop {
            let mut buf = [0u8; 1];
            xof.squeeze(&mut buf);
            let candidate = buf[0] as usize;
            if candidate <= i {
                break candidate;
            }
        };

        // Swap c[i] ← c[j]
        c.coeffs[i] = c.coeffs[j];

        // Set c[j] = ±1 based on the sign bit
        let bit_idx = i - (N - tau);
        if (signs >> bit_idx) & 1 == 1 {
            c.coeffs[j] = Q - 1; // -1 mod q
        } else {
            c.coeffs[j] = 1;
        }
    }

    c
}

//  Tests 

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87};

    // ── CoeffFromThreeBytes ──

    #[test]
    fn test_coeff_from_three_bytes_accepts_zero() {
        assert_eq!(coeff_from_three_bytes(0, 0, 0), Some(0));
    }

    #[test]
    fn test_coeff_from_three_bytes_accepts_below_q() {
        // q - 1 = 8380416 = 0x7FE001
        // b0 = 0x01, b1 = 0xE0, b2 = 0x7F → d = 0x7FE001 = 8380417-1? Let me compute:
        // Actually q = 8380417 = 0x7FE001. So q-1 = 0x7FE000
        // b0 = 0x00, b1 = 0xE0, b2 = 0x7F → d = 0 + 0xE0*256 + 0x7F*65536
        //   = 57344 + 5177344 = ... let me just test boundary
        assert_eq!(coeff_from_three_bytes(0, 0, 1), Some(65536));
        assert!(coeff_from_three_bytes(0, 0, 1).unwrap() < Q);
    }

    #[test]
    fn test_coeff_from_three_bytes_rejects_at_q() {
        // q = 8380417 = 0x7FE001
        // b0=0x01, b1=0xE0, b2=0x7F → d = 1 + 57344 + 0x7F*65536 = 1+57344+8323072 = 8380417 = q
        assert_eq!(coeff_from_three_bytes(0x01, 0xE0, 0x7F), None);
    }

    #[test]
    fn test_coeff_from_three_bytes_masks_high_bit() {
        // b2 with bit 7 set: should be masked to b2 & 0x7F
        // b2 = 0xFF → masked to 0x7F
        // d = 0xFF + 0xFF*256 + 0x7F*65536 = 255 + 65280 + 8323072 = 8388607 > q → rejected
        assert_eq!(coeff_from_three_bytes(0xFF, 0xFF, 0xFF), None);
        // b2 = 0x80 → masked to 0x00
        assert_eq!(coeff_from_three_bytes(0, 0, 0x80), Some(0));
    }

    // ── CoeffFromHalfByte ──

    #[test]
    fn test_coeff_from_half_byte_eta2() {
        // η=2: FIPS 204 Algorithm 15: return 2 - (b mod 5)
        // 0→2, 1→1, 2→0, 3→-1, 4→-2
        assert_eq!(coeff_from_half_byte(0, 2), Some(arithmetic::from_centered(2)));
        assert_eq!(coeff_from_half_byte(1, 2), Some(arithmetic::from_centered(1)));
        assert_eq!(coeff_from_half_byte(2, 2), Some(0));
        assert_eq!(coeff_from_half_byte(3, 2), Some(arithmetic::from_centered(-1)));
        assert_eq!(coeff_from_half_byte(4, 2), Some(arithmetic::from_centered(-2)));
        // z=5: 5 mod 5 = 0 → 2
        assert_eq!(coeff_from_half_byte(5, 2), Some(arithmetic::from_centered(2)));
        // z=14: 14 mod 5 = 4 → -2
        assert_eq!(coeff_from_half_byte(14, 2), Some(arithmetic::from_centered(-2)));
        // z=15: rejected
        assert_eq!(coeff_from_half_byte(15, 2), None);
    }

    #[test]
    fn test_coeff_from_half_byte_eta4() {
        // η=4: z < 9 accepted, centered = 4 - z
        assert_eq!(coeff_from_half_byte(0, 4), Some(arithmetic::from_centered(4)));
        assert_eq!(coeff_from_half_byte(4, 4), Some(0));
        assert_eq!(coeff_from_half_byte(8, 4), Some(arithmetic::from_centered(-4)));
        // z=9: rejected
        assert_eq!(coeff_from_half_byte(9, 4), None);
        assert_eq!(coeff_from_half_byte(15, 4), None);
    }

    // ── SampleInBall ──

    #[test]
    fn test_sample_in_ball_weight() {
        // c should have exactly τ non-zero coefficients
        for (seed_val, tau) in [(0u8, 39usize), (1, 49), (2, 60)] {
            let seed = [seed_val; 32];
            let c = sample_in_ball(&seed, tau);

            let nonzero = c.coeffs.iter().filter(|&&x| x != 0).count();
            assert_eq!(nonzero, tau,
                "SampleInBall(seed={}, τ={}) has {} nonzero (expected {})",
                seed_val, tau, nonzero, tau);
        }
    }

    #[test]
    fn test_sample_in_ball_values() {
        // All non-zero coefficients must be ±1 (i.e., 1 or q-1 mod q)
        let c = sample_in_ball(&[42; 32], 39);
        for &coeff in &c.coeffs {
            assert!(
                coeff == 0 || coeff == 1 || coeff == Q - 1,
                "Coefficient should be 0, 1, or -1 mod q, got {}", coeff
            );
        }
    }

    #[test]
    fn test_sample_in_ball_deterministic() {
        let seed = [0xAB; 32];
        let c1 = sample_in_ball(&seed, 39);
        let c2 = sample_in_ball(&seed, 39);
        for i in 0..N {
            assert_eq!(c1.coeffs[i], c2.coeffs[i],
                "SampleInBall should be deterministic");
        }
    }

    #[test]
    fn test_sample_in_ball_different_seeds() {
        let c1 = sample_in_ball(&[1; 32], 39);
        let c2 = sample_in_ball(&[2; 32], 39);
        // They should differ (astronomically unlikely to be equal)
        let same = (0..N).all(|i| c1.coeffs[i] == c2.coeffs[i]);
        assert!(!same, "Different seeds should produce different challenges");
    }

    #[test]
    fn test_sample_in_ball_all_tau_values() {
        // Test with all three ML-DSA τ values
        for &tau in &[ML_DSA_44.tau, ML_DSA_65.tau, ML_DSA_87.tau] {
            let c = sample_in_ball(&[0xFF; 32], tau);
            let nonzero = c.coeffs.iter().filter(|&&x| x != 0).count();
            assert_eq!(nonzero, tau, "τ={} weight mismatch", tau);

            // Roughly half should be +1 and half -1
            let pos = c.coeffs.iter().filter(|&&x| x == 1).count();
            let neg = c.coeffs.iter().filter(|&&x| x == Q - 1).count();
            assert_eq!(pos + neg, tau);
        }
    }

    // ── ExpandA ──

    #[test]
    fn test_expand_a_dimensions() {
        let rho = [0u8; 32];
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let a = expand_a(&rho, params);
            assert_eq!(a.len(), params.k, "{}: A should have k={} rows", params.name, params.k);
            for row in &a {
                assert_eq!(row.len(), params.l, "{}: A should have l={} cols", params.name, params.l);
            }
        }
    }

    #[test]
    fn test_expand_a_deterministic() {
        let rho = [0x42; 32];
        let params = &ML_DSA_44;
        let a1 = expand_a(&rho, params);
        let a2 = expand_a(&rho, params);
        for i in 0..params.k {
            for j in 0..params.l {
                for c in 0..N {
                    assert_eq!(a1[i][j].coeffs[c], a2[i][j].coeffs[c],
                        "ExpandA should be deterministic at [{i}][{j}][{c}]");
                }
            }
        }
    }

    #[test]
    fn test_expand_a_uniform() {
        // Coefficients should be uniform in [0, q).
        // Basic check: they shouldn't all be zero, and should cover a wide range.
        let rho = [0x99; 32];
        let a = expand_a(&rho, &ML_DSA_44);

        let all_coeffs: Vec<u32> = a.iter()
            .flat_map(|row| row.iter())
            .flat_map(|p| p.coeffs.iter().copied())
            .collect();

        // All should be < q
        for &c in &all_coeffs {
            assert!(c < Q, "Coefficient {} out of range", c);
        }

        // Check distribution: max should be close to q, min close to 0
        let max_c = *all_coeffs.iter().max().unwrap();
        let min_c = *all_coeffs.iter().min().unwrap();
        assert!(max_c > Q / 2, "Max coefficient {} is suspiciously small", max_c);
        assert!(min_c < Q / 2, "Min coefficient {} is suspiciously large", min_c);
    }

    #[test]
    fn test_expand_a_different_indices_differ() {
        let rho = [0; 32];
        let a = expand_a(&rho, &ML_DSA_44);

        // A[0][0] and A[0][1] should be different polynomials
        let same = (0..N).all(|c| a[0][0].coeffs[c] == a[0][1].coeffs[c]);
        assert!(!same, "Different matrix entries should differ");

        // A[0][0] and A[1][0] should also differ
        let same2 = (0..N).all(|c| a[0][0].coeffs[c] == a[1][0].coeffs[c]);
        assert!(!same2, "Different rows should differ");
    }

    // ── ExpandS ──

    #[test]
    fn test_expand_s_dimensions() {
        let rho_prime = [0u8; 64];
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let (s1, s2) = expand_s(&rho_prime, params);
            assert_eq!(s1.len(), params.l, "{}: s1 length", params.name);
            assert_eq!(s2.len(), params.k, "{}: s2 length", params.name);
        }
    }

    #[test]
    fn test_expand_s_bound() {
        // All coefficients should be in [-η, η]
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let rho_prime = [0x55; 64];
            let (s1, s2) = expand_s(&rho_prime, params);

            for (label, vec) in [("s1", &s1), ("s2", &s2)] {
                for (j, poly) in vec.iter().enumerate() {
                    for i in 0..N {
                        let centered = arithmetic::to_centered(poly.coeffs[i]);
                        assert!(
                            centered.abs() <= params.eta as i32,
                            "{}: {}[{}][{}] = {} out of [-{}, {}]",
                            params.name, label, j, i, centered, params.eta, params.eta
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_expand_s_deterministic() {
        let rho_prime = [0xAA; 64];
        let params = &ML_DSA_44;
        let (s1_a, s2_a) = expand_s(&rho_prime, params);
        let (s1_b, s2_b) = expand_s(&rho_prime, params);

        for j in 0..params.l {
            for i in 0..N {
                assert_eq!(s1_a[j].coeffs[i], s1_b[j].coeffs[i]);
            }
        }
        for j in 0..params.k {
            for i in 0..N {
                assert_eq!(s2_a[j].coeffs[i], s2_b[j].coeffs[i]);
            }
        }
    }

    #[test]
    fn test_expand_s_different_vectors_differ() {
        let rho_prime = [0; 64];
        let (s1, s2) = expand_s(&rho_prime, &ML_DSA_44);

        // s1[0] and s1[1] should differ
        let same = (0..N).all(|i| s1[0].coeffs[i] == s1[1].coeffs[i]);
        assert!(!same, "s1[0] and s1[1] should differ");

        // s1[0] and s2[0] should differ (different counters)
        let same2 = (0..N).all(|i| s1[0].coeffs[i] == s2[0].coeffs[i]);
        assert!(!same2, "s1[0] and s2[0] should differ");
    }

    #[test]
    fn test_expand_s_distribution() {
        // Coefficients should cover the full range [-η, η]
        let rho_prime = [0x77; 64];
        let (s1, _) = expand_s(&rho_prime, &ML_DSA_44);

        let mut seen = std::collections::HashSet::new();
        for poly in &s1 {
            for &c in &poly.coeffs {
                seen.insert(arithmetic::to_centered(c));
            }
        }

        // For η=2, we should see all of {-2, -1, 0, 1, 2}
        for v in -2i32..=2 {
            assert!(seen.contains(&v), "Missing value {} in s1 (η=2)", v);
        }
    }

    // ── ExpandMask ──

    #[test]
    fn test_expand_mask_dimensions() {
        let rho_pp = [0u8; 64];
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let y = expand_mask(&rho_pp, 0, params);
            assert_eq!(y.len(), params.l, "{}: y length", params.name);
        }
    }

    #[test]
    fn test_expand_mask_bound() {
        // Coefficients should be in [-(γ₁-1), γ₁]
        // (actually the range from BitUnpack(γ₁-1, γ₁) is [-(γ₁-1), γ₁])
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let rho_pp = [0xDD; 64];
            let y = expand_mask(&rho_pp, 0, params);

            let gamma1 = params.gamma1 as i32;
            for (j, poly) in y.iter().enumerate() {
                for i in 0..N {
                    let centered = arithmetic::to_centered(poly.coeffs[i]);
                    assert!(
                        centered >= -(gamma1 - 1) && centered <= gamma1,
                        "{}: y[{}][{}] = {} out of [{}, {}]",
                        params.name, j, i, centered, -(gamma1 - 1), gamma1
                    );
                }
            }
        }
    }

    #[test]
    fn test_expand_mask_deterministic() {
        let rho_pp = [0xEE; 64];
        let params = &ML_DSA_44;
        let y1 = expand_mask(&rho_pp, 0, params);
        let y2 = expand_mask(&rho_pp, 0, params);
        for j in 0..params.l {
            for i in 0..N {
                assert_eq!(y1[j].coeffs[i], y2[j].coeffs[i]);
            }
        }
    }

    #[test]
    fn test_expand_mask_different_kappa() {
        // Different κ values should produce different masks
        let rho_pp = [0; 64];
        let params = &ML_DSA_44;
        let y0 = expand_mask(&rho_pp, 0, params);
        let y1 = expand_mask(&rho_pp, params.l as u16, params);

        let same = (0..N).all(|i| y0[0].coeffs[i] == y1[0].coeffs[i]);
        assert!(!same, "Different κ should produce different masks");
    }
}
