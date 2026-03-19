//! Field arithmetic modulo q = 8380417.
//!
//! This module provides the foundational arithmetic for ML-DSA:
//!
//! **Basic operations** (§4 of FIPS 204):
//!   - Addition, subtraction, multiplication mod q
//!   - Barrett reduction (avoids expensive division)
//!   - Centered representatives: conversion between [0, q) and [-(q-1)/2, (q-1)/2]
//!   - Modular exponentiation (for computing NTT roots)
//!
//! **Rounding and decomposition** (FIPS 204 Algorithms 35–40):
//!   - Power2Round: split r into (r1, r0) where r = r1·2^d + r0
//!   - Decompose:   split r into (r1, r0) where r ≈ r1·α + r0 (α = 2γ₂)
//!   - HighBits / LowBits: extract the respective parts from Decompose
//!   - MakeHint / UseHint: hint mechanism for signature verification

use crate::params::{Q, D};

//  Barrett Reduction 
//
// Barrett reduction avoids division by q.  For x < q², we compute:
//
//   t = x − ⌊x · m / 2^k⌋ · q
//
// where k and m are chosen so that the approximation is tight enough
// to need at most one conditional subtraction.
//
// For q = 8380417:
//   q² < 2^46, so we need k ≥ 46.
//   Choosing k = 48: m = ⌊2^48 / q⌋ = 33554942
//   For any x < 2^46 < 2^48, the result after one correction is exact.

const BARRETT_SHIFT: u64 = 48;
const BARRETT_MULTIPLIER: u64 = (1u64 << BARRETT_SHIFT) / (Q as u64); // 33554942

/// Reduce x mod q using Barrett reduction.
///
/// **Precondition**: x < 2^48 (satisfied by all products of two values < q,
/// since q² ≈ 7.02 × 10¹³ < 2^46).
///
/// **Postcondition**: result is in [0, q).
#[inline(always)]
pub fn reduce(x: u64) -> u32 {
    let quotient = ((x as u128 * BARRETT_MULTIPLIER as u128) >> BARRETT_SHIFT) as u64;
    let mut r = (x - quotient * (Q as u64)) as u32;
    if r >= Q {
        r -= Q;
    }
    r
}

// Basic Field Operations 

/// Modular addition: (a + b) mod q.
///
/// **Precondition**: a, b ∈ [0, q).
/// **Postcondition**: result ∈ [0, q).
#[inline(always)]
pub fn add(a: u32, b: u32) -> u32 {
    let sum = a + b;
    if sum >= Q { sum - Q } else { sum }
}

/// Modular subtraction: (a - b) mod q.
///
/// **Precondition**: a, b ∈ [0, q).
/// **Postcondition**: result ∈ [0, q).
#[inline(always)]
pub fn sub(a: u32, b: u32) -> u32 {
    if a >= b {
        a - b
    } else {
        a + Q - b
    }
}

/// Modular multiplication: (a · b) mod q.
///
/// **Precondition**: a, b ∈ [0, q).
/// **Postcondition**: result ∈ [0, q).
///
/// Uses Barrett reduction on the 64-bit product.
#[inline(always)]
pub fn mul(a: u32, b: u32) -> u32 {
    reduce(a as u64 * b as u64)
}

/// Modular negation: (-a) mod q.
///
/// **Precondition**: a ∈ [0, q).
/// **Postcondition**: result ∈ [0, q).
#[inline(always)]
pub fn neg(a: u32) -> u32 {
    if a == 0 { 0 } else { Q - a }
}

//  Centered Representatives 

/// Convert from standard representative [0, q) to centered [-(q-1)/2, (q-1)/2].
///
/// Values 0 ..= (q-1)/2  map to themselves (non-negative).
/// Values (q-1)/2+1 ..= q-1  map to negative values (subtract q).
///
/// Example: q=8380417
///   0 → 0,  1 → 1,  4190208 → 4190208
///   4190209 → -4190208,  8380416 → -1
#[inline(always)]
pub fn to_centered(a: u32) -> i32 {
    debug_assert!(a < Q, "to_centered: input {} must be < q", a);
    let half = (Q / 2) as i32; // 4190208
    if (a as i32) > half {
        a as i32 - Q as i32
    } else {
        a as i32
    }
}

/// Convert from centered representative back to standard [0, q).
///
/// Inverse of `to_centered`.
#[inline(always)]
pub fn from_centered(a: i32) -> u32 {
    debug_assert!(
        a >= -((Q as i32 - 1) / 2) && a <= (Q as i32 / 2),
        "from_centered: {} out of range", a
    );
    if a < 0 {
        (a + Q as i32) as u32
    } else {
        a as u32
    }
}

//  Modular Exponentiation 

/// Compute base^exp mod q via binary exponentiation.
///
/// Used for:
///   - Computing NTT roots of unity (ζ = 1753 raised to various powers)
///   - Computing modular inverses via Fermat's little theorem: a^{-1} = a^{q-2}
pub fn pow_mod(base: u32, mut exp: u32) -> u32 {
    let mut result: u64 = 1;
    let mut b: u64 = base as u64;
    while exp > 0 {
        if exp & 1 == 1 {
            result = reduce(result * b) as u64;
        }
        b = reduce(b * b) as u64;
        exp >>= 1;
    }
    result as u32
}

/// Compute the modular inverse a^{-1} mod q using Fermat's little theorem.
///
/// Since q is prime: a^{-1} ≡ a^{q-2} (mod q).
///
/// **Precondition**: a ≠ 0.
#[inline]
pub fn inv_mod(a: u32) -> u32 {
    debug_assert!(a != 0, "Cannot invert zero");
    pow_mod(a, Q - 2)
}

// Power2Round
// FIPS 204, Algorithm 35.
//
// Splits r ∈ Z_q into (r1, r0) such that:
//   r ≡ r1 · 2^d + r0  (mod q)
//
// where r0 ∈ (-(2^(d-1) - 1),  2^(d-1)]  (centered mod 2^d).
//
// This is used in KeyGen to split t into t1 (public) and t0 (secret).
// d = 13 always.

/// Power2Round: decompose r into (r1, r0) where r = r1·2^d + r0.
///
/// **Input**: r ∈ [0, q)
/// **Output**: (r1, r0) where:
///   - r0 ∈ [-(2^(d-1)-1), 2^(d-1)]  (a centered remainder)
///   - r1 = (r - r0) / 2^d
///   - r1 · 2^d + r0 = r
pub fn power2round(r: u32) -> (u32, i32) {
    debug_assert!(r < Q, "power2round: input must be < q");

    // r0 = r mod± 2^d  (centered mod 2^d)
    let r0 = centered_mod_power2(r, D);

    // r1 = (r - r0) / 2^d
    let r1 = ((r as i32 - r0) >> D) as u32;

    (r1, r0)
}

/// Centered modulo 2^d: compute r mod± 2^d.
///
/// Result is in [-(2^(d-1) - 1), 2^(d-1)].
///
/// Algorithm:
///   1. r0 = r mod 2^d          (standard remainder in [0, 2^d))
///   2. if r0 > 2^(d-1): r0 -= 2^d   (center it)
fn centered_mod_power2(r: u32, d: u32) -> i32 {
    let two_d = 1i32 << d;           // 2^d = 8192
    let half = 1i32 << (d - 1);      // 2^(d-1) = 4096
    let mask = two_d - 1;            // 2^d - 1 = 8191

    let mut r0 = (r as i32) & mask;  // r mod 2^d, in [0, 2^d)
    if r0 > half {
        r0 -= two_d;
    }
    r0
}

//  Decompose 
// FIPS 204, Algorithm 36.
//
// Splits r ∈ Z_q into (r1, r0) such that:
//   r ≡ r1 · α + r0  (mod q)      (approximately)
//
// where α = 2·γ₂ and r0 ∈ (-α/2, α/2].
//
// Special case: when r - r0 = q - 1 (the "wrap-around"), set r1 = 0
// and decrement r0 by 1.  This ensures r1 stays in [0, m) where
// m = (q-1)/α.
//
// Used during signing: w1 = HighBits(w) determines the commitment,
// and r0 = LowBits(w - cs2) is checked against the rejection bound.

/// Decompose: split r into (r1, r0) where r ≈ r1·α + r0.
///
/// **Input**: r ∈ [0, q),  gamma2 (determines α = 2·gamma2)
/// **Output**: (r1, r0) where:
///   - r0 ∈ (-α/2, α/2]  (note: exclusive lower bound)
///   - r1 ∈ [0, m)  where m = (q-1)/α
///   - r1·α + r0 ≡ r (mod q)  (except the special case)
pub fn decompose(r: u32, gamma2: u32) -> (u32, i32) {
    debug_assert!(r < Q, "decompose: input must be < q");

    let alpha = 2 * gamma2;

    // r0 = r mod± α  (centered modulo α)
    let mut r0 = centered_mod_alpha(r, alpha);

    // Special case: if r - r0 would equal q-1, adjust to avoid overflow
    // of r1 beyond the valid range [0, m).
    if r as i32 - r0 == (Q - 1) as i32 {
        r0 -= 1;
        return (0, r0);
    }

    // General case
    let r1 = ((r as i32 - r0) / alpha as i32) as u32;
    (r1, r0)
}

/// Centered modulo α: compute r mod± α.
///
/// Result is in (-α/2, α/2], i.e., the range is half-open with
/// -α/2 excluded and α/2 included.
///
/// Algorithm:
///   1. r0 = r mod α                (in [0, α))
///   2. if r0 > α/2: r0 -= α        (center it)
fn centered_mod_alpha(r: u32, alpha: u32) -> i32 {
    let half = (alpha / 2) as i32;
    let mut r0 = (r % alpha) as i32;
    if r0 > half {
        r0 -= alpha as i32;
    }
    r0
}

//  HighBits / LowBits 
// FIPS 204, Algorithms 37 and 38.

/// HighBits: extract the high-order representative r1 from Decompose(r).
///
/// This determines which "bin" r falls into, and is used to compute w1
/// (the commitment that gets hashed during signing).
#[inline]
pub fn high_bits(r: u32, gamma2: u32) -> u32 {
    decompose(r, gamma2).0
}

/// LowBits: extract the low-order representative r0 from Decompose(r).
///
/// This is checked against the rejection bound γ₂ - β during signing.
#[inline]
pub fn low_bits(r: u32, gamma2: u32) -> i32 {
    decompose(r, gamma2).1
}

//  MakeHint 
// FIPS 204, Algorithm 39.
//
// MakeHint checks whether adding z to r changes the HighBits.
// During signing, this is used to compute the hint vector h:
//   h[i][j] = MakeHint(-ct0[i][j], w_minus_cs2[i][j] + ct0[i][j])
//
// The hint allows the verifier to recover w1 = HighBits(w) from
// w' = Az - ct1·2^d (which it can compute) without knowing w itself.

/// MakeHint: returns true if HighBits(r) ≠ HighBits(r + z mod q).
///
/// **Inputs**:
///   - z: a value (as centered representative, i.e. i32)
///   - r: a value in [0, q)
///   - gamma2: the decomposition parameter
///
/// **Output**: true if the hint is needed (high bits differ)
pub fn make_hint(z: i32, r: u32, gamma2: u32) -> bool {
    let r1 = high_bits(r, gamma2);

    // Compute (r + z) mod q via centered arithmetic
    let rz = ((r as i64) + (z as i64)).rem_euclid(Q as i64) as u32;
    let v1 = high_bits(rz, gamma2);

    r1 != v1
}

//  UseHint 
// FIPS 204, Algorithm 40.
//
// UseHint recovers the correct HighBits value using the hint.
//
// If h = false (no hint): return r1 as-is.
// If h = true:  adjust r1 by ±1, depending on the sign of r0.
//   - r0 > 0: r1 = (r1 + 1) mod m
//   - r0 ≤ 0: r1 = (r1 - 1) mod m
// where m = (q-1)/α.

/// UseHint: given hint h and value r, recover the correct high bits.
///
/// **Inputs**:
///   - h: the hint bit (from MakeHint)
///   - r: a value in [0, q)
///   - gamma2: the decomposition parameter
///
/// **Output**: the corrected high-bits value in [0, m)
pub fn use_hint(h: bool, r: u32, gamma2: u32) -> u32 {
    let alpha = 2 * gamma2;
    let m = (Q - 1) / alpha;   // number of high-bits "bins"

    let (r1, r0) = decompose(r, gamma2);

    if !h {
        return r1;
    }

    // Adjust r1 by ±1 based on sign of r0
    if r0 > 0 {
        // r0 is positive: the true value is slightly above the bin boundary
        (r1 + 1) % m
    } else {
        // r0 is non-positive: the true value is slightly below
        if r1 == 0 { m - 1 } else { r1 - 1 }
    }
}

//  Tests 

#[cfg(test)]
mod tests {
    use super::*;

    // ── Barrett reduction ──

    #[test]
    fn test_reduce_small_values() {
        assert_eq!(reduce(0), 0);
        assert_eq!(reduce(1), 1);
        assert_eq!(reduce(Q as u64 - 1), Q - 1);
        assert_eq!(reduce(Q as u64), 0);
        assert_eq!(reduce(Q as u64 + 1), 1);
    }

    #[test]
    fn test_reduce_large_products() {
        // Test with products close to q²
        let a = Q - 1;
        let b = Q - 1;
        let product = a as u64 * b as u64; // (q-1)² = q²-2q+1
        let expected = ((product) % (Q as u64)) as u32;
        assert_eq!(reduce(product), expected);
    }

    // ── Basic arithmetic ──

    #[test]
    fn test_add_identity() {
        for &a in &[0, 1, Q / 2, Q - 1] {
            assert_eq!(add(a, 0), a, "add({}, 0)", a);
            assert_eq!(add(0, a), a, "add(0, {})", a);
        }
    }

    #[test]
    fn test_add_wrapping() {
        assert_eq!(add(Q - 1, 1), 0);
        assert_eq!(add(Q - 1, 2), 1);
        assert_eq!(add(Q / 2, Q / 2 + 1), 0);  // q is odd, so q/2 + q/2+1 = q
    }

    #[test]
    fn test_sub_identity() {
        for &a in &[0, 1, Q / 2, Q - 1] {
            assert_eq!(sub(a, 0), a, "sub({}, 0)", a);
            assert_eq!(sub(a, a), 0, "sub({}, {})", a, a);
        }
    }

    #[test]
    fn test_sub_wrapping() {
        assert_eq!(sub(0, 1), Q - 1);
        assert_eq!(sub(0, Q - 1), 1);
    }

    #[test]
    fn test_add_sub_inverse() {
        // (a + b) - b = a for various values
        let pairs = [(0, 0), (1, Q - 1), (12345, 67890), (Q / 2, Q / 2)];
        for (a, b) in pairs {
            assert_eq!(sub(add(a, b), b), a, "add-sub inverse for ({}, {})", a, b);
        }
    }

    #[test]
    fn test_mul_identity() {
        assert_eq!(mul(0, 42), 0);
        assert_eq!(mul(42, 0), 0);
        assert_eq!(mul(1, 42), 42);
        assert_eq!(mul(42, 1), 42);
    }

    #[test]
    fn test_mul_inverse_exists() {
        // (q+1)/2 * 2 ≡ 1 (mod q) because q is odd
        let half = (Q + 1) / 2;
        assert_eq!(mul(half, 2), 1, "(q+1)/2 * 2 should be 1 mod q");
    }

    #[test]
    fn test_mul_commutativity() {
        let pairs = [(123, 456), (Q - 1, Q - 2), (1000000, 7777777 % Q)];
        for (a, b) in pairs {
            assert_eq!(mul(a, b), mul(b, a), "mul({}, {}) should be commutative", a, b);
        }
    }

    #[test]
    fn test_mul_distributive() {
        // a * (b + c) = a*b + a*c
        let (a, b, c) = (12345, 67890, 11111);
        let lhs = mul(a, add(b, c));
        let rhs = add(mul(a, b), mul(a, c));
        assert_eq!(lhs, rhs, "Distributive law failed");
    }

    #[test]
    fn test_neg() {
        assert_eq!(neg(0), 0);
        assert_eq!(neg(1), Q - 1);
        assert_eq!(neg(Q - 1), 1);
        assert_eq!(add(42, neg(42)), 0);
    }

    // ── Centered representatives ──

    #[test]
    fn test_centered_roundtrip() {
        let test_vals = [0, 1, Q / 2, Q / 2 + 1, Q - 1, 42, 1000000];
        for &v in &test_vals {
            let centered = to_centered(v);
            let back = from_centered(centered);
            assert_eq!(back, v, "centered roundtrip failed for {}", v);
        }
    }

    #[test]
    fn test_centered_range() {
        // 0 maps to 0 (non-negative)
        assert_eq!(to_centered(0), 0);
        // q/2 = 4190208 maps to 4190208 (the maximum non-negative value)
        assert_eq!(to_centered(Q / 2), (Q / 2) as i32);
        // q/2 + 1 = 4190209 maps to -4190208 (wraps to negative)
        assert_eq!(to_centered(Q / 2 + 1), -((Q / 2) as i32));
        // q - 1 maps to -1
        assert_eq!(to_centered(Q - 1), -1);
    }

    // ── Modular exponentiation ──

    #[test]
    fn test_pow_mod_basic() {
        assert_eq!(pow_mod(2, 0), 1);
        assert_eq!(pow_mod(2, 1), 2);
        assert_eq!(pow_mod(2, 10), 1024);
        assert_eq!(pow_mod(0, 5), 0);
        assert_eq!(pow_mod(1, 1000000), 1);
    }

    #[test]
    fn test_pow_mod_fermat() {
        // By Fermat's little theorem: a^(q-1) ≡ 1 (mod q) for a ≠ 0
        for &a in &[1u32, 2, 42, 12345, Q - 1] {
            assert_eq!(
                pow_mod(a, Q - 1), 1,
                "Fermat's little theorem failed for a={}", a
            );
        }
    }

    #[test]
    fn test_inv_mod() {
        for &a in &[1u32, 2, 42, 12345, Q - 1, 256] {
            let a_inv = inv_mod(a);
            assert_eq!(
                mul(a, a_inv), 1,
                "a * a^{{-1}} should be 1 for a={}", a
            );
        }
    }

    #[test]
    fn test_inv_mod_256() {
        // 256^{-1} mod q is used in the inverse NTT.
        let n_inv = inv_mod(256);
        assert_eq!(mul(256, n_inv), 1);
        // Verify the actual value
        assert_eq!(n_inv, pow_mod(256, Q - 2));
    }

    // ── Primitive root for NTT ──

    #[test]
    fn test_zeta_is_primitive_512th_root() {
        // ζ = 1753 should be a primitive 512th root of unity.
        // That means: ζ^256 ≡ -1 (mod q) and ζ^512 ≡ 1 (mod q).
        let zeta: u32 = 1753;
        let z256 = pow_mod(zeta, 256);
        let z512 = pow_mod(zeta, 512);

        // ζ^256 = q - 1 = -1 mod q
        assert_eq!(z256, Q - 1, "ζ^256 should be -1 mod q");
        // ζ^512 = 1 mod q
        assert_eq!(z512, 1, "ζ^512 should be 1 mod q");
        // ζ^128 ≠ ±1 (it's primitive, not a lower-order root)
        let z128 = pow_mod(zeta, 128);
        assert_ne!(z128, 1, "ζ should not be a 128th root");
        assert_ne!(z128, Q - 1, "ζ should not be a 256th root at order 128");
    }

    // ── Power2Round ──

    #[test]
    fn test_power2round_reconstruction() {
        // For any r: r1 · 2^d + r0 = r
        let test_vals = [0, 1, Q - 1, 123456, 4000000, 7777777, Q / 2];
        let two_d = 1i64 << D;

        for &r in &test_vals {
            let (r1, r0) = power2round(r);
            let reconstructed = r1 as i64 * two_d + r0 as i64;
            assert_eq!(
                reconstructed as u32, r,
                "Power2Round reconstruction failed for r={}: r1={}, r0={}", r, r1, r0
            );
        }
    }

    #[test]
    fn test_power2round_r0_range() {
        // r0 should be in [-(2^(d-1)-1), 2^(d-1)] = [-4095, 4096]
        let half = 1i32 << (D - 1); // 4096
        let test_vals = [0, 1, Q - 1, Q / 2, 8191, 8192, 8193, 4096, 4095];

        for &r in &test_vals {
            let (_r1, r0) = power2round(r);
            assert!(
                r0 >= -(half - 1) && r0 <= half,
                "Power2Round r0={} out of range [-{}, {}] for r={}",
                r0, half - 1, half, r
            );
        }
    }

    #[test]
    fn test_power2round_r1_range() {
        // r1 should be in [0, (q-1)/2^d] = [0, 1023]
        // Actually the max t1 value determines how many bits we need.
        // (q-1) / 2^13 = 8380416 / 8192 = 1022.xxx → max r1 = 1023
        for &r in &[0, Q - 1, Q / 2] {
            let (r1, _r0) = power2round(r);
            assert!(
                r1 <= 1023,
                "Power2Round r1={} too large for r={}", r1, r
            );
        }
    }

    // ── Decompose ──

    #[test]
    fn test_decompose_reconstruction_gamma2_44() {
        // For ML-DSA-44: γ₂ = 95232, α = 190464
        let gamma2 = (Q - 1) / 88;
        let alpha = 2 * gamma2;
        let test_vals = [0, 1, Q - 1, Q / 2, 4000000, 190463, 190464, 190465];

        for &r in &test_vals {
            let (r1, r0) = decompose(r, gamma2);
            let reconstructed = ((r1 as i64) * (alpha as i64) + (r0 as i64))
                .rem_euclid(Q as i64) as u32;
            assert_eq!(
                reconstructed, r,
                "Decompose reconstruction failed for r={}: r1={}, r0={}", r, r1, r0
            );
        }
    }

    #[test]
    fn test_decompose_reconstruction_gamma2_65() {
        // For ML-DSA-65/87: γ₂ = 261888, α = 523776
        let gamma2 = (Q - 1) / 32;
        let alpha = 2 * gamma2;
        let test_vals = [0, 1, Q - 1, Q / 2, 261887, 261888, 523775, 523776];

        for &r in &test_vals {
            let (r1, r0) = decompose(r, gamma2);
            let reconstructed = ((r1 as i64) * (alpha as i64) + (r0 as i64))
                .rem_euclid(Q as i64) as u32;
            assert_eq!(
                reconstructed, r,
                "Decompose reconstruction failed for r={}: r1={}, r0={}", r, r1, r0
            );
        }
    }

    #[test]
    fn test_decompose_r0_range() {
        // r0 should be in (-α/2, α/2]
        let gamma2 = (Q - 1) / 88;
        let alpha = 2 * gamma2;
        let half = (alpha / 2) as i32;

        for r in (0..Q).step_by(10007) {
            let (_r1, r0) = decompose(r, gamma2);
            assert!(
                r0 > -half && r0 <= half,
                "r0={} out of range (-{}, {}] for r={}", r0, half, half, r
            );
        }
    }

    #[test]
    fn test_decompose_r1_range() {
        // r1 should be in [0, m) where m = (q-1)/α
        let gamma2 = (Q - 1) / 88;
        let m = (Q - 1) / (2 * gamma2); // = 44

        for r in (0..Q).step_by(10007) {
            let (r1, _r0) = decompose(r, gamma2);
            assert!(
                r1 < m,
                "r1={} out of range [0, {}) for r={}", r1, m, r
            );
        }
    }

    #[test]
    fn test_decompose_special_case() {
        // The special case: when r - r0 = q - 1, we set r1 = 0 and r0 -= 1.
        // This happens when r mod α gives a centered value such that
        // r - centered_r0 = q - 1.
        let gamma2 = (Q - 1) / 88;
        let (r1, _r0) = decompose(Q - 1, gamma2);
        // r1 should be 0 due to the special case
        assert_eq!(r1, 0, "Decompose special case: r1 should be 0 for r=q-1");
    }

    // ── HighBits / LowBits ──

    #[test]
    fn test_highbits_lowbits_consistency() {
        let gamma2 = (Q - 1) / 88;
        for r in (0..Q).step_by(10007) {
            let h = high_bits(r, gamma2);
            let l = low_bits(r, gamma2);
            let (r1, r0) = decompose(r, gamma2);
            assert_eq!(h, r1, "HighBits mismatch for r={}", r);
            assert_eq!(l, r0, "LowBits mismatch for r={}", r);
        }
    }

    // ── MakeHint / UseHint ──

    #[test]
    fn test_hint_no_change_means_false() {
        // When z = 0, adding it doesn't change HighBits, so hint = false.
        let gamma2 = (Q - 1) / 88;
        for r in (0..Q).step_by(50000) {
            assert_eq!(
                make_hint(0, r, gamma2), false,
                "MakeHint(0, {}) should be false", r
            );
        }
    }

    #[test]
    fn test_use_hint_no_hint_returns_r1() {
        // UseHint(false, r) should just return HighBits(r).
        let gamma2 = (Q - 1) / 88;
        for r in (0..Q).step_by(50000) {
            let r1 = high_bits(r, gamma2);
            assert_eq!(
                use_hint(false, r, gamma2), r1,
                "UseHint(false, {}) should equal HighBits({})", r, r
            );
        }
    }

    #[test]
    fn test_make_use_hint_roundtrip() {
        // The core property: for the verifier to work correctly,
        // UseHint(MakeHint(z, r), r+z) should equal HighBits(r).
        //
        // This is the fundamental hint property from the Dilithium paper.
        let gamma2 = (Q - 1) / 88;

        let test_cases: Vec<(i32, u32)> = vec![
            (0, 0),
            (0, Q - 1),
            (1, 100000),
            (-1, 100000),
            (1000, 95000),
            (-1000, 95000),
            (5000, Q / 2),
            (-5000, Q / 2),
        ];

        for (z, r) in test_cases {
            let rz = ((r as i64) + (z as i64)).rem_euclid(Q as i64) as u32;
            let hint = make_hint(z, r, gamma2);
            let recovered = use_hint(hint, rz, gamma2);
            let expected = high_bits(r, gamma2);
            assert_eq!(
                recovered, expected,
                "Hint roundtrip failed: z={}, r={}, hint={}, UseHint={}, HighBits(r)={}",
                z, r, hint, recovered, expected
            );
        }
    }

    #[test]
    fn test_make_use_hint_roundtrip_exhaustive_gamma2_44() {
        // Test MakeHint/UseHint roundtrip over many values for ML-DSA-44.
        let gamma2 = (Q - 1) / 88;

        for r in (0..Q).step_by(100003) {
            for &z in &[0i32, 1, -1, 100, -100, 1000, -1000] {
                let rz = ((r as i64) + (z as i64)).rem_euclid(Q as i64) as u32;
                let hint = make_hint(z, r, gamma2);
                let recovered = use_hint(hint, rz, gamma2);
                let expected = high_bits(r, gamma2);
                assert_eq!(
                    recovered, expected,
                    "Failed for r={}, z={}", r, z
                );
            }
        }
    }

    #[test]
    fn test_make_use_hint_roundtrip_exhaustive_gamma2_65() {
        // Same roundtrip test for ML-DSA-65/87 parameters.
        let gamma2 = (Q - 1) / 32;

        for r in (0..Q).step_by(100003) {
            for &z in &[0i32, 1, -1, 100, -100, 1000, -1000] {
                let rz = ((r as i64) + (z as i64)).rem_euclid(Q as i64) as u32;
                let hint = make_hint(z, r, gamma2);
                let recovered = use_hint(hint, rz, gamma2);
                let expected = high_bits(r, gamma2);
                assert_eq!(
                    recovered, expected,
                    "Failed for r={}, z={}", r, z
                );
            }
        }
    }

    // ── Cross-module: arithmetic + params consistency ──

    #[test]
    fn test_gamma2_values_produce_valid_decompositions() {
        // Both γ₂ values used in ML-DSA should divide (q-1)/2 evenly.
        let gamma2_44 = (Q - 1) / 88;
        let gamma2_65 = (Q - 1) / 32;

        // α = 2·γ₂ should divide q-1
        assert_eq!((Q - 1) % (2 * gamma2_44), 0);
        assert_eq!((Q - 1) % (2 * gamma2_65), 0);
    }
}