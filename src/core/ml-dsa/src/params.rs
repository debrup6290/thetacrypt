//! ML-DSA parameter sets per FIPS 204 (Module-Lattice Digital Signature Algorithm).
//!
//! Reference: FIPS 204, Table 1 and Section 4.
//!
//! All three security levels share the following constants:
//!   - q = 8380417        (the modulus, a prime with q ≡ 1 mod 512)
//!   - n = 256            (polynomial degree)
//!   - d = 13             (dropped bits in Power2Round)
//!
//! The three parameter sets differ in matrix dimensions, coefficient bounds,
//! and derived sizes for public keys, secret keys, and signatures.

// ═══════════════════════ Global Constants ═══════════════════════════

/// The modulus q. This is prime, and satisfies q ≡ 1 (mod 512), which is
/// needed for the NTT over polynomials of degree 256 (we need primitive
/// 512th roots of unity in Z_q).
///
/// Decomposition: q = 2^23 - 2^13 + 1 = 8380417
pub const Q: u32 = 8380417;

/// Polynomial degree: all polynomials live in R_q = Z_q[X]/(X^256 + 1).
pub const N: usize = 256;

/// Number of bits dropped by Power2Round.
/// The public vector t is split as t = t1·2^D + t0.
/// Only t1 is stored in the public key; t0 is kept in the secret key.
pub const D: u32 = 13;

// ═══════════════════════ Parameter Set ══════════════════════════════

/// Complete parameter set for one ML-DSA security level.
///
/// Notation follows FIPS 204:
///   - (k, l): matrix A has dimensions k × l
///   - η (eta): secret key coefficients are uniform in [-η, η]
///   - γ₁ (gamma1): masking vector y coefficients in [-(γ₁-1), γ₁-1]
///   - γ₂ (gamma2): low-order rounding range
///   - τ (tau): Hamming weight of challenge polynomial c
///   - β (beta) = τ·η: the norm bound used in rejection sampling
///   - ω (omega): max number of ones in the hint vector h
///   - λ: collision strength in bits; c̃ has λ/4 bytes
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub name: &'static str,

    // ── Matrix dimensions ──
    pub k: usize,
    pub l: usize,

    // ── Coefficient bounds ──
    pub eta: u32,
    pub gamma1: u32,
    pub gamma2: u32,
    pub tau: usize,
    pub beta: u32,
    pub omega: usize,

    // ── Hash / commitment sizes ──
    /// Byte length of the commitment hash c̃ = λ/4.
    pub c_tilde_bytes: usize,

    // ── Encoded byte sizes (FIPS 204 Table 2) ──
    pub pk_bytes: usize,
    pub sk_bytes: usize,
    pub sig_bytes: usize,
}

// ═══════════════════════ ML-DSA-44 ═════════════════════════════════
// NIST Security Level 2 (equivalent to AES-128)

pub const ML_DSA_44: Params = Params {
    name: "ML-DSA-44",

    k: 4,
    l: 4,

    eta: 2,
    gamma1: 1 << 17,       // 2^17 = 131_072
    gamma2: (Q - 1) / 88,  // = 95_232
    tau: 39,
    beta: 78,               // τ·η = 39·2
    omega: 80,

    c_tilde_bytes: 32,      // λ=128 → 128/4 = 32

    pk_bytes: 1312,
    sk_bytes: 2560,
    sig_bytes: 2420,
};

// ═══════════════════════ ML-DSA-65 ═════════════════════════════════
// NIST Security Level 3 (equivalent to AES-192)

pub const ML_DSA_65: Params = Params {
    name: "ML-DSA-65",

    k: 6,
    l: 5,

    eta: 4,
    gamma1: 1 << 19,       // 2^19 = 524_288
    gamma2: (Q - 1) / 32,  // = 261_888
    tau: 49,
    beta: 196,              // τ·η = 49·4
    omega: 55,

    c_tilde_bytes: 48,      // λ=192 → 192/4 = 48

    pk_bytes: 1952,
    sk_bytes: 4032,
    sig_bytes: 3309,
};

// ═══════════════════════ ML-DSA-87 ═════════════════════════════════
// NIST Security Level 5 (equivalent to AES-256)

pub const ML_DSA_87: Params = Params {
    name: "ML-DSA-87",

    k: 8,
    l: 7,

    eta: 2,
    gamma1: 1 << 19,       // 2^19 = 524_288
    gamma2: (Q - 1) / 32,  // = 261_888
    tau: 60,
    beta: 120,              // τ·η = 60·2
    omega: 75,

    c_tilde_bytes: 64,      // λ=256 → 256/4 = 64

    pk_bytes: 2592,
    sk_bytes: 4896,
    sig_bytes: 4627,
};

// ═══════════════════════ Helper Methods ═════════════════════════════

impl Params {
    /// Number of bits to encode one coefficient of the masking vector y.
    /// Coefficients live in [0, 2·γ₁], packed using γ₁+1 representable values.
    ///   γ₁ = 2^17 → 18 bits
    ///   γ₁ = 2^19 → 20 bits
    pub fn gamma1_bits(&self) -> u32 {
        if self.gamma1 == (1 << 17) { 18 } else { 20 }
    }

    /// The number of distinct high-bits values: m = (q-1) / (2·γ₂).
    ///   ML-DSA-44: m = 8380416 / 190464 = 44
    ///   ML-DSA-65/87: m = 8380416 / 523776 = 16
    pub fn high_bits_range(&self) -> u32 {
        (Q - 1) / (2 * self.gamma2)
    }

    /// Bits needed to encode one high-bits coefficient w1.
    ///
    /// w1 takes values in {0, 1, ..., m-1} where m = high_bits_range().
    /// We need ceil(log2(m)) bits, which equals bit_length(m - 1).
    ///   ML-DSA-44: m=44, m-1=43 → bit_length(43) = 6
    ///   ML-DSA-65/87: m=16, m-1=15 → bit_length(15) = 4
    pub fn w1_bits(&self) -> u32 {
        let m = self.high_bits_range();
        debug_assert!(m >= 2, "m must be at least 2");
        // bit_length(m-1) = number of bits to represent 0..m-1
        32 - (m - 1).leading_zeros()
    }

    /// Bytes per polynomial when encoding secret key coefficients in [-η, η].
    /// Shifted to [0, 2η], each coefficient uses ceil(log2(2η+1)) bits.
    ///   η=2: values in [0,4] → 3 bits → 256·3/8 = 96 bytes
    ///   η=4: values in [0,8] → 4 bits → 256·4/8 = 128 bytes
    pub fn eta_packed_bytes(&self) -> usize {
        match self.eta {
            2 => N * 3 / 8,  // 96
            4 => N * 4 / 8,  // 128
            _ => unreachable!("Unsupported eta value"),
        }
    }

    /// Bytes per polynomial for t0 (low bits from Power2Round).
    /// t0 ∈ [-(2^(d-1)-1), 2^(d-1)], packed with d=13 bits each.
    /// 256 · 13 / 8 = 416 bytes.
    pub fn t0_packed_bytes(&self) -> usize {
        N * D as usize / 8  // 416
    }

    /// Bytes per polynomial for t1 (high bits in public key).
    /// t1 needs 10 bits per coefficient: 256 · 10 / 8 = 320 bytes.
    pub fn t1_packed_bytes(&self) -> usize {
        N * 10 / 8  // 320
    }
}

// ═══════════════════════════ Tests ══════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── q properties ──

    #[test]
    fn test_q_decomposition() {
        assert_eq!(Q, (1u32 << 23) - (1u32 << 13) + 1);
        assert_eq!(Q, 8380417);
    }

    #[test]
    fn test_q_is_prime() {
        let q = Q as u64;
        let mut d = 2u64;
        while d * d <= q {
            assert!(q % d != 0, "q is divisible by {}", d);
            d += 1;
        }
    }

    #[test]
    fn test_q_ntt_compatible() {
        assert_eq!(Q % 512, 1, "q must satisfy q ≡ 1 (mod 512)");
    }

    // ── Derived parameter consistency ──

    #[test]
    fn test_beta_equals_tau_times_eta() {
        for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
            assert_eq!(
                p.beta,
                p.tau as u32 * p.eta,
                "{}: β should equal τ·η",
                p.name
            );
        }
    }

    #[test]
    fn test_gamma2_exact_values() {
        assert_eq!(ML_DSA_44.gamma2, 95232, "ML-DSA-44 γ₂");
        assert_eq!(ML_DSA_65.gamma2, 261888, "ML-DSA-65 γ₂");
        assert_eq!(ML_DSA_87.gamma2, 261888, "ML-DSA-87 γ₂");
    }

    #[test]
    fn test_high_bits_range_values() {
        assert_eq!(ML_DSA_44.high_bits_range(), 44);
        assert_eq!(ML_DSA_65.high_bits_range(), 16);
        assert_eq!(ML_DSA_87.high_bits_range(), 16);
    }

    // ── Encoded sizes ──

    #[test]
    fn test_pk_size_formula() {
        for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
            let expected = 32 + p.k * p.t1_packed_bytes();
            assert_eq!(
                p.pk_bytes, expected,
                "{}: pk size mismatch", p.name
            );
        }
    }

    #[test]
    fn test_sk_size_formula() {
        for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
            let expected = 32 + 32 + 64
                + p.l * p.eta_packed_bytes()
                + p.k * p.eta_packed_bytes()
                + p.k * p.t0_packed_bytes();
            assert_eq!(
                p.sk_bytes, expected,
                "{}: sk size mismatch (formula={}, table={})",
                p.name, expected, p.sk_bytes
            );
        }
    }

    #[test]
    fn test_sig_size_formula() {
        for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
            let z_bytes_per_poly = N * p.gamma1_bits() as usize / 8;
            let expected = p.c_tilde_bytes + p.l * z_bytes_per_poly + p.omega + p.k;
            assert_eq!(
                p.sig_bytes, expected,
                "{}: sig size mismatch (formula={}, table={})",
                p.name, expected, p.sig_bytes
            );
        }
    }

    // ── Helper methods ──

    #[test]
    fn test_gamma1_bits() {
        assert_eq!(ML_DSA_44.gamma1_bits(), 18);
        assert_eq!(ML_DSA_65.gamma1_bits(), 20);
        assert_eq!(ML_DSA_87.gamma1_bits(), 20);
    }

    #[test]
    fn test_w1_bits() {
        // ML-DSA-44: m=44, values 0..43, need 6 bits (2^5=32 < 44 ≤ 64=2^6)
        assert_eq!(ML_DSA_44.w1_bits(), 6);
        // ML-DSA-65: m=16, values 0..15, need 4 bits (2^4=16, 15 fits in 4 bits)
        assert_eq!(ML_DSA_65.w1_bits(), 4);
        // ML-DSA-87: same as 65
        assert_eq!(ML_DSA_87.w1_bits(), 4);
    }

    #[test]
    fn test_eta_packed_bytes() {
        assert_eq!(ML_DSA_44.eta_packed_bytes(), 96);
        assert_eq!(ML_DSA_65.eta_packed_bytes(), 128);
        assert_eq!(ML_DSA_87.eta_packed_bytes(), 96);
    }

    #[test]
    fn test_t0_packed_bytes() {
        assert_eq!(ML_DSA_44.t0_packed_bytes(), 416);
        assert_eq!(ML_DSA_65.t0_packed_bytes(), 416);
        assert_eq!(ML_DSA_87.t0_packed_bytes(), 416);
    }

    #[test]
    fn test_t1_packed_bytes() {
        assert_eq!(ML_DSA_44.t1_packed_bytes(), 320);
        assert_eq!(ML_DSA_65.t1_packed_bytes(), 320);
        assert_eq!(ML_DSA_87.t1_packed_bytes(), 320);
    }
}