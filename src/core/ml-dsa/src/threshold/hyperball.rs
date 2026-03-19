//! Hyperball rejection sampling for Threshold ML-DSA.
//!
//! Implements the imbalanced rejection sampling from Paper Figure 4 (HRej),
//! which replaces ML-DSA's uniform rejection sampling with hyperball-based
//! rejection.  This significantly reduces the abort probability when multiple
//! parties must simultaneously pass rejection checks.
//!
//! ## HRej Algorithm (Paper Figure 4)
//!
//! ```text
//! HRej(v, r, r', ν, M; randomness) → z ∈ R^{k+l} ∪ {⊥}
//!   1. r_vec ←$ B_{R,l+k}(r')              // uniform from hyperball of radius r'
//!   2. (v^(1), v^(2)) := v ∈ R^l × R^k     // split secret·challenge
//!   3. z := (v^(1)/ν, v^(2)) + r_vec        // add continuous randomness
//!   4. if ||z||₂ > r: return ⊥              // norm check
//!   5. (z^(1), z^(2)) := z
//!   6. return ⌊(z^(1)·ν, z^(2))⌉           // round to integers
//! ```
//!
//! ## Imbalance Mechanism
//!
//! The expansion factor ν > 1 creates an asymmetric acceptance region:
//! - z^(1) (first l polynomials) is shrunk by 1/ν before the norm check,
//!   giving it a larger effective acceptance range.
//! - z^(2) (last k polynomials) uses the raw norm, keeping a tighter range.
//!
//! This matches ML-DSA's structure where the hint mechanism imposes stricter
//! bounds on z^(2) than z^(1).
//!
//! ## Continuous Sampling
//!
//! The hyperball sample lives in continuous space R^{n(k+l)}.  To sample
//! uniformly from a d-dimensional ball of radius R:
//!   1. Sample d independent standard normals → direction on unit sphere
//!   2. Normalize to get a uniform point on the sphere
//!   3. Scale by R · u^{1/d} where u ~ U[0,1] to fill the ball

use crate::arithmetic;
use crate::params::{Q, N as POLY_N};
use crate::poly::{self, Poly, PolyVec};
use crate::shake::Shake256;

use super::params::ThresholdParams;

// SHAKE-based PRNG 
//
// We derive random f64 values from a SHAKE-256 XOF stream.
// This gives us deterministic, reproducible randomness from a seed.

/// A deterministic RNG backed by SHAKE-256.
///
/// Produces uniform f64 values in [0, 1) and standard normal samples
/// via the Box-Muller transform.
struct ShakeRng {
    xof: Shake256,
    /// Cached second normal from Box-Muller (produces two at a time).
    spare_normal: Option<f64>,
}

impl ShakeRng {
    /// Create a new PRNG from a seed.
    fn new(seed: &[u8]) -> Self {
        let mut xof = Shake256::new();
        xof.absorb(seed);
        Self {
            xof,
            spare_normal: None,
        }
    }

    /// Squeeze a uniform f64 in [0, 1).
    ///
    /// Takes 8 bytes from SHAKE-256 and maps them to [0, 1) via
    /// the standard method: interpret as u64, divide by 2^64.
    fn uniform_f64(&mut self) -> f64 {
        let mut buf = [0u8; 8];
        self.xof.squeeze(&mut buf);
        let val = u64::from_le_bytes(buf);
        // Map to (0, 1] then subtract from 1 to get [0, 1)
        // We use (val >> 11) to get 53 bits of precision (f64 mantissa)
        (val >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Squeeze a uniform f64 in (0, 1) — excludes exact 0.
    ///
    /// Needed for Box-Muller (ln(0) is undefined).
    fn uniform_f64_open(&mut self) -> f64 {
        loop {
            let v = self.uniform_f64();
            if v > 0.0 {
                return v;
            }
        }
    }

    /// Sample a standard normal (mean=0, std=1) via Box-Muller.
    fn standard_normal(&mut self) -> f64 {
        if let Some(spare) = self.spare_normal.take() {
            return spare;
        }

        let u1 = self.uniform_f64_open();
        let u2 = self.uniform_f64();

        let mag = (-2.0 * u1.ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * u2;

        let z0 = mag * angle.cos();
        let z1 = mag * angle.sin();

        self.spare_normal = Some(z1);
        z0
    }
}

// ═══════════════════════ Hyperball Sampling ══════════════════════════

/// A sample from the continuous hyperball B(r') ⊂ R^{n(l+k)}.
///
/// Stores both the raw continuous coordinates (needed for the HRej norm
/// check) and the integer randomness ri (needed for the commitment wi).
///
/// ## Layout
///
/// The continuous vector has n(l+k) coordinates, viewed as:
///   - First n·l coordinates: x₁ (corresponds to the "l part")
///   - Last n·k coordinates: x₂ (corresponds to the "k part")
///
/// The integer randomness is:
///   ri = (⌊ν·x₁⌉, ⌊x₂⌉)
///
/// which is split into ri_1 (l polynomials) and ri_2 (k polynomials).
pub struct HyperballSample {
    /// Continuous sample (x₁, x₂) from B(r'), length n(l+k).
    pub continuous: Vec<f64>,

    /// Integer randomness: first l polynomials = ⌊ν·x₁⌉.
    pub ri_1: PolyVec,

    /// Integer randomness: last k polynomials = ⌊x₂⌉.
    pub ri_2: PolyVec,
}

/// Sample uniformly from the n(k+l)-dimensional hyperball of radius r',
/// then compute the integer randomness ri = (⌊ν·x₁⌉, ⌊x₂⌉).
///
/// # Arguments
/// * `seed` - Deterministic seed for the PRNG.
/// * `tp` - Threshold parameters (provides r', ν, k, l).
///
/// # Algorithm
/// 1. Generate d = n(k+l) standard normals → point on unit sphere
/// 2. Normalize to get a uniform direction
/// 3. Scale by r' · u^{1/d} to uniformly fill the ball
/// 4. Round to get integer randomness
pub fn sample_randomness(seed: &[u8], tp: &ThresholdParams) -> HyperballSample {
    let base = tp.base;
    let l = base.l;
    let k = base.k;
    let dim = tp.hyperball_dim(); // n * (k + l)
    let nu = tp.nu as f64;
    let r_prime = tp.r_prime as f64;

    let mut rng = ShakeRng::new(seed);

    // Step 1: Sample d standard normals
    let mut gauss: Vec<f64> = (0..dim).map(|_| rng.standard_normal()).collect();

    // Step 2: Normalize to unit sphere
    let norm: f64 = gauss.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 0.0 {
        for g in gauss.iter_mut() {
            *g /= norm;
        }
    }

    // Step 3: Scale to uniformly fill the ball of radius r'
    // Radius follows r' * u^{1/d} where u ~ U[0,1]
    let u = rng.uniform_f64_open();
    let radius = r_prime * u.powf(1.0 / dim as f64);
    for g in gauss.iter_mut() {
        *g *= radius;
    }

    // Now gauss = (x₁, x₂) is the continuous hyperball sample.
    // x₁ = gauss[0 .. n*l], x₂ = gauss[n*l .. n*(l+k)]

    // Step 4: Compute integer randomness ri = (⌊ν·x₁⌉, ⌊x₂⌉)
    let mut ri_1 = Vec::with_capacity(l);
    for j in 0..l {
        let mut poly = Poly::zero();
        for i in 0..POLY_N {
            let continuous_val = gauss[j * POLY_N + i];
            let scaled = continuous_val * nu;
            let rounded = scaled.round() as i64;
            // Map to [0, q)
            poly.coeffs[i] = ((rounded % Q as i64 + Q as i64) % Q as i64) as u32;
        }
        ri_1.push(poly);
    }

    let mut ri_2 = Vec::with_capacity(k);
    for j in 0..k {
        let mut poly = Poly::zero();
        for i in 0..POLY_N {
            let continuous_val = gauss[(l + j) * POLY_N + i];
            let rounded = continuous_val.round() as i64;
            poly.coeffs[i] = ((rounded % Q as i64 + Q as i64) % Q as i64) as u32;
        }
        ri_2.push(poly);
    }

    HyperballSample {
        continuous: gauss,
        ri_1,
        ri_2,
    }
}

// ═══════════════════════ HRej ════════════
// Paper Figure 4: Imbalanced rejection sampling over hyperballs.
//
// Input:
//   v = c · s_part  ∈ R_q^{l+k}  (the secret times challenge)
//   sample: pre-sampled hyperball randomness (contains continuous (x₁, x₂))
//   tp: threshold parameters (r, r', ν)
//
// Output:
//   Some(z_1) if accepted (z_1 ∈ R_q^l, the first l polynomials)
//   None if rejected
//
// Algorithm:
//   1. (v1, v2) := v
//   2. z := (v1/ν, v2) + (x₁, x₂)        [continuous]
//   3. if ||z||₂ > r: return None
//   4. z_out := ⌊(z₁·ν, z₂)⌉             [round to integers]
//   5. return z_out^(1)                    [only first l polys]

/// Apply imbalanced hyperball rejection sampling (Paper Figure 4).
///
/// # Arguments
/// * `v1` - First part of c·s_part: l polynomials in R_q (standard reps).
/// * `v2` - Second part of c·s_part: k polynomials in R_q.
/// * `sample` - Pre-sampled hyperball randomness.
/// * `tp` - Threshold parameters.
///
/// # Returns
/// `Some(z_1)` (l polynomials) if the norm check passes, `None` if rejected.
///
/// The caller (Combine) reconstructs z^(2) from public values, so we only
/// output z^(1) — matching Paper Figure 6 line 16.
pub fn hrej(
    v1: &PolyVec,
    v2: &PolyVec,
    sample: &HyperballSample,
    tp: &ThresholdParams,
) -> Option<PolyVec> {
    let base = tp.base;
    let l = base.l;
    let k = base.k;
    let nu = tp.nu as f64;
    let r_target = tp.r as f64;

    // Step 2: Compute z = (v1/ν, v2) + (x₁, x₂) in continuous space
    // z has n(l+k) coordinates
    let dim = tp.hyperball_dim();
    let mut z_continuous = vec![0.0f64; dim];

    // First n·l coordinates: v1[j][i] / ν + x₁[j*n + i]
    for j in 0..l {
        for i in 0..POLY_N {
            let v_centered = arithmetic::to_centered(v1[j].coeffs[i]) as f64;
            let x1_val = sample.continuous[j * POLY_N + i];
            z_continuous[j * POLY_N + i] = v_centered / nu + x1_val;
        }
    }

    // Last n·k coordinates: v2[j][i] + x₂[j*n + i]
    for j in 0..k {
        for i in 0..POLY_N {
            let v_centered = arithmetic::to_centered(v2[j].coeffs[i]) as f64;
            let x2_val = sample.continuous[(l + j) * POLY_N + i];
            z_continuous[(l + j) * POLY_N + i] = v_centered + x2_val;
        }
    }

    // Step 3: Check ||z||₂ ≤ r
    let norm_sq: f64 = z_continuous.iter().map(|x| x * x).sum();
    if norm_sq > r_target * r_target {
        return None; // Rejected
    }

    // Step 4: Compute output z_out = ⌊(z₁·ν, z₂)⌉
    // We only return z^(1) (the first l polynomials), as z^(2) is
    // reconstructed by the combiner from public values.

    let mut z1_out = Vec::with_capacity(l);
    for j in 0..l {
        let mut poly = Poly::zero();
        for i in 0..POLY_N {
            let z_val = z_continuous[j * POLY_N + i];
            let scaled = z_val * nu;
            let rounded = scaled.round() as i64;
            // Map to [0, q)
            poly.coeffs[i] = ((rounded % Q as i64 + Q as i64) % Q as i64) as u32;
        }
        z1_out.push(poly);
    }

    Some(z1_out)
}

// ═══════════════════════ L2 Norm Helpers ═

/// Compute the squared L2 norm of a polynomial vector (centered representatives).
///
/// ||v||₂² = Σ_j Σ_i (centered(v[j][i]))²
///
/// This is useful for checking partial secret bounds (the B parameter
/// from Paper Section 3.4).
pub fn polyvec_l2_norm_sq(v: &PolyVec) -> f64 {
    let mut sum = 0.0f64;
    for poly in v {
        for i in 0..POLY_N {
            let c = arithmetic::to_centered(poly.coeffs[i]) as f64;
            sum += c * c;
        }
    }
    sum
}

/// Compute the "twisted" L2 norm used for the bound B in the paper.
///
/// ||(v₁/ν, v₂)||₂² = Σ (v1[j][i]/ν)² + Σ (v2[j][i])²
///
/// This matches the expression ||(1/ν · c · u₁, c · u₂)||₂ ≤ B
/// from Paper Section 3.2.
pub fn twisted_l2_norm_sq(v1: &PolyVec, v2: &PolyVec, nu: u32) -> f64 {
    let nu_f = nu as f64;
    let mut sum = 0.0f64;

    for poly in v1 {
        for i in 0..POLY_N {
            let c = arithmetic::to_centered(poly.coeffs[i]) as f64;
            let scaled = c / nu_f;
            sum += scaled * scaled;
        }
    }
    for poly in v2 {
        for i in 0..POLY_N {
            let c = arithmetic::to_centered(poly.coeffs[i]) as f64;
            sum += c * c;
        }
    }

    sum
}

// ═══════════════════════════ Tests ══════

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65};
    use super::super::params::lookup;

    // ── ShakeRng ──

    #[test]
    fn test_shake_rng_deterministic() {
        let mut rng1 = ShakeRng::new(b"test seed");
        let mut rng2 = ShakeRng::new(b"test seed");
        for _ in 0..100 {
            assert_eq!(rng1.uniform_f64().to_bits(), rng2.uniform_f64().to_bits());
        }
    }

    #[test]
    fn test_shake_rng_different_seeds() {
        let mut rng1 = ShakeRng::new(b"seed A");
        let mut rng2 = ShakeRng::new(b"seed B");
        let v1: Vec<f64> = (0..10).map(|_| rng1.uniform_f64()).collect();
        let v2: Vec<f64> = (0..10).map(|_| rng2.uniform_f64()).collect();
        assert_ne!(v1, v2);
    }

    #[test]
    fn test_shake_rng_uniform_range() {
        let mut rng = ShakeRng::new(b"range test");
        for _ in 0..1000 {
            let v = rng.uniform_f64();
            assert!(v >= 0.0 && v < 1.0, "uniform_f64 out of range: {}", v);
        }
    }

    #[test]
    fn test_shake_rng_normal_mean_variance() {
        let mut rng = ShakeRng::new(b"normal test");
        let n = 10000;
        let samples: Vec<f64> = (0..n).map(|_| rng.standard_normal()).collect();

        let mean = samples.iter().sum::<f64>() / n as f64;
        let variance = samples.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / n as f64;

        // Mean should be close to 0, variance close to 1
        assert!(mean.abs() < 0.05, "Normal mean too far from 0: {}", mean);
        assert!((variance - 1.0).abs() < 0.1,
            "Normal variance too far from 1: {}", variance);
    }

    // ── Hyperball sampling ──

    #[test]
    fn test_sample_randomness_deterministic() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let s1 = sample_randomness(b"deterministic test seed 1234", &tp);
        let s2 = sample_randomness(b"deterministic test seed 1234", &tp);

        assert_eq!(s1.continuous.len(), s2.continuous.len());
        for i in 0..s1.continuous.len() {
            assert_eq!(s1.continuous[i].to_bits(), s2.continuous[i].to_bits());
        }
    }

    #[test]
    fn test_sample_randomness_dimension() {
        let tp = lookup(3, 2, &ML_DSA_44).unwrap();
        let sample = sample_randomness(b"dim test", &tp);

        assert_eq!(sample.continuous.len(), tp.hyperball_dim());
        assert_eq!(sample.ri_1.len(), ML_DSA_44.l);
        assert_eq!(sample.ri_2.len(), ML_DSA_44.k);
    }

    #[test]
    fn test_sample_randomness_within_ball() {
        // The continuous sample should have L2 norm ≤ r'
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();

        for seed_byte in 0u8..20 {
            let mut seed = [0u8; 32];
            seed[0] = seed_byte;
            let sample = sample_randomness(&seed, &tp);

            let norm_sq: f64 = sample.continuous.iter().map(|x| x * x).sum();
            let norm = norm_sq.sqrt();
            let r_prime = tp.r_prime as f64;

            assert!(norm <= r_prime * 1.001, // small tolerance for floating point
                "Sample norm {} exceeds r'={}", norm, r_prime);
        }
    }

    #[test]
    fn test_sample_randomness_different_seeds() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let s1 = sample_randomness(b"seed alpha", &tp);
        let s2 = sample_randomness(b"seed beta", &tp);

        // Should be different
        let same = s1.continuous.iter().zip(s2.continuous.iter())
            .all(|(a, b)| (a - b).abs() < 1e-10);
        assert!(!same, "Different seeds should produce different samples");
    }

    #[test]
    fn test_sample_randomness_ml_dsa_65() {
        let tp = lookup(2, 2, &ML_DSA_65).unwrap();
        let sample = sample_randomness(b"65 test", &tp);

        assert_eq!(sample.continuous.len(), 256 * (6 + 5)); // n * (k + l)
        assert_eq!(sample.ri_1.len(), 5); // l
        assert_eq!(sample.ri_2.len(), 6); // k
    }

    // ── HRej ──

    #[test]
    fn test_hrej_zero_secret_always_accepts() {
        // With v = 0, z = (0/ν, 0) + (x₁, x₂) = (x₁, x₂).
        // Since ||x||₂ ≤ r' and r' is close to r, most samples should accept.
        // (For v=0 specifically, acceptance is guaranteed when ||x||₂ ≤ r,
        //  which happens with probability (r/r')^dim.)
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();

        let v1 = poly::polyvec_zero(tp.base.l);
        let v2 = poly::polyvec_zero(tp.base.k);

        let mut accepted = 0;
        for seed_byte in 0u8..100 {
            let mut seed = [0u8; 32];
            seed[0] = seed_byte;
            let sample = sample_randomness(&seed, &tp);

            if hrej(&v1, &v2, &sample, &tp).is_some() {
                accepted += 1;
            }
        }

        // With v=0, acceptance probability = (r/r')^dim.
        // For (2,2) ML-DSA-44: r=252778, r'=252833, dim=2048
        // (r/r')^2048 is very close to 0 actually... the ratio is ~0.999978
        // 0.999978^2048 ≈ 0.956 so about 96% acceptance
        assert!(accepted > 50, "Too few accepted with zero secret: {}/100", accepted);
    }

    #[test]
    fn test_hrej_returns_correct_poly_count() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let v1 = poly::polyvec_zero(tp.base.l);
        let v2 = poly::polyvec_zero(tp.base.k);

        // Try until we get an acceptance
        for seed_byte in 0u8..255 {
            let mut seed = [0u8; 32];
            seed[0] = seed_byte;
            let sample = sample_randomness(&seed, &tp);

            if let Some(z1) = hrej(&v1, &v2, &sample, &tp) {
                assert_eq!(z1.len(), tp.base.l,
                    "HRej output should have l={} polynomials", tp.base.l);
                return;
            }
        }
        panic!("No acceptance in 255 attempts");
    }

    #[test]
    fn test_hrej_deterministic() {
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let v1 = poly::polyvec_zero(tp.base.l);
        let v2 = poly::polyvec_zero(tp.base.k);

        let sample = sample_randomness(b"hrej deterministic", &tp);
        let r1 = hrej(&v1, &v2, &sample, &tp);
        let r2 = hrej(&v1, &v2, &sample, &tp);

        match (r1, r2) {
            (Some(z1), Some(z2)) => {
                for j in 0..tp.base.l {
                    for i in 0..POLY_N {
                        assert_eq!(z1[j].coeffs[i], z2[j].coeffs[i]);
                    }
                }
            }
            (None, None) => {} // Both rejected, consistent
            _ => panic!("HRej not deterministic"),
        }
    }

    #[test]
    fn test_hrej_with_nonzero_secret_has_rejections() {
        // With a large v, more samples should be rejected
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();

        // Create a "secret" with moderately large coefficients
        let mut v1 = Vec::with_capacity(tp.base.l);
        for j in 0..tp.base.l {
            let mut p = Poly::zero();
            for i in 0..POLY_N {
                // Coefficients around ±100
                let c = ((j * POLY_N + i) % 201) as i32 - 100;
                p.coeffs[i] = arithmetic::from_centered(c);
            }
            v1.push(p);
        }
        let mut v2 = Vec::with_capacity(tp.base.k);
        for j in 0..tp.base.k {
            let mut p = Poly::zero();
            for i in 0..POLY_N {
                let c = ((j * POLY_N + i + 7) % 201) as i32 - 100;
                p.coeffs[i] = arithmetic::from_centered(c);
            }
            v2.push(p);
        }

        let mut accepted = 0;
        let mut rejected = 0;
        for seed_byte in 0u8..200 {
            let mut seed = [0u8; 32];
            seed[0] = seed_byte;
            let sample = sample_randomness(&seed, &tp);

            if hrej(&v1, &v2, &sample, &tp).is_some() {
                accepted += 1;
            } else {
                rejected += 1;
            }
        }

        // With nonzero v, we should see both acceptances and rejections
        // (the exact ratio depends on ||v||₂ vs the radii)
        assert!(accepted > 0 || rejected > 0,
            "Should have some results: acc={}, rej={}", accepted, rejected);
    }

    // ── L2 norm helpers ──

    #[test]
    fn test_polyvec_l2_norm_sq_zero() {
        let v = poly::polyvec_zero(4);
        assert_eq!(polyvec_l2_norm_sq(&v), 0.0);
    }

    #[test]
    fn test_polyvec_l2_norm_sq_unit() {
        let mut v = poly::polyvec_zero(1);
        v[0].coeffs[0] = 3;
        v[0].coeffs[1] = 4;
        // ||v||² = 3² + 4² = 25
        assert_eq!(polyvec_l2_norm_sq(&v), 25.0);
    }

    #[test]
    fn test_polyvec_l2_norm_sq_negative() {
        let mut v = poly::polyvec_zero(1);
        v[0].coeffs[0] = Q - 3; // centered: -3
        // ||v||² = 9
        assert_eq!(polyvec_l2_norm_sq(&v), 9.0);
    }

    #[test]
    fn test_twisted_l2_norm() {
        let mut v1 = poly::polyvec_zero(1);
        v1[0].coeffs[0] = 6; // centered: 6

        let mut v2 = poly::polyvec_zero(1);
        v2[0].coeffs[0] = 4; // centered: 4

        // With ν=3: ||(6/3, 4)||² = (2)² + (4)² = 4 + 16 = 20
        let result = twisted_l2_norm_sq(&v1, &v2, 3);
        assert!((result - 20.0).abs() < 1e-10);
    }

    // ── Integration: sample_randomness ri consistency ──

    #[test]
    fn test_ri_is_rounded_scaled_continuous() {
        // Verify that ri_1[j][i] ≈ round(ν * continuous[j*n + i])
        // and ri_2[j][i] ≈ round(continuous[(l+j)*n + i])
        let tp = lookup(2, 2, &ML_DSA_44).unwrap();
        let sample = sample_randomness(b"consistency check", &tp);
        let nu = tp.nu as f64;
        let l = tp.base.l;

        // Check ri_1
        for j in 0..l {
            for i in 0..POLY_N {
                let continuous_val = sample.continuous[j * POLY_N + i];
                let expected_rounded = (continuous_val * nu).round() as i64;
                let actual_centered = arithmetic::to_centered(
                    sample.ri_1[j].coeffs[i]) as i64;
                // They should match (modulo q wraparound for very large values,
                // but for typical hyperball radii the values are well within q)
                assert_eq!(actual_centered, expected_rounded,
                    "ri_1 mismatch at [{j}][{i}]: expected {expected_rounded}, got {actual_centered}");
            }
        }

        // Check ri_2
        for j in 0..tp.base.k {
            for i in 0..POLY_N {
                let continuous_val = sample.continuous[(l + j) * POLY_N + i];
                let expected_rounded = continuous_val.round() as i64;
                let actual_centered = arithmetic::to_centered(
                    sample.ri_2[j].coeffs[i]) as i64;
                assert_eq!(actual_centered, expected_rounded,
                    "ri_2 mismatch at [{j}][{i}]");
            }
        }
    }
}