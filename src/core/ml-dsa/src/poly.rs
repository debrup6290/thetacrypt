//! Polynomials in R_q = Z_q[X]/(X^256 + 1) and the Number Theoretic Transform.
//!
//! This module implements:
//!   - The `Poly` type: 256 coefficients in [0, q)
//!   - Forward NTT (FIPS 204, Algorithm 41)
//!   - Inverse NTT (FIPS 204, Algorithm 42)
//!   - Pointwise (coefficient-wise) operations in NTT domain
//!   - Vector and matrix operations over polynomials
//!
//! ## NTT Background
//!
//! The NTT is the key to efficient polynomial multiplication.  In the
//! ring R_q = Z_q[X]/(X^256+1), multiplying two polynomials naively
//! costs O(n²) = O(65536) multiplications.  Via the NTT, it costs O(n log n).
//!
//! The NTT exploits the factorization of X^256+1 over Z_q.  Since
//! q ≡ 1 (mod 512), there exists a primitive 512th root of unity
//! ζ = 1753, meaning ζ^256 = -1 mod q.  This lets X^256+1 split into
//! 256 linear factors, and the NTT evaluates the polynomial at all
//! 256 roots.
//!
//! ## Zeta Table
//!
//! The NTT butterfly needs specific powers of ζ accessed in bit-reversed
//! order.  For index k (1 ≤ k ≤ 255), the twiddle factor is:
//!
//!   zetas[k] = ζ^{BitRev8(k)} mod q
//!
//! where BitRev8 reverses the 8-bit representation of k.

use crate::arithmetic;
use crate::params::{N, Q};

//  Poly Type 

/// A polynomial in R_q = Z_q[X]/(X^256 + 1).
///
/// Coefficients are stored as `[u32; 256]` with each value in [0, q).
/// The polynomial a₀ + a₁X + a₂X² + ... + a₂₅₅X^255 is stored as
/// `coeffs[i] = aᵢ`.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Poly {
    #[serde(with = "serde_arrays")]
    pub coeffs: [u32; N],
}

impl Default for Poly {
    fn default() -> Self {
        Self { coeffs: [0u32; N] }
    }
}

impl core::fmt::Debug for Poly {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Poly([{}, {}, {}, ... {}, {}])",
            self.coeffs[0], self.coeffs[1], self.coeffs[2],
            self.coeffs[N - 2], self.coeffs[N - 1])
    }
}

impl Poly {
    /// The zero polynomial.
    pub fn zero() -> Self {
        Self::default()
    }

    /// Coefficient-wise addition: (self + other) mod q.
    pub fn add(&self, other: &Poly) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = arithmetic::add(self.coeffs[i], other.coeffs[i]);
        }
        result
    }

    /// Coefficient-wise subtraction: (self - other) mod q.
    pub fn sub(&self, other: &Poly) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = arithmetic::sub(self.coeffs[i], other.coeffs[i]);
        }
        result
    }

    /// Pointwise (coefficient-wise) multiplication: used in NTT domain.
    ///
    /// In the NTT domain, polynomial multiplication becomes pointwise:
    ///   NTT(a · b) = NTT(a) ⊙ NTT(b)
    pub fn pointwise_mul(&self, other: &Poly) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = arithmetic::mul(self.coeffs[i], other.coeffs[i]);
        }
        result
    }

    /// Multiply all coefficients by a scalar.
    pub fn scalar_mul(&self, s: u32) -> Poly {
        let mut result = Poly::zero();
        for i in 0..N {
            result.coeffs[i] = arithmetic::mul(self.coeffs[i], s);
        }
        result
    }

    /// Check if all coefficients (as centered representatives) have
    /// infinity norm strictly less than `bound`.
    ///
    /// Returns true if |c_i|_centered < bound for all i.
    /// Used in the rejection check during signing.
    pub fn check_norm(&self, bound: u32) -> bool {
        for i in 0..N {
            let c = arithmetic::to_centered(self.coeffs[i]);
            if c.unsigned_abs() >= bound {
                return false;
            }
        }
        true
    }

    /// Reduce all coefficients mod q.
    /// Useful after accumulating many additions without intermediate reduction.
    pub fn reduce(&mut self) {
        for i in 0..N {
            self.coeffs[i] = arithmetic::reduce(self.coeffs[i] as u64);
        }
    }
}

//  Zeta Table 

/// The primitive 512th root of unity in Z_q.
///
/// ζ = 1753 satisfies ζ^256 = -1 (mod q) and ζ^512 = 1 (mod q).
const ZETA: u32 = 1753;

/// Reverse the lowest 8 bits of `x`.
///
/// Used to compute the bit-reversed index for the NTT twiddle factors.
/// Example: bit_rev8(1) = 128, bit_rev8(2) = 64, bit_rev8(3) = 192.
fn bit_rev8(x: u32) -> u32 {
    let mut result = 0u32;
    let mut v = x;
    for _ in 0..8 {
        result = (result << 1) | (v & 1);
        v >>= 1;
    }
    result
}

/// Compute the full zeta table: zetas[k] = ζ^{BitRev8(k)} mod q.
///
/// The NTT uses zetas[1..=255]; zetas[0] = ζ^0 = 1 (not directly used).
///
/// This is computed from scratch to avoid hardcoding errors.
/// In a production implementation, you'd precompute and embed this table.
fn compute_zetas() -> [u32; N] {
    let mut zetas = [0u32; N];
    for k in 0..N {
        let exp = bit_rev8(k as u32);
        zetas[k] = arithmetic::pow_mod(ZETA, exp);
    }
    zetas
}

//  Forward NTT 
// FIPS 204, Algorithm 41: NTT(f)
//
// Cooley-Tukey butterfly.  Transforms a polynomial from the "normal"
// domain to the NTT domain.  After transformation, polynomial
// multiplication becomes pointwise.
//
// Pseudocode from FIPS 204:
//   k ← 0
//   for len in [128, 64, 32, 16, 8, 4, 2, 1]:
//     for start in [0, 2·len, 4·len, ...]:
//       k ← k + 1
//       z ← zetas[k]
//       for j in [start, start+1, ..., start+len-1]:
//         t ← z · f̂[j+len]
//         f̂[j+len] ← f̂[j] - t
//         f̂[j]     ← f̂[j] + t
//   return f̂

/// Forward NTT: convert polynomial from normal domain to NTT domain.
pub fn ntt(a: &Poly) -> Poly {
    let zetas = compute_zetas();
    let mut f = a.clone();

    let mut k: usize = 0;
    let mut len: usize = 128;
    while len >= 1 {
        let mut start: usize = 0;
        while start < N {
            k += 1;
            let z = zetas[k];
            for j in start..(start + len) {
                let t = arithmetic::mul(z, f.coeffs[j + len]);
                f.coeffs[j + len] = arithmetic::sub(f.coeffs[j], t);
                f.coeffs[j] = arithmetic::add(f.coeffs[j], t);
            }
            start += 2 * len;
        }
        len >>= 1;
    }

    f
}

//  Inverse NTT 
// FIPS 204, Algorithm 42: NTT⁻¹(f̂)
//
// Gentleman-Sande butterfly.  Transforms from NTT domain back to
// normal domain.
//
// Pseudocode from FIPS 204:
//   k ← 256
//   for len in [1, 2, 4, 8, 16, 32, 64, 128]:
//     for start in [0, 2·len, 4·len, ...]:
//       k ← k - 1
//       z ← -zetas[k]     (i.e., q - zetas[k])
//       for j in [start, start+1, ..., start+len-1]:
//         t ← f[j]
//         f[j]     ← t + f[j+len]
//         f[j+len] ← z · (t - f[j+len])
//   f ← f · n⁻¹           (multiply every coefficient by 256⁻¹ mod q)
//   return f

/// Inverse NTT: convert from NTT domain back to normal domain.
pub fn inv_ntt(a: &Poly) -> Poly {
    let zetas = compute_zetas();
    let mut f = a.clone();

    let mut k: usize = N; // 256
    let mut len: usize = 1;
    while len < N {
        let mut start: usize = 0;
        while start < N {
            k -= 1;
            let z = Q - zetas[k]; // -zetas[k] mod q
            for j in start..(start + len) {
                let t = f.coeffs[j];
                f.coeffs[j] = arithmetic::add(t, f.coeffs[j + len]);
                f.coeffs[j + len] = arithmetic::mul(z, arithmetic::sub(t, f.coeffs[j + len]));
            }
            start += 2 * len;
        }
        len <<= 1;
    }

    // Multiply by n⁻¹ = 256⁻¹ mod q
    let n_inv = arithmetic::inv_mod(N as u32);
    for j in 0..N {
        f.coeffs[j] = arithmetic::mul(f.coeffs[j], n_inv);
    }

    f
}

//  Vector / Matrix Operations 

/// A vector of polynomials (length determined at runtime from params).
pub type PolyVec = Vec<Poly>;

/// A matrix of polynomials (k rows × l columns).
pub type PolyMat = Vec<Vec<Poly>>;

/// Create a zero vector of `len` polynomials.
pub fn polyvec_zero(len: usize) -> PolyVec {
    (0..len).map(|_| Poly::zero()).collect()
}

/// Coefficient-wise add two polynomial vectors.
pub fn polyvec_add(a: &PolyVec, b: &PolyVec) -> PolyVec {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(ai, bi)| ai.add(bi)).collect()
}

/// Coefficient-wise subtract two polynomial vectors.
pub fn polyvec_sub(a: &PolyVec, b: &PolyVec) -> PolyVec {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(ai, bi)| ai.sub(bi)).collect()
}

/// Forward NTT on each polynomial in a vector.
pub fn polyvec_ntt(v: &PolyVec) -> PolyVec {
    v.iter().map(|p| ntt(p)).collect()
}

/// Inverse NTT on each polynomial in a vector.
pub fn polyvec_inv_ntt(v: &PolyVec) -> PolyVec {
    v.iter().map(|p| inv_ntt(p)).collect()
}

/// Matrix-vector product in NTT domain.
///
/// Computes result[i] = Σ_j mat[i][j] ⊙ vec[j]  for i in 0..k.
///
/// All inputs must be in NTT domain; output is in NTT domain.
/// This is the core operation: t̂ = Â · ŝ₁.
pub fn mat_vec_mul_ntt(mat: &PolyMat, vec: &PolyVec) -> PolyVec {
    let k = mat.len();
    let mut result = polyvec_zero(k);
    for i in 0..k {
        for j in 0..vec.len() {
            let prod = mat[i][j].pointwise_mul(&vec[j]);
            result[i] = result[i].add(&prod);
        }
    }
    result
}

/// Inner product of two polynomial vectors in NTT domain.
///
/// Computes Σ_i a[i] ⊙ b[i].  Used for c · s operations.
pub fn polyvec_pointwise_acc(a: &PolyVec, b: &PolyVec) -> Poly {
    assert_eq!(a.len(), b.len());
    let mut result = Poly::zero();
    for i in 0..a.len() {
        result = result.add(&a[i].pointwise_mul(&b[i]));
    }
    result
}

/// Check that all polynomials in a vector have infinity norm < bound.
pub fn polyvec_check_norm(v: &PolyVec, bound: u32) -> bool {
    v.iter().all(|p| p.check_norm(bound))
}

//  Tests 

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helper: create a polynomial with deterministic coefficients ──

    fn make_poly(seed: u32) -> Poly {
        let mut p = Poly::zero();
        // Simple deterministic fill (not cryptographic, just for testing)
        let mut val = seed;
        for i in 0..N {
            val = val.wrapping_mul(1103515245).wrapping_add(12345);
            p.coeffs[i] = (val >> 16) % Q;
        }
        p
    }

    fn make_small_poly(seed: u32, bound: u32) -> Poly {
        let mut p = Poly::zero();
        let mut val = seed;
        for i in 0..N {
            val = val.wrapping_mul(1103515245).wrapping_add(12345);
            let centered = ((val >> 16) % (2 * bound + 1)) as i32 - bound as i32;
            p.coeffs[i] = arithmetic::from_centered(centered);
        }
        p
    }

    // ── Zeta table ──

    #[test]
    fn test_zetas_first_entry() {
        let zetas = compute_zetas();
        // zetas[0] = ζ^{BitRev8(0)} = ζ^0 = 1
        assert_eq!(zetas[0], 1);
    }

    #[test]
    fn test_zetas_entry_1() {
        let zetas = compute_zetas();
        // zetas[1] = ζ^{BitRev8(1)} = ζ^128
        let expected = arithmetic::pow_mod(ZETA, 128);
        assert_eq!(zetas[1], expected);
    }

    #[test]
    fn test_bit_rev8() {
        assert_eq!(bit_rev8(0), 0);
        assert_eq!(bit_rev8(1), 128);
        assert_eq!(bit_rev8(2), 64);
        assert_eq!(bit_rev8(3), 192);
        assert_eq!(bit_rev8(128), 1);
        assert_eq!(bit_rev8(255), 255);
        // Self-inverse: bit_rev8(bit_rev8(x)) == x
        for x in 0..256u32 {
            assert_eq!(bit_rev8(bit_rev8(x)), x, "bit_rev8 not self-inverse for {}", x);
        }
    }

    #[test]
    fn test_zetas_are_powers_of_zeta() {
        let zetas = compute_zetas();
        // Every entry should be a power of ζ, meaning it satisfies x^512 = 1 mod q
        for k in 0..N {
            let z512 = arithmetic::pow_mod(zetas[k], 512);
            assert_eq!(z512, 1,
                "zetas[{}] = {} is not a 512th root of unity", k, zetas[k]);
        }
    }

    // ── NTT / INTT roundtrip ──

    #[test]
    fn test_ntt_inv_ntt_roundtrip_zero() {
        let p = Poly::zero();
        let p_ntt = ntt(&p);
        let p_back = inv_ntt(&p_ntt);
        for i in 0..N {
            assert_eq!(p_back.coeffs[i], 0, "Zero poly roundtrip failed at {}", i);
        }
    }

    #[test]
    fn test_ntt_inv_ntt_roundtrip_constant() {
        // Constant polynomial: f(X) = 42
        let mut p = Poly::zero();
        p.coeffs[0] = 42;
        let p_ntt = ntt(&p);
        let p_back = inv_ntt(&p_ntt);
        assert_eq!(p_back.coeffs[0], 42, "Constant roundtrip failed at coeff 0");
        for i in 1..N {
            assert_eq!(p_back.coeffs[i], 0, "Constant roundtrip failed at coeff {}", i);
        }
    }

    #[test]
    fn test_ntt_inv_ntt_roundtrip_x() {
        // f(X) = X  (coeffs[1] = 1, rest 0)
        let mut p = Poly::zero();
        p.coeffs[1] = 1;
        let p_ntt = ntt(&p);
        let p_back = inv_ntt(&p_ntt);
        assert_eq!(p_back.coeffs[1], 1);
        for i in 0..N {
            if i != 1 {
                assert_eq!(p_back.coeffs[i], 0, "X roundtrip failed at coeff {}", i);
            }
        }
    }

    #[test]
    fn test_ntt_inv_ntt_roundtrip_random() {
        // Random-looking polynomial
        let p = make_poly(12345);
        let p_ntt = ntt(&p);
        let p_back = inv_ntt(&p_ntt);
        for i in 0..N {
            assert_eq!(p.coeffs[i], p_back.coeffs[i],
                "Random roundtrip failed at index {}", i);
        }
    }

    #[test]
    fn test_ntt_inv_ntt_roundtrip_multiple() {
        // Test with several different polynomials
        for seed in [0u32, 1, 42, 9999, 0xDEAD, 0xBEEF, Q - 1] {
            let p = make_poly(seed);
            let p_back = inv_ntt(&ntt(&p));
            for i in 0..N {
                assert_eq!(p.coeffs[i], p_back.coeffs[i],
                    "Roundtrip failed for seed={} at index {}", seed, i);
            }
        }
    }

    // ── NTT linearity ──

    #[test]
    fn test_ntt_linearity() {
        // NTT(a + b) == NTT(a) + NTT(b)
        let a = make_poly(111);
        let b = make_poly(222);

        let sum_then_ntt = ntt(&a.add(&b));
        let ntt_a = ntt(&a);
        let ntt_b = ntt(&b);
        let ntt_then_sum = ntt_a.add(&ntt_b);

        for i in 0..N {
            assert_eq!(sum_then_ntt.coeffs[i], ntt_then_sum.coeffs[i],
                "NTT linearity failed at index {}", i);
        }
    }

    #[test]
    fn test_ntt_scalar_linearity() {
        // NTT(s·a) == s · NTT(a)
        let a = make_poly(333);
        let s = 42u32;

        let scaled_then_ntt = ntt(&a.scalar_mul(s));
        let ntt_then_scaled = ntt(&a).scalar_mul(s);

        for i in 0..N {
            assert_eq!(scaled_then_ntt.coeffs[i], ntt_then_scaled.coeffs[i],
                "NTT scalar linearity failed at index {}", i);
        }
    }

    // ── NTT multiplication correctness ──

    /// Schoolbook polynomial multiplication in R_q = Z_q[X]/(X^256+1).
    /// This is O(n²) but guaranteed correct — used as the ground truth.
    fn schoolbook_mul(a: &Poly, b: &Poly) -> Poly {
        let mut result = [0u64; 2 * N]; // accumulate in u64 to avoid overflow

        // Standard polynomial multiplication
        for i in 0..N {
            for j in 0..N {
                result[i + j] += a.coeffs[i] as u64 * b.coeffs[j] as u64;
            }
        }

        // Reduce modulo X^256 + 1: coeff[i+256] contributes -coeff[i+256] to coeff[i]
        let mut out = Poly::zero();
        for i in 0..N {
            // result[i] - result[i + N], all mod q
            let pos = result[i] % Q as u64;
            let neg = result[i + N] % Q as u64;
            out.coeffs[i] = arithmetic::sub(pos as u32, neg as u32);
        }
        out
    }

    #[test]
    fn test_ntt_multiplication_simple() {
        // f(X) = 1, g(X) = X  → f·g = X
        let mut f = Poly::zero();
        f.coeffs[0] = 1;
        let mut g = Poly::zero();
        g.coeffs[1] = 1;

        let f_hat = ntt(&f);
        let g_hat = ntt(&g);
        let product_hat = f_hat.pointwise_mul(&g_hat);
        let product = inv_ntt(&product_hat);

        // Expected: X
        assert_eq!(product.coeffs[0], 0);
        assert_eq!(product.coeffs[1], 1);
        for i in 2..N {
            assert_eq!(product.coeffs[i], 0, "1 * X should give X, but coeff[{}]={}", i, product.coeffs[i]);
        }
    }

    #[test]
    fn test_ntt_multiplication_x_times_x() {
        // X · X = X² in R_q
        let mut x = Poly::zero();
        x.coeffs[1] = 1;

        let x_hat = ntt(&x);
        let x2_hat = x_hat.pointwise_mul(&x_hat);
        let x2 = inv_ntt(&x2_hat);

        assert_eq!(x2.coeffs[2], 1, "X² should have coeff[2]=1");
        for i in 0..N {
            if i != 2 {
                assert_eq!(x2.coeffs[i], 0, "X² unexpected coeff at {}", i);
            }
        }
    }

    #[test]
    fn test_ntt_multiplication_x255_times_x() {
        // X^255 · X = X^256 = -1 mod (X^256+1)
        // So result should be: coeff[0] = q-1, rest 0
        let mut x255 = Poly::zero();
        x255.coeffs[255] = 1;
        let mut x1 = Poly::zero();
        x1.coeffs[1] = 1;

        let a_hat = ntt(&x255);
        let b_hat = ntt(&x1);
        let c_hat = a_hat.pointwise_mul(&b_hat);
        let c = inv_ntt(&c_hat);

        assert_eq!(c.coeffs[0], Q - 1, "X^256 should give -1, got {}", c.coeffs[0]);
        for i in 1..N {
            assert_eq!(c.coeffs[i], 0, "X^256 unexpected coeff at {}: {}", i, c.coeffs[i]);
        }
    }

    #[test]
    fn test_ntt_mul_vs_schoolbook() {
        // The gold standard test: NTT multiplication must match schoolbook.
        let a = make_small_poly(42, 100);   // Small coefficients to avoid confusion
        let b = make_small_poly(99, 100);

        // Schoolbook (ground truth)
        let expected = schoolbook_mul(&a, &b);

        // NTT-based
        let a_hat = ntt(&a);
        let b_hat = ntt(&b);
        let c_hat = a_hat.pointwise_mul(&b_hat);
        let got = inv_ntt(&c_hat);

        for i in 0..N {
            assert_eq!(got.coeffs[i], expected.coeffs[i],
                "NTT mul vs schoolbook mismatch at coeff {}: NTT={}, schoolbook={}",
                i, got.coeffs[i], expected.coeffs[i]);
        }
    }

    #[test]
    fn test_ntt_mul_vs_schoolbook_large() {
        // Larger coefficients
        let a = make_poly(0xCAFE);
        let b = make_poly(0xBABE);

        let expected = schoolbook_mul(&a, &b);

        let a_hat = ntt(&a);
        let b_hat = ntt(&b);
        let c_hat = a_hat.pointwise_mul(&b_hat);
        let got = inv_ntt(&c_hat);

        for i in 0..N {
            assert_eq!(got.coeffs[i], expected.coeffs[i],
                "NTT mul vs schoolbook (large) mismatch at coeff {}", i);
        }
    }

    #[test]
    fn test_ntt_mul_vs_schoolbook_multiple() {
        // Test several pairs
        let seeds = [(1u32, 2u32), (42, 99), (1000, 2000), (Q - 1, 1)];
        for (sa, sb) in seeds {
            let a = make_poly(sa);
            let b = make_poly(sb);
            let expected = schoolbook_mul(&a, &b);
            let got = inv_ntt(&ntt(&a).pointwise_mul(&ntt(&b)));
            for i in 0..N {
                assert_eq!(got.coeffs[i], expected.coeffs[i],
                    "NTT vs schoolbook failed for seeds ({}, {}) at coeff {}", sa, sb, i);
            }
        }
    }

    // ── Poly basic operations ──

    #[test]
    fn test_poly_add_sub_inverse() {
        let a = make_poly(100);
        let zero = a.sub(&a);
        for i in 0..N {
            assert_eq!(zero.coeffs[i], 0, "a - a should be 0 at coeff {}", i);
        }
    }

    #[test]
    fn test_poly_add_commutative() {
        let a = make_poly(200);
        let b = make_poly(300);
        let ab = a.add(&b);
        let ba = b.add(&a);
        for i in 0..N {
            assert_eq!(ab.coeffs[i], ba.coeffs[i]);
        }
    }

    #[test]
    fn test_check_norm() {
        let mut p = Poly::zero();
        p.coeffs[0] = 5;
        p.coeffs[1] = Q - 5; // centered = -5

        assert!(p.check_norm(6), "norm should be < 6");
        assert!(!p.check_norm(5), "norm should not be < 5 (it equals 5)");
        assert!(!p.check_norm(4), "norm should not be < 4");
    }

    // ── Vector operations ──

    #[test]
    fn test_polyvec_add_sub() {
        let a = vec![make_poly(10), make_poly(20)];
        let b = vec![make_poly(30), make_poly(40)];

        let sum = polyvec_add(&a, &b);
        let diff = polyvec_sub(&sum, &b);

        for k in 0..2 {
            for i in 0..N {
                assert_eq!(diff[k].coeffs[i], a[k].coeffs[i],
                    "polyvec add-sub inverse failed");
            }
        }
    }

    #[test]
    fn test_polyvec_ntt_inv_ntt_roundtrip() {
        let v = vec![make_poly(50), make_poly(60), make_poly(70)];
        let v_ntt = polyvec_ntt(&v);
        let v_back = polyvec_inv_ntt(&v_ntt);

        for k in 0..3 {
            for i in 0..N {
                assert_eq!(v[k].coeffs[i], v_back[k].coeffs[i],
                    "polyvec NTT roundtrip failed at vec[{}][{}]", k, i);
            }
        }
    }

    // ── Matrix-vector multiply ──

    #[test]
    fn test_mat_vec_mul_identity() {
        // "Identity-like" test: 1×1 matrix with identity poly, vector of length 1.
        // A = [[1]], v = [p] → result = [1·p] = [p]
        let p = make_poly(500);
        let mut one = Poly::zero();
        one.coeffs[0] = 1;

        let mat = vec![vec![ntt(&one)]];
        let vec = vec![ntt(&p)];
        let result_ntt = mat_vec_mul_ntt(&mat, &vec);
        let result = inv_ntt(&result_ntt[0]);

        for i in 0..N {
            assert_eq!(result.coeffs[i], p.coeffs[i],
                "Identity mat-vec mul failed at coeff {}", i);
        }
    }

    #[test]
    fn test_polyvec_check_norm() {
        let mut p1 = Poly::zero();
        p1.coeffs[0] = 3;
        let mut p2 = Poly::zero();
        p2.coeffs[0] = Q - 3; // -3

        let v = vec![p1, p2];
        assert!(polyvec_check_norm(&v, 4), "Norm should be < 4");
        assert!(!polyvec_check_norm(&v, 3), "Norm should not be < 3");
    }
}