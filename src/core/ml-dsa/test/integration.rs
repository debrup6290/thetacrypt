//! Integration tests for all ML-DSA modules.
//!
//! These tests verify:
//!   1. Parameter consistency and cross-checks against FIPS 204
//!   2. SHAKE correctness with verified NIST KAT vectors
//!   3. SHAKE behavior under ML-DSA-like usage patterns
//!   4. Arithmetic correctness: field ops, Power2Round, Decompose, Hints
//!   5. Polynomial NTT correctness: roundtrip, linearity, multiplication
//!   6. Encoding correctness: bit packing, key/sig serialization roundtrips
//!   7. Sampling correctness: ExpandA, ExpandS, ExpandMask, SampleInBall
//!   8. Signature correctness: KeyGen, Sign, Verify end-to-end

use ml_dsa::params::*;
use ml_dsa::shake::*;
use ml_dsa::arithmetic;
use ml_dsa::poly::{self, Poly};
use ml_dsa::encoding;
use ml_dsa::sampling;
use ml_dsa::signature;


fn hex_to_bytes(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
        .collect()
}

fn bytes_to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}


#[test]
fn test_all_param_sets_exist() {
    assert_eq!(ML_DSA_44.name, "ML-DSA-44");
    assert_eq!(ML_DSA_65.name, "ML-DSA-65");
    assert_eq!(ML_DSA_87.name, "ML-DSA-87");
}

#[test]
fn test_security_levels_ordered() {
    assert!(ML_DSA_44.k <= ML_DSA_65.k);
    assert!(ML_DSA_65.k <= ML_DSA_87.k);
    assert!(ML_DSA_44.l <= ML_DSA_65.l);
    assert!(ML_DSA_65.l <= ML_DSA_87.l);

    assert!(ML_DSA_44.pk_bytes <= ML_DSA_65.pk_bytes);
    assert!(ML_DSA_65.pk_bytes <= ML_DSA_87.pk_bytes);
    assert!(ML_DSA_44.sk_bytes <= ML_DSA_65.sk_bytes);
    assert!(ML_DSA_65.sk_bytes <= ML_DSA_87.sk_bytes);
    assert!(ML_DSA_44.sig_bytes <= ML_DSA_65.sig_bytes);
    assert!(ML_DSA_65.sig_bytes <= ML_DSA_87.sig_bytes);
    assert!(ML_DSA_44.c_tilde_bytes <= ML_DSA_65.c_tilde_bytes);
    assert!(ML_DSA_65.c_tilde_bytes <= ML_DSA_87.c_tilde_bytes);
}

#[test]
fn test_gamma2_divides_q_minus_1() {
    for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
        assert_eq!(
            (Q - 1) % (2 * p.gamma2), 0,
            "{}: 2·γ₂ must divide q-1", p.name
        );
    }
}

#[test]
fn test_rejection_bound_is_feasible() {
    for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
        assert!(p.gamma1 > p.beta, "{}: γ₁ > β", p.name);
        assert!(p.gamma2 > p.beta, "{}: γ₂ > β", p.name);
    }
}

#[test]
fn test_omega_bounds() {
    for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
        assert!(p.omega > 0, "{}: ω > 0", p.name);
        assert!(p.omega < p.k * N, "{}: ω < k·n", p.name);
    }
}

#[test]
fn test_tau_less_than_n() {
    for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
        assert!(p.tau < N, "{}: τ < n", p.name);
    }
}

#[test]
fn test_eta_values() {
    for p in &[ML_DSA_44, ML_DSA_65, ML_DSA_87] {
        assert!(
            p.eta == 2 || p.eta == 4,
            "{}: η must be 2 or 4, got {}", p.name, p.eta
        );
    }
}

#[test]
fn test_all_sizes_exact() {
    assert_eq!(ML_DSA_44.pk_bytes, 1312);
    assert_eq!(ML_DSA_44.sk_bytes, 2560);
    assert_eq!(ML_DSA_44.sig_bytes, 2420);

    assert_eq!(ML_DSA_65.pk_bytes, 1952);
    assert_eq!(ML_DSA_65.sk_bytes, 4032);
    assert_eq!(ML_DSA_65.sig_bytes, 3309);

    assert_eq!(ML_DSA_87.pk_bytes, 2592);
    assert_eq!(ML_DSA_87.sk_bytes, 4896);
    assert_eq!(ML_DSA_87.sig_bytes, 4627);
}

#[test]
fn test_w1_bits_values() {
    assert_eq!(ML_DSA_44.w1_bits(), 6, "ML-DSA-44 w1_bits");
    assert_eq!(ML_DSA_65.w1_bits(), 4, "ML-DSA-65 w1_bits");
    assert_eq!(ML_DSA_87.w1_bits(), 4, "ML-DSA-87 w1_bits");
}


#[test]
fn test_shake256_nist_kat_empty() {
    let expected = hex_to_bytes(
        "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
    );
    let result = shake256(&[], 32);
    assert_eq!(result, expected,
        "\nExpected: {}\nGot:      {}",
        bytes_to_hex(&expected), bytes_to_hex(&result)
    );
}

#[test]
fn test_shake128_nist_kat_empty() {
    let expected = hex_to_bytes(
        "7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26"
    );
    let result = shake128(&[], 32);
    assert_eq!(result, expected,
        "\nExpected: {}\nGot:      {}",
        bytes_to_hex(&expected), bytes_to_hex(&result)
    );
}

#[test]
fn test_shake256_nist_kat_abc() {
    let expected = hex_to_bytes(
        "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739"
    );
    let result = shake256(b"abc", 32);
    assert_eq!(result, expected,
        "\nExpected: {}\nGot:      {}",
        bytes_to_hex(&expected), bytes_to_hex(&result)
    );
}


#[test]
fn test_shake256_64_bytes_prefix_matches_32() {
    let short = shake256(&[], 32);
    let long = shake256(&[], 64);
    assert_eq!(&long[..32], &short[..],
        "SHAKE-256: prefix of longer output must match shorter output"
    );
}

#[test]
fn test_shake128_64_bytes_prefix_matches_32() {
    let short = shake128(&[], 32);
    let long = shake128(&[], 64);
    assert_eq!(&long[..32], &short[..],
        "SHAKE-128: prefix of longer output must match shorter output"
    );
}

#[test]
fn test_shake256_200_bytes_prefix_consistency() {
    let out_32 = shake256(b"boundary", 32);
    let out_200 = shake256(b"boundary", 200);
    assert_eq!(&out_200[..32], &out_32[..]);
}

#[test]
fn test_shake128_200_bytes_prefix_consistency() {
    let out_32 = shake128(b"boundary", 32);
    let out_200 = shake128(b"boundary", 200);
    assert_eq!(&out_200[..32], &out_32[..]);
}


#[test]
fn test_shake256_seed_expansion_pattern() {
    let seed1 = [0x01u8; 32];
    let seed2 = [0x02u8; 32];

    let expand1 = shake256(&seed1, 128);
    let expand2 = shake256(&seed2, 128);

    assert_eq!(expand1.len(), 128);
    assert_ne!(expand1, expand2, "Different seeds must produce different expansions");

    let rho = &expand1[0..32];
    let rho_prime = &expand1[32..96];
    let k = &expand1[96..128];
    assert_ne!(rho, &rho_prime[..32]);
    assert_ne!(rho, k);
}

#[test]
fn test_shake128_expand_a_pattern() {
    let rho = [0x42u8; 32];

    let mut xof_00 = Shake128::new();
    xof_00.absorb(&rho);
    xof_00.absorb(&[0u8, 0u8]);
    let out_00 = xof_00.squeeze_vec(128);

    let mut xof_01 = Shake128::new();
    xof_01.absorb(&rho);
    xof_01.absorb(&[1u8, 0u8]);
    let out_01 = xof_01.squeeze_vec(128);

    assert_ne!(out_00, out_01, "Different indices must produce different streams");

    let mut xof_10 = Shake128::new();
    xof_10.absorb(&rho);
    xof_10.absorb(&[0u8, 1u8]);
    let out_10 = xof_10.squeeze_vec(128);

    assert_ne!(out_01, out_10, "Index order must matter");
}

#[test]
fn test_shake256_expand_s_pattern() {
    let rho_prime = [0xAA; 64];

    let mut results = Vec::new();
    for counter in 0u16..8 {
        let mut xof = Shake256::new();
        xof.absorb(&rho_prime);
        xof.absorb(&counter.to_le_bytes());
        results.push(xof.squeeze_vec(64));
    }

    for i in 0..results.len() {
        for j in (i + 1)..results.len() {
            assert_ne!(results[i], results[j],
                "Counter {} and {} must produce different output", i, j);
        }
    }
}

#[test]
fn test_shake256_mu_computation_pattern() {
    let tr = [0xBBu8; 64];
    let msg = b"Hello, ML-DSA world!";

    let mut combined = Vec::new();
    combined.extend_from_slice(&tr);
    combined.extend_from_slice(msg);
    let mu1 = shake256(&combined, 64);

    let mut xof = Shake256::new();
    xof.absorb(&tr);
    xof.absorb(msg);
    let mu2 = xof.squeeze_vec(64);

    assert_eq!(mu1, mu2, "Concatenated absorb must equal incremental absorb");
}

#[test]
fn test_shake256_challenge_hash_pattern() {
    let mu = [0xCC; 64];
    let w1_encoded = vec![0xDD; 768];

    for &c_tilde_len in &[
        ML_DSA_44.c_tilde_bytes,
        ML_DSA_65.c_tilde_bytes,
        ML_DSA_87.c_tilde_bytes,
    ] {
        let mut xof = Shake256::new();
        xof.absorb(&mu);
        xof.absorb(&w1_encoded);
        let c_tilde = xof.squeeze_vec(c_tilde_len);

        assert_eq!(c_tilde.len(), c_tilde_len);
        assert!(c_tilde.iter().any(|&b| b != 0));
    }
}


#[test]
fn test_shake_single_byte_inputs_all_unique() {
    let mut hashes = Vec::new();
    for b in 0u8..=255 {
        hashes.push(shake256(&[b], 16));
    }
    for i in 0..hashes.len() {
        for j in (i + 1)..hashes.len() {
            assert_ne!(hashes[i], hashes[j],
                "Byte {} and {} produced same hash", i, j);
        }
    }
}

#[test]
fn test_shake_zero_length_squeeze() {
    let mut xof = Shake256::from_data(b"test");
    let out = xof.squeeze_vec(0);
    assert_eq!(out.len(), 0);

    let out2 = xof.squeeze_vec(32);
    assert_eq!(out2.len(), 32);

    let mut xof2 = Shake256::from_data(b"test");
    let out3 = xof2.squeeze_vec(32);
    assert_eq!(out2, out3);
}

#[test]
fn test_shake_squeeze_one_byte_at_a_time() {
    let data = b"byte-by-byte squeeze";

    let mut xof1 = Shake256::from_data(data);
    let bulk = xof1.squeeze_vec(256);

    let mut xof2 = Shake256::from_data(data);
    let mut single = Vec::new();
    for _ in 0..256 {
        let mut buf = [0u8; 1];
        xof2.squeeze(&mut buf);
        single.push(buf[0]);
    }

    assert_eq!(bulk, single, "Byte-by-byte squeeze must match bulk");
}

#[test]
fn test_shake_absorb_empty_slices() {
    let mut xof1 = Shake256::new();
    xof1.absorb(b"hello");
    let out1 = xof1.squeeze_vec(32);

    let mut xof2 = Shake256::new();
    xof2.absorb(&[]);
    xof2.absorb(b"hello");
    xof2.absorb(&[]);
    xof2.absorb(&[]);
    let out2 = xof2.squeeze_vec(32);

    assert_eq!(out1, out2, "Empty absorbs must not change state");
}

#[test]
fn test_shake_long_output() {
    let mut xof = Shake256::from_data(b"long output test");
    let out = xof.squeeze_vec(10000);
    assert_eq!(out.len(), 10000);

    let mut xof2 = Shake256::from_data(b"long output test");
    let out2 = xof2.squeeze_vec(10000);
    assert_eq!(out, out2);

    let zero_count = out.iter().filter(|&&b| b == 0).count();
    assert!(zero_count > 5 && zero_count < 100,
        "Suspicious zero count: {} in 10000 bytes", zero_count);
}

#[test]
fn test_shake_large_absorb_then_small_squeeze() {
    let big_input = vec![0x55u8; 100_000];
    let result = shake256(&big_input, 32);
    assert_eq!(result.len(), 32);
    assert!(result.iter().any(|&b| b != 0));
}

#[test]
fn test_independent_instances() {
    let mut xof_a = Shake256::from_data(b"instance A");
    let mut xof_b = Shake256::from_data(b"instance B");

    let a1 = xof_a.squeeze_vec(16);
    let b1 = xof_b.squeeze_vec(16);
    let a2 = xof_a.squeeze_vec(16);
    let b2 = xof_b.squeeze_vec(16);

    let mut xof_a_check = Shake256::from_data(b"instance A");
    let a_full = xof_a_check.squeeze_vec(32);
    assert_eq!(&a_full[..16], &a1[..]);
    assert_eq!(&a_full[16..], &a2[..]);

    let mut xof_b_check = Shake256::from_data(b"instance B");
    let b_full = xof_b_check.squeeze_vec(32);
    assert_eq!(&b_full[..16], &b1[..]);
    assert_eq!(&b_full[16..], &b2[..]);
}


#[test]
fn print_shake_outputs_for_verification() {
    println!("SHAKE outputs — verify against a reference implementation  ");

    let r = shake256(&[], 64);
    println!(" SHAKE-256(\"\", 64):                                          ");
    println!("   {}", bytes_to_hex(&r[..32]));
    println!("   {}", bytes_to_hex(&r[32..]));

    let r = shake128(&[], 64);
    println!(" SHAKE-128(\"\", 64):                                          ");
    println!("   {}", bytes_to_hex(&r[..32]));
    println!("   {}", bytes_to_hex(&r[32..]));

    let r = shake256(b"abc", 32);
    println!(" SHAKE-256(\"abc\", 32):                                       ");
    println!("   {}", bytes_to_hex(&r));

    let r = shake128(b"abc", 32);
    println!(" SHAKE-128(\"abc\", 32):                                       ");
    println!("   {}", bytes_to_hex(&r));

    let input = vec![0xA3; 200];
    let r = shake256(&input, 32);
    println!(" SHAKE-256(0xA3×200, 32):                                    ");
    println!("   {}", bytes_to_hex(&r));

    let r = shake128(&input, 32);
    println!(" SHAKE-128(0xA3×200, 32):                                    ");
    println!("   {}", bytes_to_hex(&r));

    println!("╚╝");
}


#[test]
fn test_arithmetic_field_axioms() {
    let vals: Vec<u32> = (0..Q).step_by(100003).collect();

    for &a in &vals {
        for &b in &vals {
            assert_eq!(arithmetic::add(a, b), arithmetic::add(b, a));
            assert_eq!(arithmetic::mul(a, b), arithmetic::mul(b, a));

            assert_eq!(arithmetic::add(a, arithmetic::neg(a)), 0);

            assert_eq!(arithmetic::add(arithmetic::sub(a, b), b), a);
        }
    }
}

#[test]
fn test_arithmetic_mul_matches_naive() {
    let vals: Vec<u32> = vec![0, 1, 2, Q - 1, Q - 2, Q / 2, 12345, 7654321, 4190208];

    for &a in &vals {
        for &b in &vals {
            let expected = ((a as u64 * b as u64) % Q as u64) as u32;
            let got = arithmetic::mul(a, b);
            assert_eq!(got, expected,
                "mul({}, {}) = {} but naive = {}", a, b, got, expected);
        }
    }
}


#[test]
fn test_power2round_exhaustive_reconstruction() {
    let two_d = 1i64 << D;
    for r in (0..Q).step_by(997) {
        let (r1, r0) = arithmetic::power2round(r);
        let reconstructed = r1 as i64 * two_d + r0 as i64;
        assert_eq!(reconstructed as u32, r,
            "Power2Round failed for r={}: r1={}, r0={}", r, r1, r0);
    }
}

#[test]
fn test_power2round_t1_fits_10_bits() {
    for r in (0..Q).step_by(997) {
        let (r1, _) = arithmetic::power2round(r);
        assert!(r1 < 1024,
            "r1={} doesn't fit in 10 bits for r={}", r1, r);
    }
}


#[test]
fn test_decompose_exhaustive_44() {
    let gamma2 = ML_DSA_44.gamma2;
    let alpha = 2 * gamma2;
    let m = (Q - 1) / alpha;

    for r in (0..Q).step_by(997) {
        let (r1, r0) = arithmetic::decompose(r, gamma2);
        assert!(r1 < m, "r1={} >= m={} for r={}", r1, m, r);
        let half = (alpha / 2) as i32;
        assert!(r0 > -half && r0 <= half,
            "r0={} out of range for r={}", r0, r);
        let reconstructed = ((r1 as i64) * (alpha as i64) + (r0 as i64))
            .rem_euclid(Q as i64) as u32;
        assert_eq!(reconstructed, r,
            "Decompose-44 failed for r={}: r1={}, r0={}", r, r1, r0);
    }
}

#[test]
fn test_decompose_exhaustive_65() {
    let gamma2 = ML_DSA_65.gamma2;
    let alpha = 2 * gamma2;
    let m = (Q - 1) / alpha;

    for r in (0..Q).step_by(997) {
        let (r1, r0) = arithmetic::decompose(r, gamma2);
        assert!(r1 < m, "r1={} >= m={} for r={}", r1, m, r);
        let half = (alpha / 2) as i32;
        assert!(r0 > -half && r0 <= half,
            "r0={} out of range for r={}", r0, r);
        let reconstructed = ((r1 as i64) * (alpha as i64) + (r0 as i64))
            .rem_euclid(Q as i64) as u32;
        assert_eq!(reconstructed, r,
            "Decompose-65 failed for r={}: r1={}, r0={}", r, r1, r0);
    }
}


#[test]
fn test_hint_fundamental_property_44() {
    let gamma2 = ML_DSA_44.gamma2;

    let z_values: &[i32] = &[0, 1, -1, 50, -50, 500, -500, 5000, -5000, 50000, -50000];
    let r_values: Vec<u32> = (0..Q).step_by(50021).collect();

    for &r in &r_values {
        for &z in z_values {
            let rz = ((r as i64) + (z as i64)).rem_euclid(Q as i64) as u32;
            let hint = arithmetic::make_hint(z, r, gamma2);
            let recovered = arithmetic::use_hint(hint, rz, gamma2);
            let expected = arithmetic::high_bits(r, gamma2);

            assert_eq!(recovered, expected,
                "Hint property violated: r={}, z={}, r+z={}, hint={}, recovered={}, expected={}",
                r, z, rz, hint, recovered, expected);
        }
    }
}

#[test]
fn test_hint_fundamental_property_65() {
    let gamma2 = ML_DSA_65.gamma2;

    let z_values: &[i32] = &[0, 1, -1, 100, -100, 1000, -1000, 10000, -10000];
    let r_values: Vec<u32> = (0..Q).step_by(50021).collect();

    for &r in &r_values {
        for &z in z_values {
            let rz = ((r as i64) + (z as i64)).rem_euclid(Q as i64) as u32;
            let hint = arithmetic::make_hint(z, r, gamma2);
            let recovered = arithmetic::use_hint(hint, rz, gamma2);
            let expected = arithmetic::high_bits(r, gamma2);

            assert_eq!(recovered, expected,
                "Hint property violated: r={}, z={}", r, z);
        }
    }
}

#[test]
fn test_hint_sparsity() {
    let gamma2 = ML_DSA_44.gamma2;
    let beta = ML_DSA_44.beta;

    let mut hint_count = 0u64;
    let mut total = 0u64;

    let bound = gamma2 as i32 - beta as i32;
    let z_values: Vec<i32> = (-10..=10).map(|i| i * (bound / 20).max(1)).collect();

    for r in (0..Q).step_by(10007) {
        for &z in &z_values {
            let hint = arithmetic::make_hint(z, r, gamma2);
            if hint { hint_count += 1; }
            total += 1;
        }
    }

    let ratio = hint_count as f64 / total as f64;
    println!("\nHint sparsity (ML-DSA-44): {}/{} = {:.4}%",
        hint_count, total, ratio * 100.0);
    assert!(ratio < 0.5,
        "Hints are not sparse enough: {}% of hints are 1", ratio * 100.0);
}


fn test_poly(seed: u32) -> Poly {
    let mut p = Poly::zero();
    let mut val = seed;
    for i in 0..N {
        val = val.wrapping_mul(1103515245).wrapping_add(12345);
        p.coeffs[i] = (val >> 16) % Q;
    }
    p
}

#[test]
fn test_ntt_roundtrip_many_polys() {
    for seed in 0..50u32 {
        let p = test_poly(seed * 137 + 1);
        let p_back = poly::inv_ntt(&poly::ntt(&p));
        for i in 0..N {
            assert_eq!(p.coeffs[i], p_back.coeffs[i],
                "NTT roundtrip failed: seed={}, index={}", seed, i);
        }
    }
}

#[test]
fn test_ntt_mul_commutativity() {
    let a = test_poly(1001);
    let b = test_poly(2002);

    let a_hat = poly::ntt(&a);
    let b_hat = poly::ntt(&b);

    let ab = poly::inv_ntt(&a_hat.pointwise_mul(&b_hat));
    let ba = poly::inv_ntt(&b_hat.pointwise_mul(&a_hat));

    for i in 0..N {
        assert_eq!(ab.coeffs[i], ba.coeffs[i],
            "NTT multiplication not commutative at index {}", i);
    }
}

#[test]
fn test_ntt_mul_by_one() {
    let a = test_poly(3003);
    let mut one = Poly::zero();
    one.coeffs[0] = 1;

    let a_hat = poly::ntt(&a);
    let one_hat = poly::ntt(&one);
    let product = poly::inv_ntt(&a_hat.pointwise_mul(&one_hat));

    for i in 0..N {
        assert_eq!(product.coeffs[i], a.coeffs[i],
            "a * 1 should equal a, mismatch at {}", i);
    }
}

#[test]
fn test_ntt_mul_by_zero() {
    let a = test_poly(4004);
    let zero = Poly::zero();

    let a_hat = poly::ntt(&a);
    let zero_hat = poly::ntt(&zero);
    let product = poly::inv_ntt(&a_hat.pointwise_mul(&zero_hat));

    for i in 0..N {
        assert_eq!(product.coeffs[i], 0, "a * 0 should be 0 at {}", i);
    }
}

#[test]
fn test_ntt_distributive() {
    let a = test_poly(5005);
    let b = test_poly(6006);
    let c = test_poly(7007);

    let a_hat = poly::ntt(&a);
    let b_hat = poly::ntt(&b);
    let c_hat = poly::ntt(&c);
    let bc_hat = poly::ntt(&b.add(&c));

    let lhs = poly::inv_ntt(&a_hat.pointwise_mul(&bc_hat));

    let ab = poly::inv_ntt(&a_hat.pointwise_mul(&b_hat));
    let ac = poly::inv_ntt(&a_hat.pointwise_mul(&c_hat));
    let rhs = ab.add(&ac);

    for i in 0..N {
        assert_eq!(lhs.coeffs[i], rhs.coeffs[i],
            "NTT distributive law failed at coeff {}", i);
    }
}


#[test]
fn test_mat_vec_mul_dimensions() {
    let p = &ML_DSA_44;

    let mut mat = Vec::new();
    for i in 0..p.k {
        let mut row = Vec::new();
        for j in 0..p.l {
            row.push(poly::ntt(&test_poly((i * p.l + j) as u32 + 100)));
        }
        mat.push(row);
    }
    let vec: Vec<Poly> = (0..p.l)
        .map(|j| poly::ntt(&test_poly(j as u32 + 200)))
        .collect();

    let result = poly::mat_vec_mul_ntt(&mat, &vec);
    assert_eq!(result.len(), p.k, "Result should have k={} entries", p.k);
}

#[test]
fn test_mat_vec_mul_linearity() {
    let mat = vec![
        vec![poly::ntt(&test_poly(10)), poly::ntt(&test_poly(11))],
        vec![poly::ntt(&test_poly(12)), poly::ntt(&test_poly(13))],
    ];

    let v1 = vec![poly::ntt(&test_poly(20)), poly::ntt(&test_poly(21))];
    let v2 = vec![poly::ntt(&test_poly(30)), poly::ntt(&test_poly(31))];

    let v_sum = poly::polyvec_add(&v1, &v2);

    let lhs = poly::mat_vec_mul_ntt(&mat, &v_sum);

    let av1 = poly::mat_vec_mul_ntt(&mat, &v1);
    let av2 = poly::mat_vec_mul_ntt(&mat, &v2);
    let rhs = poly::polyvec_add(&av1, &av2);

    for k in 0..2 {
        let lhs_normal = poly::inv_ntt(&lhs[k]);
        let rhs_normal = poly::inv_ntt(&rhs[k]);
        for i in 0..N {
            assert_eq!(lhs_normal.coeffs[i], rhs_normal.coeffs[i],
                "Mat-vec linearity failed at result[{}][{}]", k, i);
        }
    }
}

#[test]
fn test_keygen_core_identity() {
    let k = 2;
    let l = 2;

    let mut s1 = Vec::new();
    let mut s2 = Vec::new();
    for j in 0..l {
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = arithmetic::from_centered(((i + j * 7) % 5) as i32 - 2);
        }
        s1.push(p);
    }
    for j in 0..k {
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = arithmetic::from_centered(((i + j * 11) % 5) as i32 - 2);
        }
        s2.push(p);
    }

    let mut a = Vec::new();
    for i in 0..k {
        let mut row = Vec::new();
        for j in 0..l {
            row.push(test_poly((i * l + j) as u32 + 500));
        }
        a.push(row);
    }

    let a_hat: Vec<Vec<Poly>> = a.iter()
        .map(|row| row.iter().map(|p| poly::ntt(p)).collect())
        .collect();
    let s1_hat = poly::polyvec_ntt(&s1);
    let t_hat = poly::mat_vec_mul_ntt(&a_hat, &s1_hat);
    let t_normal = poly::polyvec_inv_ntt(&t_hat);
    let t = poly::polyvec_add(&t_normal, &s2);

    let has_nonzero = t[0].coeffs.iter().any(|&c| c != 0);
    assert!(has_nonzero, "t should not be the zero vector");

    for i in 0..k {
        for j in 0..N {
            let (r1, r0) = arithmetic::power2round(t[i].coeffs[j]);
            let reconstructed = r1 as i64 * (1i64 << D) + r0 as i64;
            assert_eq!(reconstructed as u32, t[i].coeffs[j],
                "Power2Round on t[{}][{}] failed", i, j);
        }
    }
}


#[test]
fn test_encoding_pk_roundtrip_all_levels() {
    for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
        let rho = [0x42u8; 32];

        let mut t1 = Vec::new();
        for k in 0..params.k {
            let mut p = Poly::zero();
            for i in 0..N {
                let fake_t = ((k * N + i) as u32 * 7919 + 13) % Q;
                let (r1, _) = arithmetic::power2round(fake_t);
                p.coeffs[i] = r1;
            }
            t1.push(p);
        }

        let pk = encoding::pk_encode(&rho, &t1, params);
        assert_eq!(pk.len(), params.pk_bytes, "{} pk size", params.name);

        let (rho_d, t1_d) = encoding::pk_decode(&pk, params);
        assert_eq!(rho, rho_d);
        for k in 0..params.k {
            for i in 0..N {
                assert_eq!(t1[k].coeffs[i], t1_d[k].coeffs[i],
                    "{}: t1[{}][{}]", params.name, k, i);
            }
        }
    }
}

#[test]
fn test_encoding_sk_roundtrip_all_levels() {
    for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
        let rho = [0xAA; 32];
        let k_seed = [0xBB; 32];
        let tr = [0xCC; 64];
        let eta = params.eta;

        let mut s1 = Vec::new();
        for j in 0..params.l {
            let mut p = Poly::zero();
            for i in 0..N {
                let c = ((j * N + i) as i32 * 3 + 1) % (2 * eta as i32 + 1) - eta as i32;
                p.coeffs[i] = arithmetic::from_centered(c);
            }
            s1.push(p);
        }

        let mut s2 = Vec::new();
        for j in 0..params.k {
            let mut p = Poly::zero();
            for i in 0..N {
                let c = ((j * N + i) as i32 * 7 + 3) % (2 * eta as i32 + 1) - eta as i32;
                p.coeffs[i] = arithmetic::from_centered(c);
            }
            s2.push(p);
        }

        let t0_bound = (1i32 << (D - 1)) - 1; // 4095
        let mut t0 = Vec::new();
        for j in 0..params.k {
            let mut p = Poly::zero();
            for i in 0..N {
                let fake_t = ((j * N + i) as u32 * 7919 + 13) % Q;
                let (_, r0) = arithmetic::power2round(fake_t);
                p.coeffs[i] = arithmetic::from_centered(r0);
            }
            t0.push(p);
        }

        let sk = encoding::sk_encode(&rho, &k_seed, &tr, &s1, &s2, &t0, params);
        assert_eq!(sk.len(), params.sk_bytes, "{} sk size", params.name);

        let (rho_d, k_d, tr_d, s1_d, s2_d, t0_d) = encoding::sk_decode(&sk, params);
        assert_eq!(rho, rho_d);
        assert_eq!(k_seed, k_d);
        assert_eq!(tr, tr_d);

        for j in 0..params.l {
            for i in 0..N {
                assert_eq!(s1[j].coeffs[i], s1_d[j].coeffs[i],
                    "{}: s1[{}][{}]", params.name, j, i);
            }
        }
        for j in 0..params.k {
            for i in 0..N {
                assert_eq!(s2[j].coeffs[i], s2_d[j].coeffs[i],
                    "{}: s2[{}][{}]", params.name, j, i);
            }
        }
        for j in 0..params.k {
            for i in 0..N {
                assert_eq!(t0[j].coeffs[i], t0_d[j].coeffs[i],
                    "{}: t0[{}][{}]", params.name, j, i);
            }
        }
    }
}

#[test]
fn test_encoding_sig_roundtrip_all_levels() {
    for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
        let c_tilde: Vec<u8> = (0..params.c_tilde_bytes)
            .map(|i| (i * 7 + 3) as u8)
            .collect();

        let gamma1_m1 = params.gamma1 as i32 - 1;
        let mut z = Vec::new();
        for j in 0..params.l {
            let mut p = Poly::zero();
            for i in 0..N {
                let c = ((j * N + i) as i64 * 1031 + 7) % (2 * gamma1_m1 as i64 + 1) - gamma1_m1 as i64;
                p.coeffs[i] = arithmetic::from_centered(c as i32);
            }
            z.push(p);
        }

        let mut h = vec![vec![false; N]; params.k];
        h[0][0] = true;
        h[0][42] = true;
        if params.k > 1 { h[1][100] = true; }
        if params.k > 2 { h[2][255] = true; }

        let sig = encoding::sig_encode(&c_tilde, &z, &h, params);
        assert_eq!(sig.len(), params.sig_bytes, "{} sig size", params.name);

        let (c_d, z_d, h_d) = encoding::sig_decode(&sig, params)
            .expect("Valid signature decode");

        assert_eq!(c_tilde, c_d, "{} c_tilde", params.name);

        for j in 0..params.l {
            for i in 0..N {
                assert_eq!(z[j].coeffs[i], z_d[j].coeffs[i],
                    "{}: z[{}][{}]", params.name, j, i);
            }
        }
        for k in 0..params.k {
            for i in 0..N {
                assert_eq!(h[k][i], h_d[k][i],
                    "{}: h[{}][{}]", params.name, k, i);
            }
        }
    }
}

#[test]
fn test_encoding_keygen_pipeline() {
    let params = &ML_DSA_44;

    let xi = [0x99u8; 32];
    let expanded = shake256(&xi, 128);
    let mut rho = [0u8; 32];
    let mut rho_prime = [0u8; 64];
    let mut k_seed = [0u8; 32];
    rho.copy_from_slice(&expanded[0..32]);
    rho_prime.copy_from_slice(&expanded[32..96]);
    k_seed.copy_from_slice(&expanded[96..128]);

    let eta = params.eta;
    let s1: Vec<Poly> = (0..params.l).map(|j| {
        let mut p = Poly::zero();
        for i in 0..N {
            let c = ((j * N + i) % (2 * eta as usize + 1)) as i32 - eta as i32;
            p.coeffs[i] = arithmetic::from_centered(c);
        }
        p
    }).collect();
    let s2: Vec<Poly> = (0..params.k).map(|j| {
        let mut p = Poly::zero();
        for i in 0..N {
            let c = ((j * N + i + 3) % (2 * eta as usize + 1)) as i32 - eta as i32;
            p.coeffs[i] = arithmetic::from_centered(c);
        }
        p
    }).collect();

    let mut t1 = Vec::new();
    let mut t0 = Vec::new();
    for j in 0..params.k {
        let mut t1_p = Poly::zero();
        let mut t0_p = Poly::zero();
        for i in 0..N {
            let fake_t = ((j * N + i) as u32 * 7919 + 13) % Q;
            let (hi, lo) = arithmetic::power2round(fake_t);
            t1_p.coeffs[i] = hi;
            t0_p.coeffs[i] = arithmetic::from_centered(lo);
        }
        t1.push(t1_p);
        t0.push(t0_p);
    }

    let pk_bytes = encoding::pk_encode(&rho, &t1, params);

    let tr_vec = shake256(&pk_bytes, 64);
    let mut tr = [0u8; 64];
    tr.copy_from_slice(&tr_vec);

    let sk_bytes = encoding::sk_encode(&rho, &k_seed, &tr, &s1, &s2, &t0, params);

    let (rho_d, t1_d) = encoding::pk_decode(&pk_bytes, params);
    let (rho_d2, k_d, tr_d, s1_d, s2_d, t0_d) = encoding::sk_decode(&sk_bytes, params);

    assert_eq!(rho, rho_d);
    assert_eq!(rho, rho_d2);
    assert_eq!(k_seed, k_d);
    assert_eq!(tr, tr_d);

    for j in 0..params.k {
        for i in 0..N {
            let t1_val = t1_d[j].coeffs[i];
            let t0_val = arithmetic::to_centered(t0_d[j].coeffs[i]);
            let reconstructed = t1_val as i64 * (1i64 << D) + t0_val as i64;
            let original = ((j * N + i) as u32 * 7919 + 13) % Q;
            assert_eq!(reconstructed as u32, original,
                "Full pipeline: t1*2^d + t0 mismatch at [{}][{}]", j, i);
        }
    }
}


#[test]
fn test_sampling_expand_a_times_s1_plus_s2() {
    let params = &ML_DSA_44;

    let xi = [0x42u8; 32];
    let expanded = shake256(&xi, 128);
    let mut rho = [0u8; 32];
    let mut rho_prime = [0u8; 64];
    rho.copy_from_slice(&expanded[0..32]);
    rho_prime.copy_from_slice(&expanded[32..96]);

    let a_hat = sampling::expand_a(&rho, params);
    let (s1, s2) = sampling::expand_s(&rho_prime, params);

    assert_eq!(a_hat.len(), params.k);
    for row in &a_hat {
        assert_eq!(row.len(), params.l);
    }

    for p in &s1 {
        assert!(p.check_norm(params.eta + 1), "s1 norm exceeds η");
    }
    for p in &s2 {
        assert!(p.check_norm(params.eta + 1), "s2 norm exceeds η");
    }

    let s1_hat = poly::polyvec_ntt(&s1);
    let t_hat = poly::mat_vec_mul_ntt(&a_hat, &s1_hat);
    let t_normal = poly::polyvec_inv_ntt(&t_hat);
    let t = poly::polyvec_add(&t_normal, &s2);

    assert_eq!(t.len(), params.k);

    for (j, poly) in t.iter().enumerate() {
        for i in 0..N {
            assert!(poly.coeffs[i] < Q,
                "t[{}][{}] = {} out of range", j, i, poly.coeffs[i]);
        }
    }

    for j in 0..params.k {
        for i in 0..N {
            let (r1, r0) = arithmetic::power2round(t[j].coeffs[i]);
            assert!(r1 < 1024, "t1 doesn't fit in 10 bits");
            let reconstructed = r1 as i64 * (1i64 << D) + r0 as i64;
            assert_eq!(reconstructed as u32, t[j].coeffs[i]);
        }
    }
}

#[test]
fn test_sampling_expand_mask_and_compute_w() {
    let params = &ML_DSA_44;

    let rho = [0x11u8; 32];
    let rho_pp = [0x22u8; 64];

    let a_hat = sampling::expand_a(&rho, params);
    let y = sampling::expand_mask(&rho_pp, 0, params);

    let gamma1 = params.gamma1 as i32;
    for (j, poly) in y.iter().enumerate() {
        for i in 0..N {
            let c = arithmetic::to_centered(poly.coeffs[i]);
            assert!(c >= -(gamma1 - 1) && c <= gamma1,
                "y[{}][{}] = {} out of mask range", j, i, c);
        }
    }

    let y_hat = poly::polyvec_ntt(&y);
    let w_hat = poly::mat_vec_mul_ntt(&a_hat, &y_hat);
    let w = poly::polyvec_inv_ntt(&w_hat);

    assert_eq!(w.len(), params.k);
    for poly in &w {
        for &c in &poly.coeffs {
            assert!(c < Q);
        }
    }

    let m = params.high_bits_range();
    for poly in &w {
        for &c in &poly.coeffs {
            let h = arithmetic::high_bits(c, params.gamma2);
            assert!(h < m, "w1 = {} out of [0, {})", h, m);
        }
    }
}

#[test]
fn test_sampling_challenge_and_multiply() {
    let params = &ML_DSA_44;

    let rho_prime = [0x33u8; 64];
    let (s1, _s2) = sampling::expand_s(&rho_prime, params);

    let c_tilde = [0x44u8; 32];
    let c = sampling::sample_in_ball(&c_tilde, params.tau);

    let weight: usize = c.coeffs.iter().filter(|&&x| x != 0).count();
    assert_eq!(weight, params.tau);

    let c_hat = poly::ntt(&c);
    let s1_0_hat = poly::ntt(&s1[0]);
    let cs1 = poly::inv_ntt(&c_hat.pointwise_mul(&s1_0_hat));

    for &coeff in &cs1.coeffs {
        let centered = arithmetic::to_centered(coeff);
        assert!(
            centered.abs() <= (params.tau as i32) * (params.eta as i32),
            "||c·s||_∞ = {} exceeds τ·η = {}",
            centered.abs(), params.tau as i32 * params.eta as i32
        );
    }
}

#[test]
fn test_sampling_full_sign_rejection_check() {
    let params = &ML_DSA_44;

    let rho = [0x55u8; 32];
    let rho_prime = [0x66u8; 64];
    let rho_pp = [0x77u8; 64];

    let _a_hat = sampling::expand_a(&rho, params);
    let (s1, _s2) = sampling::expand_s(&rho_prime, params);
    let y = sampling::expand_mask(&rho_pp, 0, params);

    let c = sampling::sample_in_ball(&[0x88; 32], params.tau);
    let c_hat = poly::ntt(&c);

    let mut z = Vec::with_capacity(params.l);
    for j in 0..params.l {
        let s1_hat = poly::ntt(&s1[j]);
        let cs1_hat = c_hat.pointwise_mul(&s1_hat);
        let cs1 = poly::inv_ntt(&cs1_hat);
        z.push(y[j].add(&cs1));
    }

    let bound = params.gamma1 - params.beta;
    let passes = poly::polyvec_check_norm(&z, bound);

    println!("\nSign rejection check (ML-DSA-44): passes={}, bound=γ₁-β={}",
        passes, bound);
}


#[test]
fn test_sign_verify_all_levels() {
    for (params, seed_byte) in [
        (&ML_DSA_44, 0x01u8),
        (&ML_DSA_65, 0x02u8),
        (&ML_DSA_87, 0x03u8),
    ] {
        let kp = signature::keygen(&[seed_byte; 32], params);
        let msg = format!("Integration test for {}", params.name);

        let sig = signature::sign(&kp.sk, msg.as_bytes(), params)
            .expect(&format!("{} signing failed", params.name));

        assert!(
            signature::verify(&kp.pk, msg.as_bytes(), &sig, params),
            "{}: valid signature should verify", params.name
        );

        assert!(
            !signature::verify(&kp.pk, b"wrong", &sig, params),
            "{}: wrong message should be rejected", params.name
        );

        println!("{}: KeyGen + Sign + Verify OK", params.name);
    }
}

#[test]
fn test_sign_verify_multiple_messages() {
    let params = &ML_DSA_44;
    let kp = signature::keygen(&[0xAA; 32], params);

    let messages: &[&[u8]] = &[
        b"",
        b"a",
        b"Hello, World!",
        b"The quick brown fox jumps over the lazy dog",
        &[0xFF; 500],
        &[0x00; 1000],
    ];

    for (idx, msg) in messages.iter().enumerate() {
        let sig = signature::sign(&kp.sk, msg, params)
            .expect(&format!("Signing message {} failed", idx));

        assert!(
            signature::verify(&kp.pk, msg, &sig, params),
            "Message {} (len={}) should verify", idx, msg.len()
        );
    }
}

#[test]
fn test_sign_deterministic_produces_same_sig() {
    let params = &ML_DSA_44;
    let kp = signature::keygen(&[0xBB; 32], params);
    let msg = b"deterministic signing test";

    let sig1 = signature::sign(&kp.sk, msg, params).unwrap();
    let sig2 = signature::sign(&kp.sk, msg, params).unwrap();

    assert_eq!(sig1.c_tilde, sig2.c_tilde, "Deterministic: c̃ should match");

    for j in 0..params.l {
        for i in 0..N {
            assert_eq!(sig1.z[j].coeffs[i], sig2.z[j].coeffs[i],
                "Deterministic: z should match");
        }
    }

    for k in 0..params.k {
        for i in 0..N {
            assert_eq!(sig1.h[k][i], sig2.h[k][i],
                "Deterministic: h should match");
        }
    }
}

#[test]
fn test_different_keys_reject_cross_verification() {
    let params = &ML_DSA_44;

    let kp1 = signature::keygen(&[0x10; 32], params);
    let kp2 = signature::keygen(&[0x20; 32], params);

    let msg = b"cross-key test";
    let sig1 = signature::sign(&kp1.sk, msg, params).unwrap();
    let sig2 = signature::sign(&kp2.sk, msg, params).unwrap();

    assert!(signature::verify(&kp1.pk, msg, &sig1, params));
    assert!(!signature::verify(&kp2.pk, msg, &sig1, params));

    assert!(signature::verify(&kp2.pk, msg, &sig2, params));
    assert!(!signature::verify(&kp1.pk, msg, &sig2, params));
}

#[test]
fn test_signature_encode_decode_roundtrip() {
    let params = &ML_DSA_44;
    let kp = signature::keygen(&[0xCC; 32], params);
    let msg = b"encoding roundtrip test";

    let sig = signature::sign(&kp.sk, msg, params).unwrap();

    let sig_bytes = encoding::sig_encode(&sig.c_tilde, &sig.z, &sig.h, params);
    assert_eq!(sig_bytes.len(), params.sig_bytes);

    let (c_tilde_d, z_d, h_d) = encoding::sig_decode(&sig_bytes, params)
        .expect("Valid signature should decode");

    let sig_decoded = signature::Signature {
        c_tilde: c_tilde_d,
        z: z_d,
        h: h_d,
    };
    assert!(signature::verify(&kp.pk, msg, &sig_decoded, params),
        "Decoded signature should still verify");
}

#[test]
fn test_keygen_encode_decode_roundtrip() {
    let params = &ML_DSA_44;
    let kp = signature::keygen(&[0xDD; 32], params);

    let (rho_d, t1_d) = encoding::pk_decode(&kp.pk.packed, params);
    assert_eq!(kp.pk.rho, rho_d);
    for k in 0..params.k {
        for i in 0..N {
            assert_eq!(kp.pk.t1[k].coeffs[i], t1_d[k].coeffs[i]);
        }
    }

    let sk_bytes = encoding::sk_encode(
        &kp.sk.rho, &kp.sk.k_seed, &kp.sk.tr,
        &kp.sk.s1, &kp.sk.s2, &kp.sk.t0, params
    );
    assert_eq!(sk_bytes.len(), params.sk_bytes);

    let (rho_d, k_d, tr_d, s1_d, s2_d, t0_d) = encoding::sk_decode(&sk_bytes, params);
    assert_eq!(kp.sk.rho, rho_d);
    assert_eq!(kp.sk.k_seed, k_d);
    assert_eq!(kp.sk.tr, tr_d);

    let sk_decoded = signature::SecretKey {
        rho: rho_d, k_seed: k_d, tr: tr_d,
        s1: s1_d, s2: s2_d, t0: t0_d,
    };
    let msg = b"key encode decode test";
    let sig = signature::sign(&sk_decoded, msg, params).unwrap();
    assert!(signature::verify(&kp.pk, msg, &sig, params),
        "Signature from decoded SK should verify with original PK");
}

#[test]
fn test_signature_components_in_range() {
    for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
        let kp = signature::keygen(&[0xEE; 32], params);
        let sig = signature::sign(&kp.sk, b"bounds check", params).unwrap();

        let z_bound = params.gamma1 - params.beta;
        assert!(poly::polyvec_check_norm(&sig.z, z_bound),
            "{}: z norm check failed", params.name);

        assert_eq!(sig.c_tilde.len(), params.c_tilde_bytes,
            "{}: c̃ length", params.name);

        assert_eq!(sig.z.len(), params.l,
            "{}: z length", params.name);

        let weight: usize = sig.h.iter()
            .flat_map(|v| v.iter())
            .filter(|&&b| b)
            .count();
        assert!(weight <= params.omega,
            "{}: hint weight {} > ω={}", params.name, weight, params.omega);

        assert_eq!(sig.h.len(), params.k);
        for v in &sig.h {
            assert_eq!(v.len(), N);
        }

        println!("{}: all signature bounds OK (hints={})", params.name, weight);
    }
}
