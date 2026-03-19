use ml_dsa::{params, signature, poly};
use std::time::Instant;

fn main() {
    println!("\n");
    println!("ML-DSA (FIPS 204) — From-Scratch Implementation");
    println!("\n");

    let levels: &[&params::Params] = &[
        &params::ML_DSA_44,
        &params::ML_DSA_65,
        &params::ML_DSA_87,
    ];

    for p in levels {
        println!("  {}", p.name);
        println!("  Matrix: {} x {},  eta: {},  tau: {},  gamma1: 2^{}",
            p.k, p.l, p.eta, p.tau,
            if p.gamma1 == (1 << 17) { 17 } else { 19 }
        );
        println!("  PK: {} bytes,  SK: {} bytes,  Sig: {} bytes",
            p.pk_bytes, p.sk_bytes, p.sig_bytes);

        // KeyGen
        let seed = [0x42u8; 32];
        let t0 = Instant::now();
        let kp = signature::keygen(&seed, p);
        let keygen_ms = t0.elapsed().as_millis();
        println!("  KeyGen:  OK  ({} ms,  pk = {} bytes)", keygen_ms, kp.pk.packed.len());

        // Sign
        let msg = b"Hello, post-quantum world! This is a test of ML-DSA.";
        let t0 = Instant::now();
        let sig = signature::sign(&kp.sk, msg, p).expect("Signing should succeed");
        let sign_ms = t0.elapsed().as_millis();

        let hint_count: usize = sig.h.iter()
            .flat_map(|v| v.iter())
            .filter(|&&b| b)
            .count();
        println!("  Sign:    OK  ({} ms,  hints = {})", sign_ms, hint_count);

        // Verify (correct)
        let t0 = Instant::now();
        let valid = signature::verify(&kp.pk, msg, &sig, p);
        let verify_ms = t0.elapsed().as_millis();
        println!("  Verify:  {}  ({} ms)",
            if valid { "PASS" } else { "FAIL" }, verify_ms);

        // Verify (tampered message)
        let tampered = signature::verify(&kp.pk, b"tampered!", &sig, p);
        println!("  Verify (tampered msg):  {}",
            if !tampered { "correctly rejected" } else { "FAIL - should reject!" });

        // Verify (wrong key)
        let kp2 = signature::keygen(&[0x99; 32], p);
        let wrong_key = signature::verify(&kp2.pk, msg, &sig, p);
        println!("  Verify (wrong key):     {}",
            if !wrong_key { "correctly rejected" } else { "FAIL - should reject!" });

        // z norm check
        let bound = p.gamma1 - p.beta;
        let z_ok = poly::polyvec_check_norm(&sig.z, bound);
        println!("  z norm < gamma1-beta:   {}", if z_ok { "OK" } else { "FAIL" });

        println!();
    }

    println!("All done! ML-DSA implementation complete.");
}