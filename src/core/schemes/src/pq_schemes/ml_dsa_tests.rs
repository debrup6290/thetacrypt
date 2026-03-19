#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use theta_proto::scheme_types::{Group, ThresholdScheme};
    use crate::keys::key_generator::KeyGenerator;
    use crate::rand::{RngAlgorithm, RNG};
    use crate::interface::{Signature, ThresholdSignature};
    use crate::pq_schemes::ml_dsa::MlDsaSignature;

    // Import protocol types — these live in theta_protocols but we test
    // the crypto layer here directly via the ml_dsa crate.
    use ml_dsa::threshold::{rss, params as thresh_params, protocol, combine};

    #[test]
    fn test_ml_dsa_44_keygen_produces_correct_count() {
        let mut rng = RNG::new(RngAlgorithm::OsRng);
        let keys = KeyGenerator::generate_keys(
            2, 3, &mut rng,
            &ThresholdScheme::MlDsa44,
            &Group::Lattice,
            &None,
        ).expect("keygen failed");
        assert_eq!(keys.len(), 3, "should produce one key per party");
    }

    #[test]
    fn test_ml_dsa_65_threshold_sign_and_verify_2_of_3() {
        // Generate 2-of-3 ML-DSA-65 keys
        let mut rng = RNG::new(RngAlgorithm::OsRng);
        let party_keys = KeyGenerator::generate_keys(
            2, 3, &mut rng,
            &ThresholdScheme::MlDsa65,
            &Group::Lattice,
            &None,
        ).expect("keygen failed");

        let pub_key = party_keys[0].get_public_key();

        // Extract the underlying ml_dsa types from party keys 0 and 1
        let pk0 = match &party_keys[0] {
            crate::keys::keys::PrivateKeyShare::MlDsa65(k) => k,
            _ => panic!("wrong key type"),
        };
        let pk1 = match &party_keys[1] {
            crate::keys::keys::PrivateKeyShare::MlDsa65(k) => k,
            _ => panic!("wrong key type"),
        };

        let tp = pk0.to_threshold_params().expect("threshold params");
        let vk = pk0.to_verification_key().expect("vk");
        let party_key0 = pk0.to_party_key().expect("party key 0");
        let party_key1 = pk1.to_party_key().expect("party key 1");

        let msg = b"post-quantum threshold signing test";
        let act = vec![0usize, 1usize];

        // Attempt up to 10 times (probabilistic acceptance)
        let mut sig_bytes: Option<Vec<u8>> = None;
        for attempt in 0u8..10 {
            let seed0 = vec![attempt; 32];
            let seed1 = vec![attempt + 128; 32];

            // Round 1
            let (mut state0, r1_0) = protocol::share_sign_1(&vk, &party_key0, &tp, &seed0);
            let (mut state1, r1_1) = protocol::share_sign_1(&vk, &party_key1, &tp, &seed1);

            // Round 2
            let r2_0 = protocol::share_sign_2(
                &mut state0, &vk, &act, msg,
                vec![r1_0.clone_msg(), r1_1.clone_msg()],
            ).expect("round 2 party 0");
            let r2_1 = protocol::share_sign_2(
                &mut state1, &vk, &act, msg,
                vec![r1_0.clone_msg(), r1_1.clone_msg()],
            ).expect("round 2 party 1");

            // Round 3
            let r3_0 = protocol::share_sign_3(
                &state0, &vk, &party_key0,
                &[r2_0.clone_msg(), r2_1.clone_msg()],
            ).expect("round 3 party 0");
            let r3_1 = protocol::share_sign_3(
                &state1, &vk, &party_key1,
                &[r2_0.clone_msg(), r2_1.clone_msg()],
            ).expect("round 3 party 1");

            // Combine
            if let Some(sig) = combine::combine(
                &vk, &act, msg,
                &[r2_0.clone_msg(), r2_1.clone_msg()],
                &[r3_0.clone_msg(), r3_1.clone_msg()],
                &tp,
            ) {
                sig_bytes = Some(sig.to_bytes(tp.base));
                break;
            }
        }

        let sig_bytes = sig_bytes.expect("signing failed after 10 attempts");

        // Wrap in Thetacrypt Signature enum and verify through unified interface
        let sig = Signature::MlDsa65(MlDsaSignature { bytes: sig_bytes });
        let valid = ThresholdSignature::verify(&sig, &pub_key, msg)
            .expect("verify returned error");
        assert!(valid, "signature should verify");

        // Confirm tampered message fails
        let invalid = ThresholdSignature::verify(&sig, &pub_key, b"tampered")
            .expect("verify returned error");
        assert!(!invalid, "tampered message should not verify");
    }
}
