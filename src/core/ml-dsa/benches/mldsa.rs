use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};

use ml_dsa::params::{Params, ML_DSA_44, ML_DSA_65, ML_DSA_87};
use ml_dsa::signature;

fn bench_suite(c: &mut Criterion, params: &'static Params, label: &'static str) {
    let msg: Vec<u8> = vec![0u8; 32];

    // KeyGen benchmark: vary seed each iteration
    c.bench_function(&format!("keygen_{}", label), |b| {
        let mut ctr: u64 = 0;
        b.iter(|| {
            ctr = ctr.wrapping_add(1);
            let mut seed = [0u8; 32];
            seed[0..8].copy_from_slice(&ctr.to_le_bytes());
            let kp = signature::keygen(black_box(&seed), black_box(params));
            black_box(kp);
        })
    });

    // Setup once for sign/verify benches
    let seed = [0x42u8; 32];
    let kp = signature::keygen(&seed, params);

    c.bench_function(&format!("sign_{}", label), |b| {
        b.iter(|| {
            let sig = signature::sign(black_box(&kp.sk), black_box(&msg), black_box(params))
                .expect("sign failed");
            black_box(sig);
        })
    });

    let sig = signature::sign(&kp.sk, &msg, params).expect("sign failed");

    c.bench_function(&format!("verify_{}", label), |b| {
        b.iter(|| {
            let ok = signature::verify(black_box(&kp.pk), black_box(&msg), black_box(&sig), black_box(params));
            black_box(ok);
        })
    });

    // Batch verify throughput (128 signatures)
    c.bench_function(&format!("verify_batch128_{}", label), |b| {
        b.iter_batched(
            || {
                let mut sigs = Vec::with_capacity(128);
                for i in 0..128u64 {
                    let mut m = msg.clone();
                    m[0..8].copy_from_slice(&i.to_le_bytes());
                    sigs.push(signature::sign(&kp.sk, &m, params).expect("sign failed"));
                }
                sigs
            },
            |sigs| {
                let mut ok = true;
                for (i, s) in sigs.iter().enumerate() {
                    let mut m = msg.clone();
                    m[0..8].copy_from_slice(&(i as u64).to_le_bytes());
                    ok &= signature::verify(&kp.pk, &m, s, params);
                }
                black_box(ok);
            },
            BatchSize::SmallInput,
        )
    });
}

fn mldsa_benches(c: &mut Criterion) {
    bench_suite(c, &ML_DSA_44, "44");
    bench_suite(c, &ML_DSA_65, "65");
    bench_suite(c, &ML_DSA_87, "87");
}

criterion_group!(benches, mldsa_benches);
criterion_main!(benches);