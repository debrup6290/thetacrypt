# ml-dsa

From-scratch Rust implementation of FIPS 204 ML-DSA and the threshold extension
from [Efficient Threshold ML-DSA](https://eprint.iacr.org/2026/013).

## What's here

- `src/` — standard ML-DSA: arithmetic, sampling, encoding, signing, verification
- `src/threshold/` — threshold protocol: RSS key sharing, 3-round signing protocol, combine

## Usage

Generate keys and run a 2-of-3 signing session:

    cargo test -p theta_schemes ml_dsa

## Security levels

| Variant   | NIST Level | λ   |
|-----------|-----------|-----|
| ML-DSA-44 | 2         | 128 |
| ML-DSA-65 | 3         | 192 |
| ML-DSA-87 | 5         | 256 |

## Notes

- Signatures produced are standard FIPS 204 and verify with any unmodified ML-DSA verifier
- Threshold keygen uses Replicated Short Secret Sharing (RSS) over the secret key polynomials
- Signing is probabilistic — the protocol retries internally across K parallel instances