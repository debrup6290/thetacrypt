//! Threshold ML-DSA — Distributed signing per PQShield (CCS '25 / USENIX Security '26).
//!
//! Implements the threshold variant of ML-DSA where T-of-N parties cooperatively
//! produce a **standard** FIPS 204 signature, verifiable by any existing ML-DSA
//! implementation with no verifier modifications.
//!
//! ## Architecture
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────┐
//! │  combine.rs   — Assemble partial sigs → standard ML-DSA sig │
//! ├─────────────────────────────────────────────────────────────┤
//! │  protocol.rs  — ShareSign 3-round state machine             │
//! │                 Round 1: Commit (hash of w_i)               │
//! │                 Round 2: Reveal (open w_i, verify commits)  │
//! │                 Round 3: Respond (partial sig z_i via HRej) │
//! ├─────────────────────────────────────────────────────────────┤
//! │  hyperball.rs — Gaussian / hyperball sampling, HRej          │
//! │  rss.rs       — Replicated Short Secret Sharing + RSSRecover │
//! │  params.rs    — Threshold parameters (Tables 2–3 from paper) │
//! └─────────────────────────────────────────────────────────────┘
//! ```
//!
//! ## Protocol Flow
//!
//! 1. **KeyGen** (`rss.rs`): Generate RSS shares of the secret key.
//!    Each party receives C(N, N-T+1) short secret vectors.
//!
//! 2. **ShareSign** (`protocol.rs`): Three-round interactive signing:
//!    - Round 1: Each party commits to its masking vector w_i.
//!    - Round 2: Parties reveal w_i and verify commitments.
//!    - Round 3: Each party computes a partial signature z_i using
//!      hyperball rejection sampling (HRej).
//!
//! 3. **Combine** (`combine.rs`): A combiner (any party or external)
//!    aggregates partial signatures into a standard ML-DSA signature.
//!
//! 4. **Verify**: Uses the unchanged `signature::verify()` from the
//!    base ML-DSA implementation.
//!
//! ## References
//!
//! - PQShield paper: <https://eprint.iacr.org/2026/013>
//! - Full version: <https://inria.hal.science/hal-05442192v1/document>

pub mod params;
pub mod rss;
pub mod hyperball;
pub mod protocol;
pub mod combine;