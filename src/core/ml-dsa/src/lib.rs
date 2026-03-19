//! # ML-DSA (FIPS 204) — Module-Lattice Digital Signature Algorithm
//!
//! A from-scratch Rust implementation of ML-DSA (formerly CRYSTALS-Dilithium).
//!
//! ## Modules
//! - `params`    — Security parameter sets (ML-DSA-44, 65, 87)
//! - `arithmetic`— Modular arithmetic mod q = 8380417
//! - `poly`      — Polynomial ring Rq = Zq[X]/(X^256+1) and NTT
//! - `encoding`  — BitPack / BitUnpack / byte encoding per FIPS 204
//! - `shake`     — SHAKE-128/256 (Keccak) XOF
//! - `sampling`  — ExpandA, ExpandS, ExpandMask, SampleInBall
//! - `signature` — KeyGen, Sign, Verify

pub mod params;
pub mod shake;
pub mod arithmetic;
pub mod poly;
pub mod encoding;
pub mod sampling;
pub mod signature;
pub mod threshold;