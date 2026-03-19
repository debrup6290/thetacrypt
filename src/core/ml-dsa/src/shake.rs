//! SHAKE-128 and SHAKE-256 extendable-output functions (XOF).
//!
//! Built on the Keccak-f[1600] permutation, implemented from scratch
//! per FIPS 202 (SHA-3 Standard).
//!
//! ML-DSA uses SHAKE extensively:
//!   - SHAKE-128: seed expansion for the public matrix A (ExpandA)
//!   - SHAKE-256: seed expansion for secrets, masking, challenge hashing
//!
//! ## Architecture
//!
//! ```text
//! ┌─────────────────────────────────────────────┐
//! │  Keccak-f[1600] permutation (24 rounds)     │
//! │  State: 5×5 × 64-bit lanes = 1600 bits      │
//! └──────────────────┬──────────────────────────┘
//!                    │
//! ┌──────────────────▼──────────────────────────┐
//! │  Sponge construction                         │
//! │  Absorb: XOR input into rate portion → f     │
//! │  Squeeze: read rate portion → f → repeat     │
//! └──────────────────┬──────────────────────────┘
//!                    │
//! ┌──────────────────▼──────────────────────────┐
//! │  SHAKE-128 (rate=168, capacity=32, pad=0x1F) │
//! │  SHAKE-256 (rate=136, capacity=64, pad=0x1F) │
//! └─────────────────────────────────────────────┘
//! ```

//  Keccak-f[1600] Constants 

/// Number of 64-bit lanes in the state: 5 × 5 = 25.
const LANES: usize = 25;

/// Number of rounds for Keccak-f[1600].
const ROUNDS: usize = 24;

/// Round constants for the ι (iota) step.
/// Derived from a degree-8 LFSR; see FIPS 202 §3.2.5.
const RC: [u64; ROUNDS] = [
    0x0000_0000_0000_0001,
    0x0000_0000_0000_8082,
    0x8000_0000_0000_808A,
    0x8000_0000_8000_8000,
    0x0000_0000_0000_808B,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8009,
    0x0000_0000_0000_008A,
    0x0000_0000_0000_0088,
    0x0000_0000_8000_8009,
    0x0000_0000_8000_000A,
    0x0000_0000_8000_808B,
    0x8000_0000_0000_008B,
    0x8000_0000_0000_8089,
    0x8000_0000_0000_8003,
    0x8000_0000_0000_8002,
    0x8000_0000_0000_0080,
    0x0000_0000_0000_800A,
    0x8000_0000_8000_000A,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8080,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8008,
];

/// Rotation offsets for the ρ (rho) step.
/// Indexed as ROTATION[x][y] where state index = x + 5*y.
///
/// Computed via the FIPS 202 algorithm:
///   (x,y) starts at (1,0), offset(t) = (t+1)(t+2)/2 mod 64,
///   then (x,y) ← (y, (2x+3y) mod 5) for each step.
const ROTATION: [[u32; 5]; 5] = [
    //  y=0  y=1  y=2  y=3  y=4
    [    0,  36,   3,  41,  18],  // x=0
    [    1,  44,  10,  45,   2],  // x=1
    [   62,   6,  43,  15,  61],  // x=2
    [   28,  55,  25,  21,  56],  // x=3
    [   27,  20,  39,   8,  14],  // x=4
];

//  Keccak-f[1600] Permutation 

/// The Keccak-f[1600] permutation: 24 rounds of θ, ρ, π, χ, ι.
///
/// State is 25 lanes of 64 bits each, laid out as A[x + 5*y]
/// where x is the column (0..5) and y is the row (0..5).
fn keccak_f(state: &mut [u64; LANES]) {
    for round in 0..ROUNDS {
        // ── θ (theta) ──
        // C[x] = A[x,0] ⊕ A[x,1] ⊕ A[x,2] ⊕ A[x,3] ⊕ A[x,4]
        let mut c = [0u64; 5];
        for x in 0..5 {
            c[x] = state[x]
                ^ state[x + 5]
                ^ state[x + 10]
                ^ state[x + 15]
                ^ state[x + 20];
        }
        // D[x] = C[x-1] ⊕ rot(C[x+1], 1)
        let mut d = [0u64; 5];
        for x in 0..5 {
            d[x] = c[(x + 4) % 5] ^ c[(x + 1) % 5].rotate_left(1);
        }
        // A[x,y] ^= D[x]
        for x in 0..5 {
            for y in 0..5 {
                state[x + 5 * y] ^= d[x];
            }
        }

        // ── ρ (rho) and π (pi), combined ──
        // B[y, 2x+3y mod 5] = rot(A[x,y], ρ[x][y])
        let mut b = [0u64; LANES];
        for x in 0..5 {
            for y in 0..5 {
                let new_x = y;
                let new_y = (2 * x + 3 * y) % 5;
                b[new_x + 5 * new_y] =
                    state[x + 5 * y].rotate_left(ROTATION[x][y]);
            }
        }

        // ── χ (chi) ──
        // A[x,y] = B[x,y] ⊕ (¬B[x+1,y] ∧ B[x+2,y])
        for x in 0..5 {
            for y in 0..5 {
                state[x + 5 * y] = b[x + 5 * y]
                    ^ (!b[(x + 1) % 5 + 5 * y] & b[(x + 2) % 5 + 5 * y]);
            }
        }

        // ── ι (iota) ──
        state[0] ^= RC[round];
    }
}

//  Sponge Construction 

/// Internal state of the Keccak sponge.
///
/// Supports two phases:
///   1. **Absorb**: input bytes are XORed into the rate portion of the state.
///      When the rate is full, Keccak-f is applied and the position resets.
///   2. **Squeeze**: after finalization (pad + Keccak-f), output bytes are
///      read from the rate portion, applying Keccak-f for each new block.
///
/// Phase transition is one-way: once finalized for squeezing, no more
/// absorbing is allowed.
struct Sponge {
    /// The 1600-bit Keccak state, as 25 × 64-bit lanes (little-endian byte order).
    state: [u64; LANES],

    /// Rate in bytes. Determines how many bytes of the state are used for
    /// absorb/squeeze per block.
    ///   SHAKE-128: rate = 168  (capacity = 32 bytes = 256 bits)
    ///   SHAKE-256: rate = 136  (capacity = 64 bytes = 512 bits)
    rate: usize,

    /// Current byte offset within the rate portion.
    /// During absorb: how many bytes have been XORed since the last permutation.
    /// During squeeze: how many bytes have been read since the last permutation.
    offset: usize,

    /// Whether finalize() has been called.
    finalized: bool,
}

impl Sponge {
    /// Create a new sponge with the given rate (in bytes).
    fn new(rate: usize) -> Self {
        debug_assert!(rate <= 200, "Rate cannot exceed state size (200 bytes)");
        debug_assert!(rate > 0, "Rate must be positive");
        Self {
            state: [0u64; LANES],
            rate,
            offset: 0,
            finalized: false,
        }
    }

    /// XOR a single byte into the state at the current position.
    /// State is little-endian: byte `i` maps to lane `i/8`, shift `(i%8)*8`.
    #[inline(always)]
    fn xor_byte(&mut self, pos: usize, byte: u8) {
        let lane = pos / 8;
        let shift = (pos % 8) * 8;
        self.state[lane] ^= (byte as u64) << shift;
    }

    /// Read a single byte from the state at the given position.
    #[inline(always)]
    fn read_byte(&self, pos: usize) -> u8 {
        let lane = pos / 8;
        let shift = (pos % 8) * 8;
        ((self.state[lane] >> shift) & 0xFF) as u8
    }

    /// Absorb a chunk of input bytes into the sponge.
    ///
    /// Can be called multiple times before finalization.
    /// Panics if called after squeezing has started.
    fn absorb(&mut self, data: &[u8]) {
        assert!(!self.finalized, "Cannot absorb after finalization");

        for &byte in data {
            self.xor_byte(self.offset, byte);
            self.offset += 1;

            if self.offset == self.rate {
                keccak_f(&mut self.state);
                self.offset = 0;
            }
        }
    }

    /// Finalize the absorb phase with SHAKE domain separation padding.
    ///
    /// SHAKE padding (FIPS 202 §6.2):
    ///   - XOR 0x1F at current position (domain separator for SHAKE: 1111 in binary,
    ///     which is 0x0F for the 4-bit suffix, then the multi-rate pad adds 1 → 0x1F)
    ///   - XOR 0x80 at the last byte of the rate block
    ///   - Apply Keccak-f
    fn finalize(&mut self) {
        if self.finalized {
            return;
        }

        // SHAKE suffix: 0x1F (= 4-bit domain separator 1111 + start of pad10*1)
        self.xor_byte(self.offset, 0x1F);

        // Final bit of pad10*1: set the high bit of the last rate byte
        self.xor_byte(self.rate - 1, 0x80);

        // Apply the permutation
        keccak_f(&mut self.state);
        self.offset = 0;
        self.finalized = true;
    }

    /// Squeeze output bytes from the sponge.
    ///
    /// Automatically calls finalize() on first squeeze.
    /// Can be called repeatedly to generate arbitrarily many output bytes.
    fn squeeze(&mut self, output: &mut [u8]) {
        if !self.finalized {
            self.finalize();
        }

        for byte in output.iter_mut() {
            if self.offset == self.rate {
                keccak_f(&mut self.state);
                self.offset = 0;
            }
            *byte = self.read_byte(self.offset);
            self.offset += 1;
        }
    }
}

//  Public API 

/// SHAKE-128 extendable-output function.
///
/// Security: 128-bit collision resistance, 128-bit preimage resistance.
/// Rate: 168 bytes (1344 bits).  Capacity: 32 bytes (256 bits).
///
/// Used in ML-DSA for expanding the matrix A (ExpandA via RejNTTPoly).
pub struct Shake128 {
    sponge: Sponge,
}

impl Shake128 {
    /// Create a new, empty SHAKE-128 instance.
    pub fn new() -> Self {
        Self { sponge: Sponge::new(168) }
    }

    /// Create a SHAKE-128 instance and absorb initial data.
    pub fn from_data(data: &[u8]) -> Self {
        let mut s = Self::new();
        s.absorb(data);
        s
    }

    /// Absorb input bytes. Can be called multiple times.
    pub fn absorb(&mut self, data: &[u8]) {
        self.sponge.absorb(data);
    }

    /// Squeeze output bytes into the provided buffer.
    pub fn squeeze(&mut self, output: &mut [u8]) {
        self.sponge.squeeze(output);
    }

    /// Convenience: squeeze `n` bytes and return as a Vec.
    pub fn squeeze_vec(&mut self, n: usize) -> Vec<u8> {
        let mut out = vec![0u8; n];
        self.sponge.squeeze(&mut out);
        out
    }
}

/// SHAKE-256 extendable-output function.
///
/// Security: 256-bit collision resistance, 256-bit preimage resistance.
/// Rate: 136 bytes (1088 bits).  Capacity: 64 bytes (512 bits).
///
/// Used in ML-DSA for secret expansion (ExpandS), masking (ExpandMask),
/// challenge generation, and general hashing (H, H').
pub struct Shake256 {
    sponge: Sponge,
}

impl Shake256 {
    /// Create a new, empty SHAKE-256 instance.
    pub fn new() -> Self {
        Self { sponge: Sponge::new(136) }
    }

    /// Create a SHAKE-256 instance and absorb initial data.
    pub fn from_data(data: &[u8]) -> Self {
        let mut s = Self::new();
        s.absorb(data);
        s
    }

    /// Absorb input bytes. Can be called multiple times.
    pub fn absorb(&mut self, data: &[u8]) {
        self.sponge.absorb(data);
    }

    /// Squeeze output bytes into the provided buffer.
    pub fn squeeze(&mut self, output: &mut [u8]) {
        self.sponge.squeeze(output);
    }

    /// Convenience: squeeze `n` bytes and return as a Vec.
    pub fn squeeze_vec(&mut self, n: usize) -> Vec<u8> {
        let mut out = vec![0u8; n];
        self.sponge.squeeze(&mut out);
        out
    }
}

//  Convenience Functions 

/// One-shot SHAKE-128: hash `input` and return `outlen` bytes.
pub fn shake128(input: &[u8], outlen: usize) -> Vec<u8> {
    let mut xof = Shake128::from_data(input);
    xof.squeeze_vec(outlen)
}

/// One-shot SHAKE-256: hash `input` and return `outlen` bytes.
pub fn shake256(input: &[u8], outlen: usize) -> Vec<u8> {
    let mut xof = Shake256::from_data(input);
    xof.squeeze_vec(outlen)
}

//  Tests 

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: convert hex string to bytes.
    fn hex(s: &str) -> Vec<u8> {
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
            .collect()
    }

    /// Helper: convert bytes to hex string (for diagnostics).
    fn to_hex(bytes: &[u8]) -> String {
        bytes.iter().map(|b| format!("{:02x}", b)).collect()
    }

    // ── Keccak-f[1600] permutation tests ──

    #[test]
    fn test_keccak_f_zero_state_lane0() {
        // The most widely cited Keccak-f[1600] zero-state test:
        // After applying keccak_f to all zeros, lane[0] must be this value.
        // Source: XKCP (eXtended Keccak Code Package) test vectors.
        let mut state = [0u64; LANES];
        keccak_f(&mut state);
        assert_eq!(
            state[0], 0xF1258F7940E1DDE7,
            "Keccak-f zero-state lane[0] mismatch"
        );
    }

    #[test]
    fn test_keccak_f_deterministic() {
        // Applying keccak_f twice to the same input must give the same result.
        let mut state1 = [0u64; LANES];
        let mut state2 = [0u64; LANES];
        keccak_f(&mut state1);
        keccak_f(&mut state2);
        assert_eq!(state1, state2);
    }

    #[test]
    fn test_keccak_f_not_identity() {
        // keccak_f on non-zero state should change the state.
        let mut state = [0u64; LANES];
        state[0] = 1;
        let original = state;
        keccak_f(&mut state);
        assert_ne!(state, original, "keccak_f should not be the identity");
    }

    #[test]
    fn test_keccak_f_zero_state_print_lanes() {
        // Print all 25 lanes for manual verification against XKCP reference.
        // Run with: cargo test test_keccak_f_zero_state_print -- --nocapture
        let mut state = [0u64; LANES];
        keccak_f(&mut state);
        println!("\n=== Keccak-f[1600] on zero state ===");
        for i in 0..LANES {
            println!("  lane[{:2}] = 0x{:016X}", i, state[i]);
        }
        // Verify the first two lanes (universally agreed upon):
        assert_eq!(state[0], 0xF1258F7940E1DDE7, "lane[0]");
        assert_eq!(state[1], 0x84D5CCF933C0478A, "lane[1]");
    }

    // ── SHAKE-256 Known Answer Tests ──
    // Source: NIST CAVP (Cryptographic Algorithm Validation Program)

    #[test]
    fn test_shake256_empty_32() {
        // SHAKE-256("", 32) — most authoritative SHAKE test vector.
        // Source: NIST ShortMsgKAT_SHAKE256.txt, Len=0
        let expected = hex(
            "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
        );
        let result = shake256(&[], 32);
        assert_eq!(result, expected,
            "\nSHAKE-256(\"\", 32):\n  got:    {}\n  expect: {}",
            to_hex(&result), to_hex(&expected)
        );
    }

    #[test]
    fn test_shake256_abc_32() {
        // SHAKE-256("abc", 32)
        // Source: NIST ShortMsgKAT_SHAKE256.txt, Len=24 bits
        let expected = hex(
            "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739"
        );
        let result = shake256(b"abc", 32);
        assert_eq!(result, expected,
            "\nSHAKE-256(\"abc\", 32):\n  got:    {}\n  expect: {}",
            to_hex(&result), to_hex(&expected)
        );
    }

    // ── SHAKE-128 Known Answer Tests ──

    #[test]
    fn test_shake128_empty_32() {
        // SHAKE-128("", 32) — authoritative SHAKE-128 test vector.
        // Source: NIST ShortMsgKAT_SHAKE128.txt, Len=0
        let expected = hex(
            "7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26"
        );
        let result = shake128(&[], 32);
        assert_eq!(result, expected,
            "\nSHAKE-128(\"\", 32):\n  got:    {}\n  expect: {}",
            to_hex(&result), to_hex(&expected)
        );
    }

    // ── Extended output: self-consistency tests ──
    // (We verify correctness by checking that extended output is consistent
    //  with the verified first-32-byte KAT, not by hard-coding longer vectors.)

    #[test]
    fn test_shake256_extended_output_consistent() {
        // Squeeze 64 bytes at once. The first 32 must match the verified KAT.
        let result_64 = shake256(&[], 64);
        let result_32 = shake256(&[], 32);

        assert_eq!(
            &result_64[..32], &result_32[..],
            "First 32 bytes of 64-byte squeeze must match 32-byte squeeze"
        );
        // Bytes 32-63 must be non-trivial
        assert!(
            result_64[32..].iter().any(|&b| b != 0),
            "Extended output should not be all zeros"
        );
    }

    #[test]
    fn test_shake128_extended_output_consistent() {
        let result_64 = shake128(&[], 64);
        let result_32 = shake128(&[], 32);

        assert_eq!(
            &result_64[..32], &result_32[..],
            "First 32 bytes of 64-byte squeeze must match 32-byte squeeze"
        );
    }

    #[test]
    fn test_shake256_extended_output_print() {
        // Print extended output for manual verification against another implementation.
        // Run with: cargo test test_shake256_extended -- --nocapture
        let result = shake256(&[], 64);
        println!("\n=== SHAKE-256(\"\", 64) ===");
        println!("  {}", to_hex(&result));
    }

    #[test]
    fn test_shake128_abc_print() {
        // Print SHAKE-128("abc") output for manual verification.
        // Run with: cargo test test_shake128_abc -- --nocapture
        let result = shake128(b"abc", 32);
        println!("\n=== SHAKE-128(\"abc\", 32) ===");
        println!("  {}", to_hex(&result));
    }

    // ── Incremental absorb consistency ──

    #[test]
    fn test_incremental_absorb_shake256() {
        // Absorbing in chunks must equal absorbing all at once.
        let data = b"The quick brown fox jumps over the lazy dog";

        // One-shot
        let oneshot = shake256(data, 64);

        // Byte-by-byte
        let mut xof = Shake256::new();
        for &byte in data.iter() {
            xof.absorb(&[byte]);
        }
        let bytewise = xof.squeeze_vec(64);
        assert_eq!(oneshot, bytewise, "Byte-by-byte absorb must match one-shot");

        // Chunked (various sizes)
        let mut xof2 = Shake256::new();
        xof2.absorb(&data[..7]);
        xof2.absorb(&data[7..20]);
        xof2.absorb(&data[20..]);
        let chunked = xof2.squeeze_vec(64);
        assert_eq!(oneshot, chunked, "Chunked absorb must match one-shot");
    }

    #[test]
    fn test_incremental_absorb_shake128() {
        let data = b"SHAKE-128 incremental test";

        let oneshot = shake128(data, 48);

        let mut xof = Shake128::new();
        xof.absorb(&data[..10]);
        xof.absorb(&data[10..]);
        let chunked = xof.squeeze_vec(48);

        assert_eq!(oneshot, chunked);
    }

    // ── Incremental squeeze consistency ──

    #[test]
    fn test_incremental_squeeze_shake256() {
        // Squeezing in chunks must produce the same stream as squeezing all at once.
        let data = b"squeeze test";

        // Squeeze 200 bytes at once (crosses the rate boundary of 136)
        let mut xof1 = Shake256::from_data(data);
        let all_at_once = xof1.squeeze_vec(200);

        // Squeeze 200 bytes in small chunks
        let mut xof2 = Shake256::from_data(data);
        let mut piecewise = Vec::new();
        for chunk_size in &[10, 50, 3, 70, 1, 66] {
            let chunk = xof2.squeeze_vec(*chunk_size);
            piecewise.extend_from_slice(&chunk);
        }

        assert_eq!(all_at_once, piecewise, "Incremental squeeze must match bulk squeeze");
    }

    #[test]
    fn test_incremental_squeeze_shake128() {
        let data = b"shake128 squeeze";

        let mut xof1 = Shake128::from_data(data);
        let all_at_once = xof1.squeeze_vec(300); // crosses 168 rate boundary

        let mut xof2 = Shake128::from_data(data);
        let mut piecewise = Vec::new();
        for chunk_size in &[50, 100, 20, 130] {
            let chunk = xof2.squeeze_vec(*chunk_size);
            piecewise.extend_from_slice(&chunk);
        }

        assert_eq!(all_at_once, piecewise);
    }

    // ── Large input / output ──

    #[test]
    fn test_large_absorb() {
        // Absorb more than one rate block (168 bytes for SHAKE-128)
        let data = vec![0xAB; 500];

        let result1 = shake128(&data, 32);

        let mut xof = Shake128::new();
        xof.absorb(&data);
        let result2 = xof.squeeze_vec(32);

        assert_eq!(result1, result2);
        assert_ne!(result1, vec![0u8; 32], "Output should not be all zeros");
    }

    #[test]
    fn test_squeeze_across_many_blocks() {
        // Squeeze enough output to trigger multiple Keccak-f invocations.
        // SHAKE-256 rate = 136, so 1000 bytes = ~7.4 blocks.
        let mut xof = Shake256::from_data(b"multi-block squeeze");
        let out = xof.squeeze_vec(1000);

        // Verify determinism
        let mut xof2 = Shake256::from_data(b"multi-block squeeze");
        let out2 = xof2.squeeze_vec(1000);
        assert_eq!(out, out2);

        // Sanity: output shouldn't be all zeros
        assert!(out.iter().any(|&b| b != 0));
    }

    // ── Basic properties ──

    #[test]
    fn test_deterministic() {
        for _ in 0..5 {
            let a = shake256(b"determinism", 32);
            let b = shake256(b"determinism", 32);
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_different_inputs() {
        let a = shake256(b"input A", 32);
        let b = shake256(b"input B", 32);
        assert_ne!(a, b, "Different inputs should produce different hashes");
    }

    // ── Rate boundary edge cases ──

    #[test]
    fn test_absorb_exactly_one_rate_block() {
        // Absorb exactly `rate` bytes — triggers exactly one keccak_f during absorb.
        let data_128 = vec![0x42; 168]; // SHAKE-128 rate
        let result = shake128(&data_128, 32);
        assert_ne!(result, vec![0u8; 32]);

        let data_256 = vec![0x42; 136]; // SHAKE-256 rate
        let result2 = shake256(&data_256, 32);
        assert_ne!(result2, vec![0u8; 32]);
    }

    #[test]
    fn test_absorb_rate_plus_one() {
        // Absorb rate+1 bytes — exercises the boundary precisely.
        let data = vec![0xFF; 169]; // 168 + 1 for SHAKE-128
        let r1 = shake128(&data, 32);
        assert_ne!(r1, vec![0u8; 32]);

        let data2 = vec![0xFF; 137]; // 136 + 1 for SHAKE-256
        let r2 = shake256(&data2, 32);
        assert_ne!(r2, vec![0u8; 32]);
    }
}