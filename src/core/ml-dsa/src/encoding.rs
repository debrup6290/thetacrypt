//! Byte encoding and decoding of polynomials per FIPS 204.
//!
//! ML-DSA needs to serialize polynomials with different coefficient ranges
//! into compact byte strings.  This module implements:
//!
//! **Bit packing** (FIPS 204, Algorithms 16–19):
//!   - `SimpleBitPack`:   pack coefficients in [0, b]
//!   - `BitPack`:         pack coefficients in [a, b] (signed range)
//!   - `SimpleBitUnpack`: unpack to unsigned range
//!   - `BitUnpack`:       unpack to signed range
//!
//! **Key and signature encoding** (FIPS 204, §§7.1–7.3):
//!   - `pk_encode` / `pk_decode`:   public key = ρ || t1
//!   - `sk_encode` / `sk_decode`:   secret key = ρ || K || tr || s1 || s2 || t0
//!   - `sig_encode` / `sig_decode`: signature = c̃ || z || h
//!
//! **Hint encoding** (FIPS 204, Algorithms 22–23):
//!   - `hint_bit_pack` / `hint_bit_unpack`: sparse encoding of the hint vector

use crate::arithmetic;
use crate::params::{N, Q, D, Params};
use crate::poly::Poly;

//  Bit Packing 

/// SimpleBitPack: pack polynomial coefficients in [0, 2^bits - 1] using
/// `bits` bits per coefficient, little-endian byte order.
///
/// FIPS 204, Algorithm 16.
///
/// Output length: ceil(256 * bits / 8) bytes.
pub fn simple_bit_pack(p: &Poly, bits: u32) -> Vec<u8> {
    let total_bits = N * bits as usize;
    let total_bytes = (total_bits + 7) / 8;
    let mut buf = vec![0u8; total_bytes];

    let mut bit_pos: usize = 0;
    for i in 0..N {
        let val = p.coeffs[i];
        // Write `bits` bits of val starting at bit_pos in buf
        write_bits(&mut buf, bit_pos, val as u64, bits);
        bit_pos += bits as usize;
    }

    buf
}

/// SimpleBitUnpack: unpack bytes into polynomial with coefficients
/// in [0, 2^bits - 1].
///
/// FIPS 204, Algorithm 17.
pub fn simple_bit_unpack(data: &[u8], bits: u32) -> Poly {
    let mut p = Poly::zero();
    let mask = (1u64 << bits) - 1;

    let mut bit_pos: usize = 0;
    for i in 0..N {
        p.coeffs[i] = (read_bits(data, bit_pos, bits) & mask) as u32;
        bit_pos += bits as usize;
    }

    p
}

/// BitPack: pack polynomial coefficients in [a, b] by mapping to [0, a+b]
/// via the transformation: packed = a - coeff (where coeff is centered).
///
/// FIPS 204, Algorithm 18.
///
/// The coefficients are stored as standard representatives in [0, q),
/// but conceptually represent centered values in [-a, b].  We shift
/// them to [0, a+b] for packing.
///
/// Parameters:
///   - `a`: the magnitude of the most negative value (positive number)
///   - `b`: the magnitude of the most positive value (positive number)
///
/// Each coefficient c (centered in [-a, b]) is packed as (b - c) in [0, a+b].
pub fn bit_pack(p: &Poly, a: u32, b: u32) -> Vec<u8> {
    let range = a + b;
    let bits = bit_length(range);

    let mut shifted = Poly::zero();
    for i in 0..N {
        let centered = arithmetic::to_centered(p.coeffs[i]);
        // FIPS 204 Algorithm 17: packed_value = b - centered
        // Maps centered in [-a, b] to unsigned in [0, a+b]:
        //   b → b - b = 0
        //   -a → b - (-a) = a + b
        shifted.coeffs[i] = (b as i32 - centered) as u32;
    }

    simple_bit_pack(&shifted, bits)
}

/// BitUnpack: unpack bytes to polynomial with centered coefficients in [-a, b].
///
/// FIPS 204, Algorithm 19.
///
/// Reverses BitPack: reads unsigned values in [0, a+b], then maps back
/// to centered via FIPS 204: centered = b - packed_value.
///
/// Returns coefficients as standard representatives in [0, q).
pub fn bit_unpack(data: &[u8], a: u32, b: u32) -> Poly {
    let range = a + b;
    let bits = bit_length(range);

    let shifted = simple_bit_unpack(data, bits);
    let mut p = Poly::zero();

    for i in 0..N {
        let packed = shifted.coeffs[i];
        // FIPS 204 Algorithm 18: centered = b - packed_value
        let centered = b as i32 - packed as i32;
        p.coeffs[i] = arithmetic::from_centered(centered);
    }

    p
}

//  Bit I/O Helpers 

/// Write `bits` least-significant bits of `val` into `buf` starting at
/// bit position `bit_pos`.  Little-endian byte order.
fn write_bits(buf: &mut [u8], bit_pos: usize, val: u64, bits: u32) {
    let mut remaining = bits;
    let mut v = val;
    let mut pos = bit_pos;

    while remaining > 0 {
        let byte_idx = pos / 8;
        let bit_offset = pos % 8;
        let can_write = core::cmp::min(remaining, (8 - bit_offset) as u32);

        let mask = (1u64 << can_write) - 1;
        buf[byte_idx] |= ((v & mask) as u8) << bit_offset;

        v >>= can_write;
        pos += can_write as usize;
        remaining -= can_write;
    }
}

/// Read `bits` bits from `data` starting at bit position `bit_pos`.
fn read_bits(data: &[u8], bit_pos: usize, bits: u32) -> u64 {
    let mut result: u64 = 0;
    let mut remaining = bits;
    let mut pos = bit_pos;
    let mut shift = 0u32;

    while remaining > 0 {
        let byte_idx = pos / 8;
        let bit_offset = pos % 8;
        let can_read = core::cmp::min(remaining, (8 - bit_offset) as u32);

        let mask = ((1u16 << can_read) - 1) as u8;
        let val = (data[byte_idx] >> bit_offset) & mask;

        result |= (val as u64) << shift;

        pos += can_read as usize;
        shift += can_read;
        remaining -= can_read;
    }

    result
}

/// Number of bits needed to represent values in [0, max_val].
///
/// bit_length(0) = 1 (we still need 1 bit to store the value 0).
/// bit_length(1) = 1
/// bit_length(2) = 2
/// bit_length(255) = 8
pub fn bit_length(max_val: u32) -> u32 {
    if max_val == 0 {
        return 1;
    }
    32 - max_val.leading_zeros()
}

//  Public Key Encoding 
// FIPS 204, Algorithm 22 (pkEncode) / Algorithm 23 (pkDecode)
//
// Public key = ρ (32 bytes) || t1[0] || t1[1] || ... || t1[k-1]
//
// t1 coefficients are in [0, 2^10 - 1] = [0, 1023], packed at 10 bits each.
// Each polynomial: 256 * 10 / 8 = 320 bytes.

/// Encode public key: pk = ρ || pack(t1[0]) || ... || pack(t1[k-1]).
pub fn pk_encode(rho: &[u8; 32], t1: &[Poly], params: &Params) -> Vec<u8> {
    debug_assert_eq!(t1.len(), params.k);

    let mut pk = Vec::with_capacity(params.pk_bytes);
    pk.extend_from_slice(rho);
    for poly in t1 {
        pk.extend_from_slice(&simple_bit_pack(poly, 10));
    }

    debug_assert_eq!(pk.len(), params.pk_bytes);
    pk
}

/// Decode public key: extract ρ and t1.
pub fn pk_decode(pk: &[u8], params: &Params) -> ([u8; 32], Vec<Poly>) {
    debug_assert_eq!(pk.len(), params.pk_bytes);

    let mut rho = [0u8; 32];
    rho.copy_from_slice(&pk[..32]);

    let bytes_per_poly = N * 10 / 8; // 320
    let mut t1 = Vec::with_capacity(params.k);
    let mut offset = 32;
    for _ in 0..params.k {
        t1.push(simple_bit_unpack(&pk[offset..offset + bytes_per_poly], 10));
        offset += bytes_per_poly;
    }

    (rho, t1)
}

//  Secret Key Encoding 
// FIPS 204, Algorithm 24 (skEncode) / Algorithm 25 (skDecode)
//
// Secret key layout:
//   ρ (32 bytes) || K (32 bytes) || tr (64 bytes)
//   || s1[0] || ... || s1[l-1]       (each η-bit packed)
//   || s2[0] || ... || s2[k-1]       (each η-bit packed)
//   || t0[0] || ... || t0[k-1]       (each D=13 bit packed)
//
// s1, s2 have centered coefficients in [-η, η], packed via BitPack(η, η).
// t0 has centered coefficients in [-(2^(d-1)-1), 2^(d-1)], packed via
//   BitPack(2^(d-1)-1, 2^(d-1)).

/// Encode secret key.
pub fn sk_encode(
    rho: &[u8; 32],
    k_seed: &[u8; 32],
    tr: &[u8; 64],
    s1: &[Poly],
    s2: &[Poly],
    t0: &[Poly],
    params: &Params,
) -> Vec<u8> {
    let mut sk = Vec::with_capacity(params.sk_bytes);

    // Fixed-length prefix
    sk.extend_from_slice(rho);
    sk.extend_from_slice(k_seed);
    sk.extend_from_slice(tr);

    // s1: l polynomials, coefficients in [-η, η]
    let eta = params.eta;
    for poly in s1 {
        sk.extend_from_slice(&bit_pack(poly, eta, eta));
    }

    // s2: k polynomials, coefficients in [-η, η]
    for poly in s2 {
        sk.extend_from_slice(&bit_pack(poly, eta, eta));
    }

    // t0: k polynomials, coefficients in [-(2^(d-1)-1), 2^(d-1)]
    let t0_neg = (1u32 << (D - 1)) - 1; // 2^12 - 1 = 4095
    let t0_pos = 1u32 << (D - 1);       // 2^12 = 4096
    for poly in t0 {
        sk.extend_from_slice(&bit_pack(poly, t0_neg, t0_pos));
    }

    debug_assert_eq!(sk.len(), params.sk_bytes,
        "SK encoding size mismatch: got {}, expected {}", sk.len(), params.sk_bytes);
    sk
}

/// Decode secret key: extract all components.
pub fn sk_decode(
    sk: &[u8],
    params: &Params,
) -> ([u8; 32], [u8; 32], [u8; 64], Vec<Poly>, Vec<Poly>, Vec<Poly>) {
    debug_assert_eq!(sk.len(), params.sk_bytes);

    let mut offset = 0;

    // ρ
    let mut rho = [0u8; 32];
    rho.copy_from_slice(&sk[offset..offset + 32]);
    offset += 32;

    // K
    let mut k_seed = [0u8; 32];
    k_seed.copy_from_slice(&sk[offset..offset + 32]);
    offset += 32;

    // tr
    let mut tr = [0u8; 64];
    tr.copy_from_slice(&sk[offset..offset + 64]);
    offset += 64;

    // s1
    let eta = params.eta;
    let s_packed_bytes = params.eta_packed_bytes();
    let mut s1 = Vec::with_capacity(params.l);
    for _ in 0..params.l {
        s1.push(bit_unpack(&sk[offset..offset + s_packed_bytes], eta, eta));
        offset += s_packed_bytes;
    }

    // s2
    let mut s2 = Vec::with_capacity(params.k);
    for _ in 0..params.k {
        s2.push(bit_unpack(&sk[offset..offset + s_packed_bytes], eta, eta));
        offset += s_packed_bytes;
    }

    // t0
    let t0_neg = (1u32 << (D - 1)) - 1;
    let t0_pos = 1u32 << (D - 1);
    let t0_packed_bytes = params.t0_packed_bytes();
    let mut t0 = Vec::with_capacity(params.k);
    for _ in 0..params.k {
        t0.push(bit_unpack(&sk[offset..offset + t0_packed_bytes], t0_neg, t0_pos));
        offset += t0_packed_bytes;
    }

    debug_assert_eq!(offset, params.sk_bytes);
    (rho, k_seed, tr, s1, s2, t0)
}

//  Signature Encoding 
// FIPS 204, Algorithm 26 (sigEncode) / Algorithm 27 (sigDecode)
//
// Signature layout:
//   c̃ (c_tilde_bytes) || z[0] || ... || z[l-1] || hint_encoding
//
// z coefficients are in [-(γ₁-1), γ₁-1], packed via BitPack(γ₁-1, γ₁).
// Hint uses a special sparse encoding.

/// Encode signature.
pub fn sig_encode(
    c_tilde: &[u8],
    z: &[Poly],
    h: &[Vec<bool>],
    params: &Params,
) -> Vec<u8> {
    let mut sig = Vec::with_capacity(params.sig_bytes);

    // c̃
    sig.extend_from_slice(c_tilde);

    // z: l polynomials in [-(γ₁-1), γ₁-1]
    let gamma1_minus1 = params.gamma1 - 1;
    for poly in z {
        sig.extend_from_slice(&bit_pack(poly, gamma1_minus1, gamma1_minus1));
    }

    // Hint
    sig.extend_from_slice(&hint_bit_pack(h, params));

    debug_assert_eq!(sig.len(), params.sig_bytes,
        "Sig encoding size mismatch: got {}, expected {}", sig.len(), params.sig_bytes);
    sig
}

/// Decode signature: extract c̃, z, and hint h.
/// Returns None if the hint encoding is malformed.
pub fn sig_decode(
    sig: &[u8],
    params: &Params,
) -> Option<(Vec<u8>, Vec<Poly>, Vec<Vec<bool>>)> {
    if sig.len() != params.sig_bytes {
        return None;
    }

    let mut offset = 0;

    // c̃
    let c_tilde = sig[offset..offset + params.c_tilde_bytes].to_vec();
    offset += params.c_tilde_bytes;

    // z
    let gamma1_minus1 = params.gamma1 - 1;
    let z_bits = bit_length(2 * gamma1_minus1);
    let z_packed_bytes = (N * z_bits as usize + 7) / 8;
    let mut z = Vec::with_capacity(params.l);
    for _ in 0..params.l {
        z.push(bit_unpack(&sig[offset..offset + z_packed_bytes], gamma1_minus1, gamma1_minus1));
        offset += z_packed_bytes;
    }

    // Hint
    let hint_bytes = params.omega + params.k;
    let h = hint_bit_unpack(&sig[offset..offset + hint_bytes], params)?;
    offset += hint_bytes;

    debug_assert_eq!(offset, params.sig_bytes);
    Some((c_tilde, z, h))
}

//  Hint Encoding 
// FIPS 204, Algorithm 20 (HintBitPack) / Algorithm 21 (HintBitUnpack)
//
// The hint vector h has at most ω ones across all k polynomials.
// The sparse encoding stores:
//   - The indices where h[i][j] = 1, packed consecutively (one byte each)
//   - Followed by k "end-of-list" markers at positions [ω, ω+1, ..., ω+k-1]
//     that record the cumulative count of ones up through each polynomial.
//
// Total size: ω + k bytes.

/// Encode the hint vector h into ω + k bytes.
///
/// h[i] is a Vec<bool> of length N for polynomial i.
pub fn hint_bit_pack(h: &[Vec<bool>], params: &Params) -> Vec<u8> {
    let total = params.omega + params.k;
    let mut buf = vec![0u8; total];

    let mut idx = 0; // Running index into the "data" portion
    for i in 0..params.k {
        for j in 0..N {
            if h[i][j] {
                buf[idx] = j as u8;
                idx += 1;
            }
        }
        // Store the cumulative count as the end-of-list marker
        buf[params.omega + i] = idx as u8;
    }

    // Remaining slots (idx..omega) stay 0
    buf
}

/// Decode the hint vector h from ω + k bytes.
///
/// Returns None if the encoding is malformed (indices out of order,
/// counts exceed ω, etc.).
pub fn hint_bit_unpack(data: &[u8], params: &Params) -> Option<Vec<Vec<bool>>> {
    let mut h = vec![vec![false; N]; params.k];

    let mut idx: usize = 0;
    for i in 0..params.k {
        let end = data[params.omega + i] as usize;

        // End marker must be non-decreasing and ≤ ω
        if end < idx || end > params.omega {
            return None;
        }

        // Read indices for polynomial i
        let mut prev: Option<u8> = None;
        while idx < end {
            let j = data[idx] as usize;

            // Index must be in [0, N)
            if j >= N {
                return None;
            }

            // Indices within each polynomial must be strictly increasing
            if let Some(p) = prev {
                if data[idx] <= p {
                    return None;
                }
            }
            prev = Some(data[idx]);

            h[i][j] = true;
            idx += 1;
        }
    }

    // Any remaining bytes between idx and omega must be zero
    for pos in idx..params.omega {
        if data[pos] != 0 {
            return None;
        }
    }

    Some(h)
}

//  Bit-level unpack (reverse of simple_bit_pack) 

/// Unpack 256 unsigned integers, each stored using `bits` bits, from `bytes`.
pub fn unpack_bits_to_vec(bytes: &[u8], bits: u32) -> Vec<u32> {
    let mut result = Vec::with_capacity(N);
    let mut bit_pos = 0usize;
    for _ in 0..N {
        let mut val = 0u32;
        for b in 0..bits as usize {
            let byte_idx = bit_pos / 8;
            let bit_idx  = bit_pos % 8;
            if byte_idx < bytes.len() {
                val |= (((bytes[byte_idx] >> bit_idx) & 1) as u32) << b;
            }
            bit_pos += 1;
        }
        result.push(val);
    }
    result
}

//  t1 decode (10 bits per coefficient) 

pub fn unpack_t1(bytes: &[u8]) -> crate::poly::Poly {
    let vals = unpack_bits_to_vec(bytes, 10);
    let mut p = crate::poly::Poly::default();
    for i in 0..N { p.coeffs[i] = vals[i]; }
    p
}

//  Hint encode / decode (FIPS 204) 

/// Encode the hint matrix h into (omega + k) bytes.
pub fn encode_hint(h: &[Vec<bool>], params: &Params) -> Vec<u8> {
    let mut out = vec![0u8; params.omega + params.k];
    let mut index = 0usize;
    for (i, row) in h.iter().enumerate() {
        for (j, &bit) in row.iter().enumerate() {
            if bit {
                out[index] = j as u8;
                index += 1;
            }
        }
        out[params.omega + i] = index as u8;
    }
    out
}

/// Decode hint bytes back to a k×256 boolean matrix. Returns None if malformed.
pub fn decode_hint(bytes: &[u8], params: &Params) -> Option<Vec<Vec<bool>>> {
    if bytes.len() < params.omega + params.k { return None; }
    let mut h = vec![vec![false; N]; params.k];
    let mut index = 0usize;
    for i in 0..params.k {
        let end = bytes[params.omega + i] as usize;
        if end < index || end > params.omega { return None; }
        for pos in index..end {
            let j = bytes[pos] as usize;
            if j >= N { return None; }
            // Within the same row hints must be strictly increasing
            if pos > index && bytes[pos] <= bytes[pos - 1] { return None; }
            h[i][j] = true;
        }
        index = end;
    }
    Some(h)
}

//  z coefficient encode / decode 

/// Encode z polynomial: each centered coefficient c → c + gamma1, packed in gamma1_bits bits.
pub fn encode_z(poly: &crate::poly::Poly, gamma1: u32, bits: u32) -> Vec<u8> {
    let mut temp = crate::poly::Poly::default();
    for i in 0..N {
        let c = crate::arithmetic::to_centered(poly.coeffs[i]);
        temp.coeffs[i] = (c + gamma1 as i32) as u32;
    }
    simple_bit_pack(&temp, bits)
}

/// Decode z polynomial from packed bytes.
pub fn decode_z(bytes: &[u8], gamma1: u32, bits: u32) -> crate::poly::Poly {
    let vals = unpack_bits_to_vec(bytes, bits);
    let mut p = crate::poly::Poly::default();
    for i in 0..N {
        let c = vals[i] as i32 - gamma1 as i32;
        p.coeffs[i] = crate::arithmetic::from_centered(c);
    }
    p
}

//  Tests 

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ML_DSA_44, ML_DSA_65, ML_DSA_87};

    // ── Bit helpers ──

    #[test]
    fn test_bit_length() {
        assert_eq!(bit_length(0), 1);
        assert_eq!(bit_length(1), 1);
        assert_eq!(bit_length(2), 2);
        assert_eq!(bit_length(3), 2);
        assert_eq!(bit_length(4), 3);
        assert_eq!(bit_length(255), 8);
        assert_eq!(bit_length(256), 9);
        assert_eq!(bit_length(1023), 10);
        assert_eq!(bit_length(1024), 11);
    }

    #[test]
    fn test_write_read_bits_roundtrip() {
        // Write various values at various positions and read them back.
        let mut buf = vec![0u8; 32];

        write_bits(&mut buf, 0, 42, 8);
        assert_eq!(read_bits(&buf, 0, 8), 42);

        write_bits(&mut buf, 8, 1023, 10);
        assert_eq!(read_bits(&buf, 8, 10), 1023);

        // Write at a non-byte-aligned position
        let mut buf2 = vec![0u8; 32];
        write_bits(&mut buf2, 3, 0b11011, 5);
        assert_eq!(read_bits(&buf2, 3, 5), 0b11011);
    }

    #[test]
    fn test_write_read_bits_adjacent() {
        // Write two values back-to-back at non-aligned positions
        let mut buf = vec![0u8; 32];
        write_bits(&mut buf, 0, 7, 3);       // 3 bits
        write_bits(&mut buf, 3, 15, 4);      // 4 bits
        write_bits(&mut buf, 7, 255, 8);     // 8 bits

        assert_eq!(read_bits(&buf, 0, 3), 7);
        assert_eq!(read_bits(&buf, 3, 4), 15);
        assert_eq!(read_bits(&buf, 7, 8), 255);
    }

    // ── SimpleBitPack / SimpleBitUnpack ──

    #[test]
    fn test_simple_pack_unpack_10_bits() {
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i as u32 * 3 + 7) % 1024; // [0, 1023]
        }

        let packed = simple_bit_pack(&p, 10);
        assert_eq!(packed.len(), 256 * 10 / 8, "10-bit pack should be 320 bytes");

        let unpacked = simple_bit_unpack(&packed, 10);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i],
                "10-bit roundtrip failed at {}: expected {}, got {}",
                i, p.coeffs[i], unpacked.coeffs[i]);
        }
    }

    #[test]
    fn test_simple_pack_unpack_13_bits() {
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i as u32 * 17 + 3) % 8192; // [0, 8191]
        }

        let packed = simple_bit_pack(&p, 13);
        assert_eq!(packed.len(), 256 * 13 / 8); // 416 bytes

        let unpacked = simple_bit_unpack(&packed, 13);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i],
                "13-bit roundtrip failed at {}", i);
        }
    }

    #[test]
    fn test_simple_pack_unpack_3_bits() {
        // η=2: values in [0, 4], 3 bits
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i as u32) % 5; // [0, 4]
        }

        let packed = simple_bit_pack(&p, 3);
        assert_eq!(packed.len(), 256 * 3 / 8); // 96 bytes

        let unpacked = simple_bit_unpack(&packed, 3);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i],
                "3-bit roundtrip failed at {}", i);
        }
    }

    #[test]
    fn test_simple_pack_unpack_4_bits() {
        // η=4: values in [0, 8], 4 bits
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i as u32) % 9;
        }

        let packed = simple_bit_pack(&p, 4);
        assert_eq!(packed.len(), 256 * 4 / 8); // 128 bytes

        let unpacked = simple_bit_unpack(&packed, 4);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i]);
        }
    }

    #[test]
    fn test_simple_pack_unpack_18_bits() {
        // γ₁ = 2^17: z values use 18 bits
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i as u32 * 1031 + 5) % (1 << 18);
        }

        let packed = simple_bit_pack(&p, 18);
        assert_eq!(packed.len(), 256 * 18 / 8); // 576

        let unpacked = simple_bit_unpack(&packed, 18);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i]);
        }
    }

    #[test]
    fn test_simple_pack_unpack_20_bits() {
        // γ₁ = 2^19: z values use 20 bits
        let mut p = Poly::zero();
        for i in 0..N {
            p.coeffs[i] = (i as u32 * 4099 + 13) % (1 << 20);
        }

        let packed = simple_bit_pack(&p, 20);
        assert_eq!(packed.len(), 256 * 20 / 8); // 640

        let unpacked = simple_bit_unpack(&packed, 20);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i]);
        }
    }

    // ── BitPack / BitUnpack (signed) ──

    #[test]
    fn test_bit_pack_unpack_eta2() {
        // η=2: centered coefficients in [-2, 2]
        let mut p = Poly::zero();
        for i in 0..N {
            let centered = (i % 5) as i32 - 2; // cycles through -2,-1,0,1,2
            p.coeffs[i] = arithmetic::from_centered(centered);
        }

        let packed = bit_pack(&p, 2, 2);
        assert_eq!(packed.len(), ML_DSA_44.eta_packed_bytes()); // 96

        let unpacked = bit_unpack(&packed, 2, 2);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i],
                "η=2 bit pack roundtrip failed at {}", i);
        }
    }

    #[test]
    fn test_bit_pack_unpack_eta4() {
        // η=4: centered coefficients in [-4, 4]
        let mut p = Poly::zero();
        for i in 0..N {
            let centered = (i % 9) as i32 - 4;
            p.coeffs[i] = arithmetic::from_centered(centered);
        }

        let packed = bit_pack(&p, 4, 4);
        assert_eq!(packed.len(), ML_DSA_65.eta_packed_bytes()); // 128

        let unpacked = bit_unpack(&packed, 4, 4);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i],
                "η=4 bit pack roundtrip failed at {}", i);
        }
    }

    #[test]
    fn test_bit_pack_unpack_t0() {
        // t0: coefficients in [-(2^12 - 1), 2^12] = [-4095, 4096]
        let a = (1u32 << (D - 1)) - 1; // 4095
        let b = 1u32 << (D - 1);       // 4096

        let mut p = Poly::zero();
        for i in 0..N {
            // Cycle through the range
            let total_range = (a + b + 1) as usize; // 8192
            let centered = (i % total_range) as i32 - a as i32;
            p.coeffs[i] = arithmetic::from_centered(centered);
        }

        let packed = bit_pack(&p, a, b);
        assert_eq!(packed.len(), ML_DSA_44.t0_packed_bytes()); // 416

        let unpacked = bit_unpack(&packed, a, b);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i],
                "t0 bit pack roundtrip failed at {}: centered {}",
                i, arithmetic::to_centered(p.coeffs[i]));
        }
    }

    #[test]
    fn test_bit_pack_unpack_z_gamma1_17() {
        // z: coefficients in [-(γ₁-1), γ₁-1] for γ₁ = 2^17
        let gamma1 = 1u32 << 17;
        let a = gamma1 - 1; // 131071

        let mut p = Poly::zero();
        for i in 0..N {
            let centered = ((i as i32 * 1031) % (2 * a as i32 + 1)) - a as i32;
            p.coeffs[i] = arithmetic::from_centered(centered);
        }

        let packed = bit_pack(&p, a, a);
        let unpacked = bit_unpack(&packed, a, a);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i],
                "z (γ₁=2^17) roundtrip failed at {}", i);
        }
    }

    // ── Public Key encoding ──

    #[test]
    fn test_pk_encode_decode_roundtrip() {
        let rho = [0x42u8; 32];

        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let mut t1 = Vec::new();
            for k in 0..params.k {
                let mut p = Poly::zero();
                for i in 0..N {
                    p.coeffs[i] = ((k * N + i) as u32 * 5 + 3) % 1024;
                }
                t1.push(p);
            }

            let pk = pk_encode(&rho, &t1, params);
            assert_eq!(pk.len(), params.pk_bytes, "{} pk size", params.name);

            let (rho_dec, t1_dec) = pk_decode(&pk, params);
            assert_eq!(rho, rho_dec);
            for k in 0..params.k {
                for i in 0..N {
                    assert_eq!(t1[k].coeffs[i], t1_dec[k].coeffs[i],
                        "{}: t1[{}][{}] mismatch", params.name, k, i);
                }
            }
        }
    }

    // ── Secret Key encoding ──

    #[test]
    fn test_sk_encode_decode_roundtrip() {
        let rho = [0xAA; 32];
        let k_seed = [0xBB; 32];
        let tr = [0xCC; 64];

        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let eta = params.eta;

            // Create s1 with coefficients in [-η, η]
            let mut s1 = Vec::new();
            for j in 0..params.l {
                let mut p = Poly::zero();
                for i in 0..N {
                    let c = ((j * N + i) % (2 * eta as usize + 1)) as i32 - eta as i32;
                    p.coeffs[i] = arithmetic::from_centered(c);
                }
                s1.push(p);
            }

            // Create s2
            let mut s2 = Vec::new();
            for j in 0..params.k {
                let mut p = Poly::zero();
                for i in 0..N {
                    let c = ((j * N + i + 7) % (2 * eta as usize + 1)) as i32 - eta as i32;
                    p.coeffs[i] = arithmetic::from_centered(c);
                }
                s2.push(p);
            }

            // Create t0 with coefficients in [-(2^12-1), 2^12]
            let t0_bound = (1u32 << (D - 1)) - 1;
            let mut t0 = Vec::new();
            for j in 0..params.k {
                let mut p = Poly::zero();
                for i in 0..N {
                    let range = (2 * t0_bound + 2) as usize; // 8192
                    let c = ((j * N + i + 13) % range) as i32 - t0_bound as i32;
                    p.coeffs[i] = arithmetic::from_centered(c);
                }
                t0.push(p);
            }

            let sk = sk_encode(&rho, &k_seed, &tr, &s1, &s2, &t0, params);
            assert_eq!(sk.len(), params.sk_bytes, "{} sk size", params.name);

            let (rho_d, k_d, tr_d, s1_d, s2_d, t0_d) = sk_decode(&sk, params);

            assert_eq!(rho, rho_d, "{} rho", params.name);
            assert_eq!(k_seed, k_d, "{} K", params.name);
            assert_eq!(tr, tr_d, "{} tr", params.name);

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

    // ── Hint encoding ──

    #[test]
    fn test_hint_pack_unpack_empty() {
        // No hints set (all false)
        let params = &ML_DSA_44;
        let h = vec![vec![false; N]; params.k];

        let packed = hint_bit_pack(&h, params);
        assert_eq!(packed.len(), params.omega + params.k);

        let unpacked = hint_bit_unpack(&packed, params).expect("Valid hint");
        for i in 0..params.k {
            for j in 0..N {
                assert_eq!(h[i][j], unpacked[i][j]);
            }
        }
    }

    #[test]
    fn test_hint_pack_unpack_some_hints() {
        let params = &ML_DSA_44;
        let mut h = vec![vec![false; N]; params.k];

        // Set a few hints in different polynomials
        h[0][0] = true;
        h[0][100] = true;
        h[0][255] = true;
        h[1][50] = true;
        h[3][200] = true;

        let packed = hint_bit_pack(&h, params);
        let unpacked = hint_bit_unpack(&packed, params).expect("Valid hint");

        for i in 0..params.k {
            for j in 0..N {
                assert_eq!(h[i][j], unpacked[i][j],
                    "Hint mismatch at h[{}][{}]", i, j);
            }
        }
    }

    #[test]
    fn test_hint_pack_unpack_max_hints() {
        // Set exactly ω hints
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            let mut h = vec![vec![false; N]; params.k];
            let mut count = 0;
            'outer: for i in 0..params.k {
                for j in 0..N {
                    if count >= params.omega {
                        break 'outer;
                    }
                    h[i][j] = true;
                    count += 1;
                }
            }

            let packed = hint_bit_pack(&h, params);
            let unpacked = hint_bit_unpack(&packed, params).expect("Valid hint");

            for i in 0..params.k {
                for j in 0..N {
                    assert_eq!(h[i][j], unpacked[i][j],
                        "{}: hint max test h[{}][{}]", params.name, i, j);
                }
            }
        }
    }

    #[test]
    fn test_hint_unpack_rejects_bad_order() {
        // Indices within a polynomial must be strictly increasing.
        let params = &ML_DSA_44;
        let mut data = vec![0u8; params.omega + params.k];

        // Put indices 5, 3 (wrong order) in polynomial 0
        data[0] = 5;
        data[1] = 3; // Not strictly increasing!
        data[params.omega] = 2; // End marker: 2 entries in poly 0

        let result = hint_bit_unpack(&data, params);
        assert!(result.is_none(), "Should reject non-increasing indices");
    }

    // ── Packed size validation ──

    #[test]
    fn test_packed_sizes_match_params() {
        // Verify that our packing produces exactly the right number of bytes
        // for every component, matching the sizes in params.
        for params in &[&ML_DSA_44, &ML_DSA_65, &ML_DSA_87] {
            // t1: 10 bits
            let t1_bytes = simple_bit_pack(&Poly::zero(), 10).len();
            assert_eq!(t1_bytes, params.t1_packed_bytes(), "{} t1", params.name);

            // η-packing
            let eta_bytes = bit_pack(&Poly::zero(), params.eta, params.eta).len();
            assert_eq!(eta_bytes, params.eta_packed_bytes(), "{} eta", params.name);

            // t0: D=13 bits (asymmetric: a=4095, b=4096)
            let t0_a = (1u32 << (D - 1)) - 1;
            let t0_b = 1u32 << (D - 1);
            let t0_bytes = bit_pack(&Poly::zero(), t0_a, t0_b).len();
            assert_eq!(t0_bytes, params.t0_packed_bytes(), "{} t0", params.name);

            // z: γ₁ bits
            let g1m1 = params.gamma1 - 1;
            let z_bytes = bit_pack(&Poly::zero(), g1m1, g1m1).len();
            let expected_z_bytes = N * params.gamma1_bits() as usize / 8;
            assert_eq!(z_bytes, expected_z_bytes, "{} z", params.name);
        }
    }
}
