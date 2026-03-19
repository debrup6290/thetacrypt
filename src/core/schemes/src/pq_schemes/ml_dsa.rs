//! Adapter layer: bridges threasold-ml-dsa types into Thetacrypt's type system.

use ml_dsa::params::{Params, ML_DSA_44, ML_DSA_65, ML_DSA_87};
use ml_dsa::threshold::params as thresh_params;
use ml_dsa::threshold::rss::{PartyKey, VerificationKey};
use theta_proto::scheme_types::{Group, ThresholdScheme};
use crate::interface::{SchemeError, Serializable};
use crate::keys::keys::calc_key_id;

pub fn scheme_to_params(scheme: ThresholdScheme) -> &'static Params {
    match scheme {
        ThresholdScheme::MlDsa44 => &ML_DSA_44,
        ThresholdScheme::MlDsa65 => &ML_DSA_65,
        ThresholdScheme::MlDsa87 => &ML_DSA_87,
        _ => panic!("scheme_to_params: not an ML-DSA scheme"),
    }
}

fn scheme_tag(scheme: ThresholdScheme) -> u8 {
    match scheme {
        ThresholdScheme::MlDsa44 => 44,
        ThresholdScheme::MlDsa65 => 65,
        ThresholdScheme::MlDsa87 => 87,
        _ => panic!("scheme_tag: not an ML-DSA scheme"),
    }
}

fn tag_to_scheme(tag: u8) -> Result<ThresholdScheme, SchemeError> {
    match tag {
        44 => Ok(ThresholdScheme::MlDsa44),
        65 => Ok(ThresholdScheme::MlDsa65),
        87 => Ok(ThresholdScheme::MlDsa87),
        _  => Err(SchemeError::UnknownScheme),
    }
}

#[derive(Clone, Debug)]
pub struct MlDsaPublicKey {
    id:     String,
    n:      u16,
    k:      u16,
    scheme: ThresholdScheme,
    group:  Group,
    packed: Vec<u8>,
    tr:     Vec<u8>,
}

impl MlDsaPublicKey {
    pub fn new(vk: &VerificationKey, n: usize, k: usize, scheme: ThresholdScheme) -> Self {
        Self {
            id:     calc_key_id(&vk.packed),
            n:      n as u16,
            k:      k as u16,
            scheme,
            group:  Group::Lattice,
            packed: vk.packed.clone(),
            tr:     vk.tr.to_vec(),
        }
    }
    pub fn get_key_id(&self)    -> &str            { &self.id }
    pub fn get_group(&self)     -> &Group          { &self.group }
    pub fn get_scheme(&self)    -> ThresholdScheme { self.scheme }
    pub fn get_threshold(&self) -> u16             { self.k }
    pub fn get_n(&self)         -> u16             { self.n }
    pub fn get_packed(&self)    -> &[u8]           { &self.packed }
    pub fn get_tr(&self)        -> &[u8]           { &self.tr }
    pub fn get_params(&self)    -> &'static Params { scheme_to_params(self.scheme) }

    pub fn to_ml_dsa_pk(&self) -> Option<ml_dsa::signature::PublicKey> {
        ml_dsa::signature::PublicKey::from_packed(self.packed.clone(), self.get_params())
    }
}

impl PartialEq for MlDsaPublicKey {
    fn eq(&self, other: &Self) -> bool { self.id == other.id }
}

impl Serializable for MlDsaPublicKey {
    fn to_bytes(&self) -> Result<Vec<u8>, SchemeError> {
        let mut out = Vec::new();
        out.push(scheme_tag(self.scheme));
        out.extend_from_slice(&self.n.to_le_bytes());
        out.extend_from_slice(&self.k.to_le_bytes());
        out.extend_from_slice(&(self.packed.len() as u32).to_le_bytes());
        out.extend_from_slice(&self.packed);
        out.extend_from_slice(&self.tr);
        Ok(out)
    }
    fn from_bytes(bytes: &Vec<u8>) -> Result<Self, SchemeError> {
        let min = 1 + 2 + 2 + 4 + 64;
        if bytes.len() < min { return Err(SchemeError::DeserializationFailed); }
        let mut p = 0;
        let scheme = tag_to_scheme(bytes[p])?; p += 1;
        let n      = u16::from_le_bytes(bytes[p..p+2].try_into().unwrap()); p += 2;
        let k      = u16::from_le_bytes(bytes[p..p+2].try_into().unwrap()); p += 2;
        let plen   = u32::from_le_bytes(bytes[p..p+4].try_into().unwrap()) as usize; p += 4;
        if bytes.len() < p + plen + 64 { return Err(SchemeError::DeserializationFailed); }
        let packed = bytes[p..p+plen].to_vec(); p += plen;
        let tr     = bytes[p..p+64].to_vec();
        Ok(Self { id: calc_key_id(&packed), n, k, scheme, group: Group::Lattice, packed, tr })
    }
}

#[derive(Clone, Debug)]
pub struct MlDsaPartyKey {
    id:           String,
    party_id:     usize,
    n:            u16,
    k:            u16,
    scheme:       ThresholdScheme,
    group:        Group,
    vk_bytes:     Vec<u8>,
    shares_bytes: Vec<u8>,
}

impl MlDsaPartyKey {
    pub fn new(
        party_key: &PartyKey,
        vk: &VerificationKey,
        n: usize,
        k: usize,
        scheme: ThresholdScheme,
    ) -> Result<Self, SchemeError> {
        let vk_bytes = bincode::serialize(vk)
            .map_err(|_| SchemeError::SerializationFailed)?;
        let shares_bytes = bincode::serialize(&party_key.shares)
            .map_err(|_| SchemeError::SerializationFailed)?;
        Ok(Self {
            id: calc_key_id(&vk.packed),
            party_id: party_key.party_id,
            n: n as u16,
            k: k as u16,
            scheme,
            group: Group::Lattice,
            vk_bytes,
            shares_bytes,
        })
    }
    pub fn get_key_id(&self)    -> &str            { &self.id }
    pub fn get_share_id(&self)  -> u16             { self.party_id as u16 }
    pub fn get_group(&self)     -> &Group          { &self.group }
    pub fn get_scheme(&self)    -> ThresholdScheme { self.scheme }
    pub fn get_threshold(&self) -> u16             { self.k }
    pub fn get_n(&self)         -> u16             { self.n }
    pub fn get_params(&self)    -> &'static Params { scheme_to_params(self.scheme) }

    pub fn get_public_key(&self) -> MlDsaPublicKey {
        let vk: VerificationKey = bincode::deserialize(&self.vk_bytes).unwrap();
        MlDsaPublicKey::new(&vk, self.n as usize, self.k as usize, self.scheme)
    }
    pub fn to_party_key(&self) -> Result<PartyKey, SchemeError> {
        let shares = bincode::deserialize(&self.shares_bytes)
            .map_err(|_| SchemeError::DeserializationFailed)?;
        let vk: VerificationKey = bincode::deserialize(&self.vk_bytes)
            .map_err(|_| SchemeError::DeserializationFailed)?;
        Ok(PartyKey { party_id: self.party_id, tr: vk.tr, shares })
    }
    pub fn to_verification_key(&self) -> Result<VerificationKey, SchemeError> {
        bincode::deserialize(&self.vk_bytes)
            .map_err(|_| SchemeError::DeserializationFailed)
    }
    pub fn to_threshold_params(&self) -> Result<thresh_params::ThresholdParams, SchemeError> {
        thresh_params::lookup(self.n as usize, self.k as usize, self.get_params())
            .ok_or(SchemeError::InvalidParams(None))
    }
}

impl PartialEq for MlDsaPartyKey {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id && self.party_id == other.party_id
    }
}

impl Serializable for MlDsaPartyKey {
    fn to_bytes(&self) -> Result<Vec<u8>, SchemeError> {
        let mut out = Vec::new();
        out.push(scheme_tag(self.scheme));
        out.extend_from_slice(&self.n.to_le_bytes());
        out.extend_from_slice(&self.k.to_le_bytes());
        out.extend_from_slice(&(self.party_id as u64).to_le_bytes());
        out.extend_from_slice(&(self.vk_bytes.len() as u32).to_le_bytes());
        out.extend_from_slice(&self.vk_bytes);
        out.extend_from_slice(&(self.shares_bytes.len() as u32).to_le_bytes());
        out.extend_from_slice(&self.shares_bytes);
        Ok(out)
    }
    fn from_bytes(bytes: &Vec<u8>) -> Result<Self, SchemeError> {
        let min = 1 + 2 + 2 + 8 + 4 + 4;
        if bytes.len() < min { return Err(SchemeError::DeserializationFailed); }
        let mut p = 0;
        let scheme   = tag_to_scheme(bytes[p])?; p += 1;
        let n        = u16::from_le_bytes(bytes[p..p+2].try_into().unwrap()); p += 2;
        let k        = u16::from_le_bytes(bytes[p..p+2].try_into().unwrap()); p += 2;
        let party_id = u64::from_le_bytes(bytes[p..p+8].try_into().unwrap()) as usize; p += 8;
        let vk_len   = u32::from_le_bytes(bytes[p..p+4].try_into().unwrap()) as usize; p += 4;
        if bytes.len() < p + vk_len + 4 { return Err(SchemeError::DeserializationFailed); }
        let vk_bytes     = bytes[p..p+vk_len].to_vec(); p += vk_len;
        let sh_len       = u32::from_le_bytes(bytes[p..p+4].try_into().unwrap()) as usize; p += 4;
        if bytes.len() < p + sh_len { return Err(SchemeError::DeserializationFailed); }
        let shares_bytes = bytes[p..p+sh_len].to_vec();
        let vk: VerificationKey = bincode::deserialize(&vk_bytes)
            .map_err(|_| SchemeError::DeserializationFailed)?;
        Ok(Self {
            id: calc_key_id(&vk.packed),
            party_id, n, k, scheme,
            group: Group::Lattice,
            vk_bytes, shares_bytes,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MlDsaSignatureShare {
    pub party_id: usize,
    pub scheme:   ThresholdScheme,
}
impl MlDsaSignatureShare {
    pub fn get_id(&self)     -> u16            { self.party_id as u16 }
    pub fn get_group(&self)  -> &Group         { &Group::Lattice }
    pub fn get_scheme(&self) -> ThresholdScheme { self.scheme }
    pub fn get_label(&self)  -> &[u8]          { &[] }
}
impl Serializable for MlDsaSignatureShare {
    fn to_bytes(&self) -> Result<Vec<u8>, SchemeError> {
        let mut out = (self.party_id as u64).to_le_bytes().to_vec();
        out.push(self.scheme as u8);
        Ok(out)
    }
    fn from_bytes(bytes: &Vec<u8>) -> Result<Self, SchemeError> {
        if bytes.len() < 9 { return Err(SchemeError::DeserializationFailed); }
        let party_id = u64::from_le_bytes(bytes[0..8].try_into().unwrap()) as usize;
        let scheme = ThresholdScheme::from_i32(bytes[8] as i32)
            .ok_or(SchemeError::DeserializationFailed)?;
        Ok(Self { party_id, scheme })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MlDsaSignature {
    pub bytes: Vec<u8>,
}

impl Serializable for MlDsaSignature {
    fn to_bytes(&self) -> Result<Vec<u8>, SchemeError> { Ok(self.bytes.clone()) }
    fn from_bytes(bytes: &Vec<u8>) -> Result<Self, SchemeError> {
        Ok(Self { bytes: bytes.clone() })
    }
}
