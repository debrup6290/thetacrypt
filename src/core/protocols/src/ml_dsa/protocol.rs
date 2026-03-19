use std::collections::HashMap;
use std::sync::Arc;

use log::{info, warn};
use theta_network::types::message::NetMessage;
use theta_schemes::keys::keys::PrivateKeyShare;
use theta_schemes::pq_schemes::ml_dsa::MlDsaPartyKey;
use tonic::async_trait;

use ml_dsa::threshold::protocol::{
    share_sign_1, share_sign_2, share_sign_3,
    Round1Msg, Round2Msg, Round3Msg, SigningState,
};
use ml_dsa::threshold::combine;
use ml_dsa::threshold::params::ThresholdParams;
use ml_dsa::threshold::rss::VerificationKey;

use crate::interface::{ProtocolError, ThresholdRoundProtocol};
use super::message_types::{MlDsaMessage, MlDsaR1Payload, MlDsaR2Payload, MlDsaR3Payload};

pub struct MlDsaProtocol {
    round:        u8,
    party_key:    Arc<MlDsaPartyKey>,
    tp:           ThresholdParams,
    vk:           VerificationKey,
    msg:          Vec<u8>,
    act:          Vec<usize>,
    signing_seed: Vec<u8>,

    own_state:    Option<SigningState>,
    own_r1:       Option<Round1Msg>,

    all_r1: HashMap<usize, Round1Msg>,
    all_r2: HashMap<usize, Round2Msg>,
    all_r3: HashMap<usize, Round3Msg>,

    result:   Option<Vec<u8>>,
    finished: bool,
}

impl MlDsaProtocol {
    pub fn new(
        private_key: Arc<PrivateKeyShare>,
        msg: Vec<u8>,
        act: Vec<usize>,
        signing_seed: Vec<u8>,
    ) -> Result<Self, ProtocolError> {
        let party_key = match private_key.as_ref() {
            PrivateKeyShare::MlDsa44(k) => Arc::new(k.clone()),
            PrivateKeyShare::MlDsa65(k) => Arc::new(k.clone()),
            PrivateKeyShare::MlDsa87(k) => Arc::new(k.clone()),
            _ => return Err(ProtocolError::InternalError),
        };

        let tp = party_key.to_threshold_params()
            .map_err(|_| ProtocolError::InternalError)?;
        let vk = party_key.to_verification_key()
            .map_err(|_| ProtocolError::InternalError)?;

        Ok(Self {
            round: 1,
            party_key,
            tp,
            vk,
            msg,
            act,
            signing_seed,
            own_state: None,
            own_r1: None,
            all_r1: HashMap::new(),
            all_r2: HashMap::new(),
            all_r3: HashMap::new(),
            result: None,
            finished: false,
        })
    }
}

impl ThresholdRoundProtocol<NetMessage> for MlDsaProtocol {
    type ProtocolMessage = MlDsaMessage;

    fn do_round(&mut self) -> Result<Self::ProtocolMessage, ProtocolError> {
        match self.round {
            1 => {
                let pk = self.party_key.to_party_key()
                    .map_err(|_| ProtocolError::InternalError)?;
                let (state, r1_msg) = share_sign_1(
                    &self.vk, &pk, &self.tp, &self.signing_seed,
                );
                let payload = MlDsaR1Payload {
                    party_id: pk.party_id,
                    commitments: r1_msg.commitments.clone(),
                };
                self.all_r1.insert(pk.party_id, r1_msg.clone_msg());
                self.own_r1 = Some(r1_msg);
                self.own_state = Some(state);
                self.round = 2;
                Ok(MlDsaMessage::Round1(payload))
            }
            2 => {
                let state = self.own_state.as_mut()
                    .ok_or(ProtocolError::InternalError)?;
                let r1_msgs: Vec<Round1Msg> = self.all_r1.values()
                    .map(|m| m.clone_msg())
                    .collect();
                let r2_msg = share_sign_2(
                    state, &self.vk, &self.act, &self.msg, r1_msgs,
                ).ok_or(ProtocolError::InternalError)?;
                let payload = MlDsaR2Payload {
                    party_id: r2_msg.party_id,
                    w_values_bytes: bincode::serialize(&r2_msg.w_values)
                        .map_err(|_| ProtocolError::InternalError)?,
                };
                self.all_r2.insert(r2_msg.party_id, r2_msg.clone_msg());
                self.round = 3;
                Ok(MlDsaMessage::Round2(payload))
            }
            3 => {
                let state = self.own_state.as_ref()
                    .ok_or(ProtocolError::InternalError)?;
                let pk = self.party_key.to_party_key()
                    .map_err(|_| ProtocolError::InternalError)?;
                let r2_msgs: Vec<Round2Msg> = self.all_r2.values()
                    .map(|m| m.clone_msg())
                    .collect();
                let r3_msg = share_sign_3(state, &self.vk, &pk, &r2_msgs)
                    .ok_or(ProtocolError::InternalError)?;
                let payload = MlDsaR3Payload {
                    party_id: r3_msg.party_id,
                    z1_values_bytes: bincode::serialize(&r3_msg.z1_values)
                        .map_err(|_| ProtocolError::InternalError)?,
                };
                self.all_r3.insert(r3_msg.party_id, r3_msg.clone_msg());
                self.round = 4;
                Ok(MlDsaMessage::Round3(payload))
            }
            _ => Err(ProtocolError::InvalidRound),
        }
    }

    fn is_ready_for_next_round(&self) -> bool {
        match self.round {
            2 => self.all_r1.len() >= self.tp.threshold,
            3 => self.all_r2.len() >= self.tp.threshold,
            _ => false,
        }
    }

    fn is_ready_to_finalize(&self) -> bool {
        self.round == 4 && self.all_r3.len() >= self.tp.threshold
    }

    fn finalize(&mut self) -> Result<Vec<u8>, ProtocolError> {
        let r2_msgs: Vec<Round2Msg> = self.all_r2.values()
            .map(|m| m.clone_msg())
            .collect();
        let r3_msgs: Vec<Round3Msg> = self.all_r3.values()
            .map(|m| m.clone_msg())
            .collect();

        let sig = combine::combine(
            &self.vk, &self.act, &self.msg,
            &r2_msgs, &r3_msgs, &self.tp,
        ).ok_or(ProtocolError::SchemeError(
            theta_schemes::interface::SchemeError::Aborted(
                "All K instances rejected — retry with fresh seed".into()
            )
        ))?;

        let sig_bytes = sig.to_bytes(&self.tp.base);
        self.result = Some(sig_bytes.clone());
        self.finished = true;
        info!("MlDsaProtocol: signature assembled successfully.");
        Ok(sig_bytes)
    }

    fn update(&mut self, message: Self::ProtocolMessage) -> Result<(), ProtocolError> {
        match message {
            MlDsaMessage::Round1(payload) => {
                if !self.all_r1.contains_key(&payload.party_id) {
                    self.all_r1.insert(payload.party_id, Round1Msg {
                        party_id: payload.party_id,
                        commitments: payload.commitments,
                    });
                }
            }
            MlDsaMessage::Round2(payload) => {
                if !self.all_r2.contains_key(&payload.party_id) {
                    let w_values = bincode::deserialize(&payload.w_values_bytes)
                        .map_err(|_| ProtocolError::InvalidShare)?;
                    self.all_r2.insert(payload.party_id, Round2Msg {
                        party_id: payload.party_id,
                        w_values,
                    });
                }
            }
            MlDsaMessage::Round3(payload) => {
                if !self.all_r3.contains_key(&payload.party_id) {
                    let z1_values = bincode::deserialize(&payload.z1_values_bytes)
                        .map_err(|_| ProtocolError::InvalidShare)?;
                    self.all_r3.insert(payload.party_id, Round3Msg {
                        party_id: payload.party_id,
                        z1_values,
                    });
                }
            }
            MlDsaMessage::Default => {}
        }
        Ok(())
    }
}
