use std::collections::HashMap;
use std::sync::Arc;

use log::info;
use theta_network::types::message::NetMessage;
use theta_schemes::keys::keys::PrivateKeyShare;
use theta_schemes::pq_schemes::ml_dsa::{MlDsaPartyKey, MlDsaSignature};
use theta_schemes::interface::{Serializable, Signature};


use ml_dsa::threshold::protocol::{
    share_sign_1, share_sign_2, share_sign_3,
    Round1Msg, Round2Msg, Round3Msg, SigningState,
};
use ml_dsa::threshold::combine;
use theta_proto::scheme_types::ThresholdScheme;
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
    signing_seed: Vec<u8>,

    // ── T-of-N fields ──
    n_parties:    usize,
    threshold:    usize,
    my_party_id:  usize,
    /// Active set — computed deterministically after round 2.
    act:          Vec<usize>,

    own_state:    Option<SigningState>,
    own_r1:       Option<Round1Msg>,

    all_r1: HashMap<usize, Round1Msg>,
    all_r2: HashMap<usize, Round2Msg>,
    all_r3: HashMap<usize, Round3Msg>,

    result:   Option<Vec<u8>>,
    finished: bool,
}

impl MlDsaProtocol {
    /// Create a new ML-DSA threshold protocol instance.
    ///
    /// All N servers run rounds 1 and 2.  After round 2 the active set
    /// (`act`) is computed deterministically as the first T party-ids
    /// (sorted) that contributed round-2 messages.  Only those T servers
    /// execute round 3; the remaining N-T servers send a no-op Default
    /// message and wait for the T round-3 messages before finalizing.
    pub fn new(
        private_key: Arc<PrivateKeyShare>,
        msg: Vec<u8>,
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

        let n_parties   = party_key.get_n() as usize;
        let threshold   = party_key.get_threshold() as usize;
        let my_party_id = party_key.get_share_id() as usize;

        Ok(Self {
            round: 1,
            party_key,
            tp,
            vk,
            msg,
            signing_seed,
            n_parties,
            threshold,
            my_party_id,
            act: Vec::new(), // computed after round 2
            own_state: None,
            own_r1: None,
            all_r1: HashMap::new(),
            all_r2: HashMap::new(),
            all_r3: HashMap::new(),
            result: None,
            finished: false,
        })
    }

    /// Deterministically choose the active set from the parties that
    /// contributed round-2 messages: sort their ids and take the first T.
    fn compute_active_set(&self) -> Vec<usize> {
        let mut ids: Vec<usize> = self.all_r2.keys().copied().collect();
        ids.sort();
        ids.truncate(self.threshold);
        ids
    }
}

impl ThresholdRoundProtocol<NetMessage> for MlDsaProtocol {
    type ProtocolMessage = MlDsaMessage;

    fn do_round(&mut self) -> Result<Self::ProtocolMessage, ProtocolError> {
        match self.round {
            // ── Round 1: Commit ─────────────────────────────────────
            // Every server participates.
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

            // ── Round 2: Reveal ─────────────────────────────────────
            // Every server participates.  We pass a placeholder `act`
            // to share_sign_2 -- it only stores it in state and does not
            // use it.  The real `act` is written before round 3.
            2 => {
                let state = self.own_state.as_mut()
                    .ok_or(ProtocolError::InternalError)?;
                let r1_msgs: Vec<Round1Msg> = self.all_r1.values()
                    .map(|m| m.clone_msg())
                    .collect();
                // Placeholder act -- overwritten before round 3.
                let placeholder_act: Vec<usize> = (0..self.n_parties).collect();
                let r2_msg = share_sign_2(
                    state, &self.vk, &placeholder_act, &self.msg, r1_msgs,
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

            // ── Round 3: Respond ────────────────────────────────────
            // Only the T servers in `act` produce a real message.
            // The remaining N-T servers send Default (a no-op).
            3 => {
                // Compute the active set deterministically.
                self.act = self.compute_active_set();
                info!(
                    "MlDsaProtocol: party {} -- active set = {:?}",
                    self.my_party_id, self.act
                );

                // Write the real act into the core signing state so that
                // share_sign_3 -> rss_recover sees the correct T-sized set.
                if let Some(state) = self.own_state.as_mut() {
                    state.set_act(self.act.clone());
                }

                // If this server is NOT in the active set, skip round 3.
                if !self.act.contains(&self.my_party_id) {
                    info!(
                        "MlDsaProtocol: party {} not in active set, sending Default for round 3",
                        self.my_party_id
                    );
                    self.round = 4;
                    return Ok(MlDsaMessage::Default);
                }

                // This server IS active -- execute share_sign_3.
                let state = self.own_state.as_ref()
                    .ok_or(ProtocolError::InternalError)?;
                let pk = self.party_key.to_party_key()
                    .map_err(|_| ProtocolError::InternalError)?;

                // Only pass round-2 messages from active-set members.
                let r2_msgs: Vec<Round2Msg> = self.act.iter()
                    .filter_map(|id| self.all_r2.get(id))
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
            // Wait for ALL N round-1 messages so every server agrees on
            // the same set of commitments before revealing.
            2 => self.all_r1.len() >= self.n_parties,
            // Wait for ALL N round-2 messages so the active-set
            // computation is deterministic across all servers.
            3 => self.all_r2.len() >= self.n_parties,
            _ => false,
        }
    }

    fn is_ready_to_finalize(&self) -> bool {
        // Need exactly T round-3 messages (only active servers send them).
        self.round == 4 && self.all_r3.len() >= self.threshold
    }

    fn finalize(&mut self) -> Result<Vec<u8>, ProtocolError> {
        // Only use round-2 data from the active set.
        let r2_msgs: Vec<Round2Msg> = self.act.iter()
            .filter_map(|id| self.all_r2.get(id))
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

        let scheme = self.party_key.get_scheme();
        let wrapped = match scheme {
            ThresholdScheme::MlDsa44 => Signature::MlDsa44(MlDsaSignature { bytes: sig_bytes }),
            ThresholdScheme::MlDsa65 => Signature::MlDsa65(MlDsaSignature { bytes: sig_bytes }),
            ThresholdScheme::MlDsa87 => Signature::MlDsa87(MlDsaSignature { bytes: sig_bytes }),
            _ => return Err(ProtocolError::InternalError),
        };
        let result = wrapped.to_bytes()
            .map_err(|_| ProtocolError::InternalError)?;
        self.result = Some(result.clone());
        self.finished = true;
        info!("MlDsaProtocol: signature assembled successfully.");
        Ok(result)
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