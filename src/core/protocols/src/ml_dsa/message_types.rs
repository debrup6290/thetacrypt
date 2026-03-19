use serde::{Deserialize, Serialize};
use theta_network::types::message::{Channel, NetMessage, NetMessageMetadata};
use crate::interface::{ProtocolError, ProtocolMessageWrapper};

#[derive(Serialize, Deserialize, Clone)]
pub struct MlDsaR1Payload {
    pub party_id: usize,
    pub commitments: Vec<Vec<u8>>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct MlDsaR2Payload {
    pub party_id: usize,
    pub w_values_bytes: Vec<u8>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct MlDsaR3Payload {
    pub party_id: usize,
    pub z1_values_bytes: Vec<u8>,
}

#[derive(Serialize, Deserialize)]
pub enum MlDsaMessage {
    Round1(MlDsaR1Payload),
    Round2(MlDsaR2Payload),
    Round3(MlDsaR3Payload),
    Default,
}

impl Default for MlDsaMessage {
    fn default() -> Self { MlDsaMessage::Default }
}

impl ProtocolMessageWrapper<NetMessage> for MlDsaMessage {
    fn unwrap(wrapped: NetMessage) -> Result<Box<Self>, ProtocolError> {
        let bytes = wrapped.get_message_data().to_owned();
        let result = serde_json::from_str::<MlDsaMessage>(
            &String::from_utf8(bytes).map_err(|_| ProtocolError::InternalError)?
        );
        match result {
            Ok(msg) => Ok(Box::new(msg)),
            Err(_)  => Err(ProtocolError::InternalError),
        }
    }

    fn wrap(&self, instance_id: &String) -> Result<NetMessage, String> {
        let data = serde_json::to_string(self)
            .map_err(|e| e.to_string())?
            .into_bytes();
        let metadata = NetMessageMetadata::new(Channel::Gossip);
        Ok(NetMessage::new(instance_id.clone(), metadata, data))
    }

    fn is_default(&self) -> bool {
        matches!(self, MlDsaMessage::Default)
    }
}
