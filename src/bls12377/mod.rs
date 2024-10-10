mod engine;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;
mod g1;
mod g2;

pub use engine::*;
pub use fq::*;
pub use fq12::*;
pub use fq2::*;
pub use fq6::*;
pub use fr::*;
pub use g1::*;
pub use g2::*;

// x = 0x8508c00000000001
const BLS_X: [u8; 64] = [
    1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
];
