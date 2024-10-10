use super::fq::Fq;
use super::Fr;
use crate::serde::CompressedFlagConfig;
use crate::{
    impl_binops_additive, impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, new_curve_impl,
};
use core::cmp;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::PrimeField;
use ff::WithSmallOrderMulGroup;
use group::{ff::Field, prime::PrimeCurveAffine, Curve, Group, GroupEncoding};
use pasta_curves::arithmetic::{Coordinates, CurveAffine, CurveExt};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

new_curve_impl!(
    (pub),
    G1,
    G1Affine,
    Fq,
    Fr,
    (G1_X, G1_Y),
    G1_A,
    G1_B,
    "bw6_g1",
    |_| unimplemented!(),
    CompressedFlagConfig::TwoSpare,
    standard_sign
);

new_curve_impl!(
    (pub),
    G2,
    G2Affine,
    Fq,
    Fr,
    (G2_X, G2_Y),
    G2_A,
    G2_B,
    "bw6_g2",
    |_| unimplemented!(),
    CompressedFlagConfig::TwoSpare,
    standard_sign
);

fn mul_by_x<C: CurveExt>(p: &C) -> C {
    let mut acc = C::identity();
    let bits = (0..64).map(|i| (super::X >> i) & 1 == 1).rev();
    for (i, b) in bits.enumerate() {
        (i != 0).then(|| acc = acc.double());
        b.then(|| acc += p);
    }
    acc
}

fn mul_u8<C: CurveExt>(p: &C, e: u8) -> C {
    let mut acc = C::identity();
    let bits = (0..8).map(|i| (e >> i) & 1 == 1).rev();
    for (i, b) in bits.enumerate() {
        (i != 0).then(|| acc = acc.double());
        b.then(|| acc += p);
    }
    acc
}

fn is_torsion_free<C: CurveExt>(p: &C, endo: Box<dyn Fn(&C) -> C>) -> Choice {
    let phi = endo(&p);

    let v = mul_by_x(&phi) - phi;
    let v = mul_by_x(&v);
    let v = mul_by_x(&v) + phi;

    (mul_by_x(p) + p).add(v).is_identity()
}

impl group::cofactor::CofactorGroup for G1 {
    type Subgroup = G1;

    fn clear_cofactor(&self) -> Self {
        let p0 = *self;
        let p1 = mul_by_x(&p0);
        let p2 = mul_by_x(&p1);
        let p3 = mul_by_x(&p2);

        let u1 = mul_u8(&p3, 103);
        let u1 = u1 - mul_u8(&p2, 83);
        let u1 = u1 - mul_u8(&p1, 40);
        let u1 = u1 + mul_u8(&p0, 136);

        let u2 = mul_u8(&p2, 7);
        let u2 = u2 + mul_u8(&p1, 89);
        let u2 = u2 + mul_u8(&p0, 130);

        let u2 = u2.endo();

        u1 + u2
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        let endo = Box::new(|p: &Self| p.endo());
        is_torsion_free(self, endo)
    }
}

impl group::cofactor::CofactorGroup for G2 {
    type Subgroup = G2;

    fn clear_cofactor(&self) -> Self {
        let p0 = *self;
        let p1 = mul_by_x(&p0);
        let p2 = mul_by_x(&p1);
        let p3 = mul_by_x(&p2);

        let u1 = mul_u8(&p3, 103);
        let u1 = u1 - mul_u8(&p2, 83);
        let u1 = u1 - mul_u8(&p1, 143);
        let u1 = u1 + mul_u8(&p0, 27);

        let u2 = mul_u8(&p2, 7);
        let u2 = u2 - mul_u8(&p1, 117);
        let u2 = u2 - mul_u8(&p0, 109);

        let u2 = u2.endo().endo();
        let u2 = u2.endo().endo();
        u1 + u2
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        let endo = Box::new(|p: &Self| p.endo().endo());
        is_torsion_free(self, endo)
    }
}

const G1_X: Fq = Fq([
    0xd6e42d7614c2d770,
    0x4bb886eddbc3fc21,
    0x64648b044098b4d2,
    0x1a585c895a422985,
    0xf1a9ac17cf8685c9,
    0x352785830727aea5,
    0xddf8cb12306266fe,
    0x6913b4bfbc9e949a,
    0x3a4b78d67ba5f6ab,
    0x0f481c06a8d02a04,
    0x91d4e7365c43edac,
    0x00f4d17cd48beca5,
]);

const G1_Y: Fq = Fq([
    0x97e805c4bd16411f,
    0x870d844e1ee6dd08,
    0x1eba7a37cb9eab4d,
    0xd544c4df10b9889a,
    0x8fe37f21a33897be,
    0xe9bf99a43a0885d2,
    0xd7ee0c9e273de139,
    0xaa6a9ec7a38dd791,
    0x8f95d3fcf765da8e,
    0x42326e7db7357c99,
    0xe217e407e218695f,
    0x009d1eb23b7cf684,
]);

const G1_A: Fq = Fq::ZERO;
const G1_B: Fq = Fq([
    0xf29a000000007ab6,
    0x8c391832e000739b,
    0x77738a6b6870f959,
    0xbe36179047832b03,
    0x84f3089e56574722,
    0xc5a3614ac0b1d984,
    0x5c81153f4906e9fe,
    0x4d28be3a9f55c815,
    0xd72c1d6f77d5f5c5,
    0x73a18e069ac04458,
    0xf9dfaa846595555f,
    0x00d0f0a60a5be58c,
]);

const G2_X: Fq = Fq([
    0x3d902a84cd9f4f78,
    0x864e451b8a9c05dd,
    0xc2b3c0d6646c5673,
    0x17a7682def1ecb9d,
    0xbe31a1e0fb768fe3,
    0x4df125e09b92d1a6,
    0x0943fce635b02ee9,
    0xffc8e7ad0605e780,
    0x8165c00a39341e95,
    0x8ccc2ae90a0f094f,
    0x73a8b8cc0ad09e0c,
    0x011027e203edd9f4,
]);

const G2_Y: Fq = Fq([
    0x9a159be4e773f67c,
    0x6b957244aa8f4e6b,
    0xa27b70c9c945a38c,
    0xacb6a09fda11d0ab,
    0x3abbdaa9bb6b1291,
    0xdbdf642af5694c36,
    0xb6360bb9560b369f,
    0xac0bd1e822b8d6da,
    0xfa355d17afe6945f,
    0x8d6a0fc1fbcad35e,
    0x72a63c7874409840,
    0x0114976e5b0db280,
]);

const G2_A: Fq = Fq::ZERO;
const G2_B: Fq = Fq([
    0x136efffffffe16c9,
    0x82cf5a6dcffe3319,
    0x6458c05f1f0e0741,
    0xd10ae605e52a4eda,
    0x41ca591c0266e100,
    0x7d0fd59c3626929f,
    0x9967dc004d00c112,
    0x1ccff9c033379af5,
    0x9ad6ec10a23f63af,
    0x5cec11251a72c235,
    0x8d18b1ae789ba83e,
    0x0024f5d6c91bd3ec,
]);

#[cfg(test)]
mod testg1 {
    use super::*;
    use group::UncompressedEncoding;
    crate::curve_testing_suite!(G1, "clear_cofactor");
    crate::curve_testing_suite!(G1);
    crate::curve_testing_suite!(G1, "endo_consistency");
    crate::curve_testing_suite!(
        G1,
        "constants",
        Fq::MODULUS,
        G1_A,
        G1_B,
        G1_X,
        G1_Y,
        Fr::MODULUS
    );
}

#[cfg(test)]
mod testg2 {
    use super::*;
    use group::UncompressedEncoding;
    crate::curve_testing_suite!(G2, "clear_cofactor");
    crate::curve_testing_suite!(G2);
    crate::curve_testing_suite!(G2, "endo_consistency");
    crate::curve_testing_suite!(
        G2,
        "constants",
        Fq::MODULUS,
        G2_A,
        G2_B,
        G2_X,
        G2_Y,
        Fr::MODULUS
    );
}
