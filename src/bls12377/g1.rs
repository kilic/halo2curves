use super::fq::Fq;
use super::Fr;
use crate::serde::{Compressed, CompressedFlagConfig};
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
    (GENERATOR_X, GENERATOR_Y),
    A,
    B,
    "bls12377_g1",
    |_| unimplemented!(),
    crate::serde::CompressedFlagConfig::ThreeSpare
);

impl Compressed<G1Affine> for G1Compressed {
    const CONFIG: CompressedFlagConfig = CompressedFlagConfig::ThreeSpare;
    fn sign(c: &G1Affine) -> Choice {
        c.y.lexicographically_largest() & !c.is_identity()
    }
    fn resolve(x: Fq, sign_set: Choice) -> CtOption<G1Affine> {
        G1Affine::y2(x).sqrt().map(|y| {
            let y = Fq::conditional_select(&y, &-y, sign_set ^ y.lexicographically_largest());
            G1Affine { x, y }
        })
    }
}

impl group::cofactor::CofactorGroup for G1 {
    type Subgroup = G1;

    fn clear_cofactor(&self) -> Self {
        self - self.mul_by_x()
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        self.mul_by_x().mul_by_x().neg().endo().ct_eq(self)
    }
}

// 0x008848defe740a67c8fc6225bf87ff5485951e2caa9d41bb188282c8bd37cb5cd5481512ffcd394eeab9b16eb21be9ef
const GENERATOR_X: Fq = Fq([
    0x260f33b9772451f4,
    0xc54dd773169d5658,
    0x5c1551c469a510dd,
    0x761662e4425e1698,
    0xc97d78cc6f065272,
    0x00a41206b361fd4d,
]);

// 0x01914a69c5102eff1f674f5d30afeec4bd7fb348ca3e52d96d182ad44fb82305c2fe3d3634a9591afd82de55559c8ea6
const GENERATOR_Y: Fq = Fq([
    0x8193961fb8cb81f3,
    0x00638d4c5f44adb8,
    0xfafaf3dad4daf54a,
    0xc27849e2d655cd18,
    0x2ec3ddb401d52814,
    0x007da93326303c71,
]);

const A: Fq = Fq::ZERO;
const B: Fq = Fq::ONE;

impl G1 {
    fn mul_by_x(&self) -> G1 {
        let mut acc = G1::identity();
        for (i, b) in super::BLS_X.into_iter().enumerate() {
            (i != 0).then(|| acc = acc.double());
            (b == 1).then(|| acc += self);
        }
        acc
    }
}
