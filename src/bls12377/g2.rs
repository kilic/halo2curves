use super::fq::Fq;
use super::fq2::Fq2;
use super::fr::Fr;
use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::ff_ext::ExtField;
use crate::group::Curve;
use crate::group::{prime::PrimeCurveAffine, Group, GroupEncoding};
use crate::serde::{Compressed, CompressedFlagConfig};
use crate::{
    impl_binops_additive, impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, new_curve_impl,
};
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

const G2_GENERATOR_X: Fq2 = Fq2 {
    // 0x018480be71c785fec89630a2a3841d01c565f071203e50317ea501f557db6b9b71889f52bb53540274e3e48f7c005196
    c0: Fq([
        0x68904082f268725b,
        0x668f2ea74f45328b,
        0xebca7a65802be84f,
        0x1e1850f4c1ada3e6,
        0x830dc22d588ef1e9,
        0x01862a81767c0982,
    ]),
    // 0x00ea6040e700403170dc5a51b1b140d5532777ee6651cecbe7223ece0799c9de5cf89984bff76fe6b26bfefa6ea16afe
    c1: Fq([
        0x5f02a915c91c7f39,
        0xf8c553ba388da2a7,
        0xd51a416dbd198850,
        0xe943c6f38ae3073a,
        0xffe24aa8259a4981,
        0x011853391e73dfdd,
    ]),
};

const G2_GENERATOR_Y: Fq2 = Fq2 {
    // 0x00690d665d446f7bd960736bcbb2efb4de03ed7274b49a58e458c282f832d204f2cf88886d8c7c2ef094094409fd4ddf
    c0: Fq([
        0xd5b19b897881430f,
        0x05be9118a5b371ed,
        0x6063f91f86c131ee,
        0x3244a61be8f4ec19,
        0xa02e425b9f9a3a12,
        0x018af8c04f3360d2,
    ]),
    // 0x00f8169fd28355189e549da3151a70aa61ef11ac3d591bf12463b01acee304c24279b83f5e52270bd9a1cdd185eb8f93
    c1: Fq([
        0x57601ac71a5b96f5,
        0xe99acc1714f2440e,
        0x2339612f10118ea9,
        0x8321e68a3b1cd722,
        0x2b543b050cc74917,
        0x00590182b396c112,
    ]),
};

const G2_B: Fq2 = Fq2 {
    c0: Fq::ZERO,
    // 010222f6db0fd6f343bd03737460c589dc7b4f91cd5fd889129207b63c6bf8000dd39e5c1ccccccd1c9ed9999999999a
    c1: Fq([
        0x8072266666666685,
        0x8df55926899999a9,
        0x7fe4561ad64f34cf,
        0xb95da6d8b6e4f01b,
        0x4b747cccfc142743,
        0x0039c3fa70f49f43,
    ]),
};

const G2_A: Fq2 = Fq2::ZERO;

new_curve_impl!(
    (pub),
    G2,
    G2Affine,
    Fq2,
    Fr,
    (G2_GENERATOR_X, G2_GENERATOR_Y),
    G2_A,
    G2_B,
    "bls12377_g2",
    |_| unimplemented!(),
    crate::serde::CompressedFlagConfig::ThreeSpare

);

impl crate::serde::endian::EndianRepr for Fq2 {
    const ENDIAN: crate::serde::endian::Endian = Fq::ENDIAN;

    fn to_bytes(&self) -> Vec<u8> {
        self.to_bytes().to_vec()
    }

    fn from_bytes(bytes: &[u8]) -> subtle::CtOption<Self> {
        Fq2::from_bytes(bytes[..Fq2::SIZE].try_into().unwrap())
    }
}

impl Compressed<G2Affine> for G2Compressed {
    const CONFIG: CompressedFlagConfig = CompressedFlagConfig::ThreeSpare;
    fn sign(c: &G2Affine) -> Choice {
        c.y.lexicographically_largest() & !c.is_identity()
    }
    fn resolve(x: Fq2, sign_set: Choice) -> CtOption<G2Affine> {
        G2Affine::y2(x).sqrt().map(|y| {
            let y = Fq2::conditional_select(&y, &-y, sign_set ^ y.lexicographically_largest());
            G2Affine { x, y }
        })
    }
}

impl group::cofactor::CofactorGroup for G2 {
    type Subgroup = G2;

    /// Clears the cofactor, using [Budroni-Pintore](https://ia.cr/2017/419).
    /// This is equivalent to multiplying by $h\_\textrm{eff} = 3(z^2 - 1) \cdot
    /// h_2$, where $h_2$ is the cofactor of $\mathbb{G}\_2$ and $z$ is the
    /// parameter of BLS12-381.
    fn clear_cofactor(&self) -> G2 {
        let t1 = self.mul_by_x(); // [x] P
        let t2 = self.psi(); // psi(P)

        self.double().psi2() // psi^2(2P)
            + (t1 + t2).mul_by_x() // psi^2(2P) + [x^2] P + [x] psi(P)
            - t1 // psi^2(2P) + [x^2 - x] P + [x] psi(P)
            - t2 // psi^2(2P) + [x^2 - x] P + [x - 1] psi(P)
            - self // psi^2(2P) + [x^2 - x - 1] P + [x - 1] psi(P)
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    /// Returns true if this point is free of an $h$-torsion component, and so it
    /// exists within the $q$-order subgroup $\mathbb{G}_2$. This should always return true
    /// unless an "unchecked" API was used.
    fn is_torsion_free(&self) -> Choice {
        // Algorithm from Section 4 of https://eprint.iacr.org/2021/1130
        // Updated proof of correctness in https://eprint.iacr.org/2022/352
        //
        // Check that psi(P) == [x] P
        self.psi().ct_eq(&self.mul_by_x())
    }
}

impl G2 {
    /// Multiply `self` by `crate::BLS_X`, using double and add.

    fn mul_by_x(&self) -> G2 {
        let mut acc = G2::identity();
        for (i, b) in super::BLS_X.into_iter().enumerate() {
            (i != 0).then(|| acc = acc.double());
            (b == 1).then(|| acc += self);
        }
        acc
    }

    fn psi(&self) -> G2 {
        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fq2 {
            c0: Fq([
                0x5892506da58478da,
                0x133366940ac2a74b,
                0x9b64a150cdf726cf,
                0x5cc426090a9c587e,
                0x5cf848adfdcd640c,
                0x004702bf3ac02380,
            ]),
            c1: Fq::ZERO,
        };
        // 1 / ((u+1) ^ (p-1)/2)
        let psi_coeff_y = Fq2 {
            c0: Fq([
                0x982c13d9d084771f,
                0xfd49de0c6da34a32,
                0x61a530d183ab0e53,
                0xdf8fe44106dd9879,
                0x40f29b58d88472bc,
                0x0158723199046d5d,
            ]),
            c1: Fq::ZERO,
        };

        // x = frobenius(x)/((u+1)^((p-1)/3))
        let mut x = self.x;
        x.frobenius_map(1);
        x.mul_assign(&psi_coeff_x);

        // y = frobenius(y)/(u+1)^((p-1)/2)
        let mut y = self.y;
        y.frobenius_map(1);
        y.mul_assign(&psi_coeff_y);

        // z = frobenius(z)
        let mut z = self.z;
        z.frobenius_map(1);

        G2 { x, y, z }
    }

    fn psi2(&self) -> G2 {
        G2 {
            // x = frobenius^2(x)/2^((p-1)/3); note that q^2 is the order of the field.
            x: self.x.mul_by_base(&Fq::ZETA),
            // y = -frobenius^2(y); note that q^2 is the order of the field.
            y: self.y.neg(),
            // z = z
            z: self.z,
        }
    }
}
