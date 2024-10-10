use super::fq::Fq;
use super::fq2::Fq2;
use super::fq6::Fq6;
use crate::ff_ext::{
    quadratic::{QuadExtField, QuadExtFieldArith, QuadSparseMul},
    ExtField,
};
use ff::Field;

pub type Fq12 = QuadExtField<Fq6>;

impl QuadExtFieldArith for Fq12 {
    type Base = Fq6;
}

impl QuadSparseMul for Fq12 {
    type Base = Fq2;
}

impl ExtField for Fq12 {
    const NON_RESIDUE: Self = Fq12::zero(); // no needs

    fn frobenius_map(&mut self, power: usize) {
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        self.c1
            .mul_assign_by_base(&FROBENIUS_COEFF_FQ12_C1[power % 12]);
    }
}

crate::impl_binops_additive!(Fq12, Fq12);
crate::impl_binops_multiplicative!(Fq12, Fq12);
crate::impl_binops_calls!(Fq12);
crate::impl_sum_prod!(Fq12);
crate::impl_cyclotomic_square!(Fq2, Fq12);

pub const FROBENIUS_COEFF_FQ12_C1: [Fq2; 12] = [
    // z ^ ((p ^ 0 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x02cdffffffffff68,
            0x51409f837fffffb1,
            0x9f7db3a98a7d3ff2,
            0x7b4e97b76e7c6305,
            0x4cf495bf803c84e8,
            0x008d6661e2fdf49a,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 1 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x6ec47a04a3f7ca9e,
            0xa42e0cb968c1fa44,
            0x578d5187fbd2bd23,
            0x930eeb0ac79dd4bd,
            0xa24883de1e09a9ee,
            0x00daa7058067d46f,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 2 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x5892506da58478da,
            0x133366940ac2a74b,
            0x9b64a150cdf726cf,
            0x5cc426090a9c587e,
            0x5cf848adfdcd640c,
            0x004702bf3ac02380,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 3 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x982c13d9d084771f,
            0xfd49de0c6da34a32,
            0x61a530d183ab0e53,
            0xdf8fe44106dd9879,
            0x40f29b58d88472bc,
            0x0158723199046d5d,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 4 - 1) / 6)
    Fq2 {
        c0: Fq([
            0xdacd106da5847973,
            0xd8fe2454bac2a79a,
            0x1ada4fd6fd832edc,
            0xfb9868449d150908,
            0xd63eb8aeea32285e,
            0x167d6a36f873fd0,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 5 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x296799d52c8cac81,
            0x591bd15304e14fee,
            0xa17df4987d85130,
            0x4c80f9363f3fc3bc,
            0x9eaa177aba7ac8ce,
            0x7dcb2c189c98ed,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 6 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x823ac00000000099,
            0xc5cabdc0b000004f,
            0x7f75ae862f8c080d,
            0x9ed4423b9278b089,
            0x79467000ec64c452,
            0x120d3e434c71c50,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 7 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x164445fb5c083563,
            0x72dd508ac73e05bc,
            0xc76610a7be368adc,
            0x8713eee839573ed1,
            0x23f281e24e979f4c,
            0x0d39340975d3c7b,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 8 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x2c766f925a7b8727,
            0x3d7f6b0253d58b5,
            0x838ec0deec122131,
            0xbd5eb3e9f658bb10,
            0x6942bd126ed3e52e,
            0x01673786dd04ed6a,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 9 - 1) / 6)
    Fq2 {
        c0: Fq([
            0xecdcac262f7b88e2,
            0x19c17f37c25cb5cd,
            0xbd4e315e365e39ac,
            0x3a92f5b1fa177b15,
            0x85486a67941cd67e,
            0x0055c8147ec0a38d,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 10 - 1) / 6)
    Fq2 {
        c0: Fq([
            0xaa3baf925a7b868e,
            0x3e0d38ef753d5865,
            0x4191258bc861923,
            0x1e8a71ae63e00a87,
            0xeffc4d11826f20dc,
            0x004663a2a83dd119,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ ((p ^ 11 - 1) / 6)
    Fq2 {
        c0: Fq([
            0x5ba1262ad3735380,
            0xbdef8bf12b1eb012,
            0x14db82e63230f6cf,
            0xcda1e0bcc1b54fd3,
            0x2790ee45b226806c,
            0x01306f19ff2877fd,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
];

#[cfg(test)]
mod test {
    macro_rules! test_fq12 {
        ($test:ident, $size: expr) => {
            paste::paste! {
            #[test]
            fn [< $test test >]() {
                use rand::SeedableRng;
                use rand_xorshift::XorShiftRng;
                let mut rng = XorShiftRng::from_seed(crate::tests::SEED);
                crate::bls12377::fq12::test::$test(&mut rng, $size);
            }
            }
        };
    }
    use super::*;
    use crate::{arith_test, frobenius_test, setup_f12_test_funcs, test};
    use ff::Field;
    use rand::RngCore;

    arith_test!(Fq12);
    // TODO Compile problems with derive_serde feature
    // serde_test!(fq12);

    // F12 specific
    setup_f12_test_funcs!(Fq12, Fq6, Fq2);
    test_fq12!(f12_mul_by_014_, 500);
    test_fq12!(f12_mul_by_034_, 500);
    frobenius_test!(Fq12, Fq, 8);
}
