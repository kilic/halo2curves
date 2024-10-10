use super::fq::Fq;
use super::fq3::Fq3;
use crate::ff_ext::{
    quadratic::{QuadExtField, QuadExtFieldArith, QuadSparseMul},
    ExtField,
};

pub type Fq6 = QuadExtField<Fq3>;

impl QuadExtFieldArith for Fq6 {
    type Base = Fq3;
}

impl QuadSparseMul for Fq6 {
    type Base = Fq;
}

impl ExtField for Fq6 {
    const NON_RESIDUE: Self = Fq6::zero(); // no needs

    fn frobenius_map(&mut self, power: usize) {
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        self.c1
            .mul_assign_by_base(&FROBENIUS_COEFF_FQ6_C1[power % 6]);
    }
}

crate::impl_binops_additive!(Fq6, Fq6);
crate::impl_binops_multiplicative!(Fq6, Fq6);
crate::impl_binops_calls!(Fq6);
crate::impl_sum_prod!(Fq6);
crate::impl_cyclotomic_square!(Fq, Fq6);

pub const FROBENIUS_COEFF_FQ6_C1: [Fq; 6] = [
    Fq::one(),
    Fq([
        0x8cfcb51bd8404a93,
        0x495e69d68495a383,
        0xd23cbc9234705263,
        0x8d2b4c2b5fcf4f52,
        0x6a798a5d20c612ce,
        0x3e825d90eb6c2443,
        0x772b249f2c9525fe,
        0x521b2ed366e4b9bb,
        0x84abb49bd7c4471d,
        0x907062359c0f17e3,
        0x3385e55030cc6f12,
        0x3f11a3a41a2606,
    ]),
    Fq([
        0x7f96b51bd840c549,
        0xd59782096496171f,
        0x49b046fd9ce14bbc,
        0x4b6163bba7527a56,
        0xef6c92fb771d59f1,
        0x0425bedbac1dfdc7,
        0xd3ac39de759c0ffd,
        0x9f43ed0e063a81d0,
        0x5bd7d20b4f9a3ce2,
        0x0411f03c36cf5c3c,
        0x2d658fd49661c472,
        0x01100249ae760b93,
    ]),
    Fq([
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
    ]),
    Fq([
        0x67a04ae427bfb5f8,
        0x9d32d491eb6a5cff,
        0x43d03c1cb68051d4,
        0x0b75ca96f69859a5,
        0x0763497f5325ec60,
        0x48076b5c278dd94d,
        0x8ca3965ff91efd06,
        0x1e6077657ea02f5d,
        0xcdd6c153a8c37724,
        0x28b5b634e5c22ea4,
        0x9e01e3efd42e902c,
        0x00e3d6815769a804,
    ]),
    Fq([
        0x75064ae427bf3b42,
        0x10f9bc5f0b69e963,
        0xcc5cb1b14e0f587b,
        0x4d3fb306af152ea1,
        0x827040e0fccea53d,
        0x82640a1166dbffc8,
        0x30228120b0181307,
        0xd137b92adf4a6748,
        0xf6aaa3e430ed815e,
        0xb514282e4b01ea4b,
        0xa422396b6e993acc,
        0x0012e5db4d0dc277,
    ]),
];

// #[cfg(test)]
// mod test {
//     use super::*;
//     crate::field_testing_suite!(Fq6, "field_arithmetic");
//     // extension field-specific
//     crate::field_testing_suite!(Fq6, "quadratic_sparse_mul", Fq3, Fq);
//     crate::field_testing_suite!(
//         Fq6,
//         "frobenius",
//         // Frobenius endomorphism power parameter for extension field
//         //  ϕ: E → E
//         //  (x, y) ↦ (x^p, y^p)
//         // p: modulus of base field (Here, Fq::MODULUS)
//         Fq::MODULUS_LIMBS
//     );
// }

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
                crate::bw6::fq6::test::$test(&mut rng, $size);
            }
            }
        };
    }
    use super::*;
    use crate::{arith_test, frobenius_test, setup_f12_test_funcs, test};
    use ff::Field;
    use rand::RngCore;

    arith_test!(Fq6);

    // F12 specific
    setup_f12_test_funcs!(Fq6, Fq3, Fq);
    test_fq12!(f12_mul_by_014_, 500);
    test_fq12!(f12_mul_by_034_, 500);
    frobenius_test!(Fq6, Fq, 8);
}
