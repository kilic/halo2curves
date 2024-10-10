use super::fq::Fq;
use crate::ff_ext::{
    cubic::{CubicExtField, CubicExtFieldArith, CubicSparseMul},
    ExtField,
};
use ff::Field;
use std::ops::MulAssign;

crate::impl_binops_additive!(Fq3, Fq3);
crate::impl_binops_multiplicative!(Fq3, Fq3);
crate::impl_binops_calls!(Fq3);
crate::impl_sum_prod!(Fq3);
pub type Fq3 = CubicExtField<Fq>;

impl CubicExtFieldArith for Fq3 {
    type Base = Fq;
}

impl CubicSparseMul for Fq3 {
    type Base = Fq;
}

impl ExtField for Fq3 {
    // (0, 1, 0)
    const NON_RESIDUE: Self = Fq3::new(Fq::ZERO, Fq::ONE, Fq::ZERO);

    fn frobenius_map(&mut self, power: usize) {
        self.c1.mul_assign(&FROBENIUS_COEFF_FQ3_C1[power % 3]);
        self.c2.mul_assign(&FROBENIUS_COEFF_FQ3_C2[power % 3]);
    }

    fn mul_by_nonresidue(self: &Fq3) -> Fq3 {
        let c0 = self.c2.mul_by_nonresidue();
        let c1 = self.c0;
        let c2 = self.c1;
        Self { c0, c1, c2 }
    }
}

pub const FROBENIUS_COEFF_FQ3_C1: [Fq; 3] = [
    Fq::ONE,
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
];

pub const FROBENIUS_COEFF_FQ3_C2: [Fq; 3] = [
    Fq::ONE,
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
];

#[cfg(test)]
mod test {
    use super::*;
    use crate::{arith_test, frobenius_test, setup_f6_test_funcs, test};
    use rand_core::RngCore;

    macro_rules! test_fq6 {
        ($test:ident, $size: expr) => {
            paste::paste! {
            #[test]
            fn [< $test test >]() {
                use rand::SeedableRng;
                use rand_xorshift::XorShiftRng;
                let mut rng = XorShiftRng::from_seed(crate::tests::SEED);
                crate::bw6::Fq3::test::$test(&mut rng, $size);
            }
            }
        };
    }

    arith_test!(Fq3);
    setup_f6_test_funcs!(Fq3, Fq);
    // test_fq6!(f6_mul_nonresidue_, 1000);
    // test_fq6!(f6_mul_by_1_, 1000);
    // test_fq6!(f6_mul_by_01_, 1000);
    frobenius_test!(Fq3, Fq, 10);
}
