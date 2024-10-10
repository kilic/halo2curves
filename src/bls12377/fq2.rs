use super::fq::Fq;
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::quadratic::{QuadExtField, QuadExtFieldArith, SQRT};
use crate::ff_ext::{ExtField, Legendre};
use core::convert::TryInto;
use std::cmp::Ordering;
use subtle::{Choice, CtOption};

crate::impl_binops_additive!(Fq2, Fq2);
crate::impl_binops_multiplicative!(Fq2, Fq2);
crate::impl_binops_calls!(Fq2);
crate::impl_sum_prod!(Fq2);
crate::impl_tower2!(Fq, Fq2);
crate::impl_tower2_from_uniform_bytes!(Fq, Fq2, 128);

pub type Fq2 = QuadExtField<Fq>;
impl QuadExtFieldArith for Fq2 {
    type Base = Fq;
    const SQRT: SQRT<Fq> = SQRT::Algorithm10 {
        precompute_e: Fq2 {
            c0: Fq::ZERO,
            c1: Fq([
                0xb808d0c529b417d3,
                0x990ec6027c53dba3,
                0xe0543cf6b3ef02dd,
                0x2cb660d9ce2c51b1,
                0x403085782b4db08c,
                0x0044e3a38500e2ac,
            ]),
        },
        precompute_f: Fq2 {
            c0: Fq([
                0x88fd3ffffffffd07,
                0x7f37c04d4ffffe74,
                0xfe81201ffa68f7bb,
                0x4e661ca22778db8c,
                0xba8be6fd148d4f4f,
                0x0114c5a35730b618,
            ]),
            c1: Fq::ZERO,
        },
        q_minus_1_over_4: &[
            0x2142300000000000,
            0x05c2d7510c000000,
            0xc7bcd88bee825200,
            0xc688b67cc03d44e3,
            0xb18ec1701b28524e,
            0x006b8e9185f1443a,
        ],
    };

    fn square_assign(el: &mut QuadExtField<Self::Base>) {
        let t0 = el.c0 - el.c1;

        let t1 = el.c1.double(); // 2c1
        let t1 = t1.double() + t1; // 6c1
        let t1 = t0 + t1; // c0 + 5c1
        let t1 = t0 * t1; // (c0 - c1) (c0 + 5c1)

        el.c1 = (el.c0 * el.c1).double(); // 2c0c1
        el.c0 = t1 - el.c1.double();
    }
}

impl ExtField for Fq2 {
    // (0, 1)
    const NON_RESIDUE: Self = Fq2::new(Fq::ZERO, Fq::ONE);

    fn mul_by_nonresidue(&self) -> Self {
        Self {
            c0: self.c1.mul_by_nonresidue(),
            c1: self.c0,
        }
    }

    fn frobenius_map(&mut self, power: usize) {
        if power % 2 != 0 {
            self.conjugate();
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::{
        arith_test, constants_test, f2_test, frobenius_test, legendre_test, serde_test, test,
    };
    use rand_core::RngCore;

    constants_test!(Fq2);

    arith_test!(Fq2);
    legendre_test!(Fq2);
    test!(arith, Fq2, sqrt_test, 1000);

    serde_test!(Fq2);

    f2_test!(Fq2, Fq);
    frobenius_test!(Fq2, Fq, 20);

    #[test]
    fn test_fq2_mul_nonresidue() {
        let e = Fq2::random(rand_core::OsRng);
        let a0 = e.mul_by_nonresidue();
        let a1 = e * Fq2::NON_RESIDUE;
        assert_eq!(a0, a1);
    }
}
