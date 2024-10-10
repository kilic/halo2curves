use super::fq::Fq;
use super::fq2::Fq2;
use crate::ff_ext::{
    cubic::{CubicExtField, CubicExtFieldArith, CubicSparseMul},
    ExtField,
};
use ff::Field;

crate::impl_binops_additive!(Fq6, Fq6);
crate::impl_binops_multiplicative!(Fq6, Fq6);
crate::impl_binops_calls!(Fq6);
crate::impl_sum_prod!(Fq6);
pub type Fq6 = CubicExtField<Fq2>;

impl CubicExtFieldArith for Fq6 {
    type Base = Fq2;
}

impl CubicSparseMul for Fq6 {
    type Base = Fq2;
}

impl ExtField for Fq6 {
    // (0, 1, 0)
    const NON_RESIDUE: Self = Fq6::new(Fq2::ZERO, Fq2::ONE, Fq2::ZERO);

    fn frobenius_map(&mut self, power: usize) {
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        self.c2.frobenius_map(power);
        self.c1.mul_assign(&FROBENIUS_COEFF_FQ6_C1[power % 6]);
        self.c2.mul_assign(&FROBENIUS_COEFF_FQ6_C2[power % 6]);
    }

    fn mul_by_nonresidue(self: &Fq6) -> Fq6 {
        let c0 = self.c2.mul_by_nonresidue();
        let c1 = self.c0;
        let c2 = self.c1;
        Self { c0, c1, c2 }
    }
}

pub const FROBENIUS_COEFF_FQ6_C1: [Fq2; 6] = [
    // z ^ (( p ^ 0 - 1) / 3)
    Fq2 {
        c0: Fq([
            0x2cdffffffffff68,
            0x51409f837fffffb1,
            0x9f7db3a98a7d3ff2,
            0x7b4e97b76e7c6305,
            0x4cf495bf803c84e8,
            0x008d6661e2fdf49a,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( p ^ 1 - 1) / 3)
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
    // z ^ (( p ^ 2 - 1) / 3)
    Fq2 {
        c0: Fq([
            0xdacd106da5847973,
            0xd8fe2454bac2a79a,
            0x1ada4fd6fd832edc,
            0xfb9868449d150908,
            0xd63eb8aeea32285e,
            0x0167d6a36f873fd0,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( p ^ 3 - 1) / 3)
    Fq2 {
        c0: Fq([
            0x823ac00000000099,
            0xc5cabdc0b000004f,
            0x7f75ae862f8c080d,
            0x9ed4423b9278b089,
            0x79467000ec64c452,
            0x0120d3e434c71c50,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( p ^ 4 - 1) / 3)
    Fq2 {
        c0: Fq([
            0x2c766f925a7b8727,
            0x03d7f6b0253d58b5,
            0x838ec0deec122131,
            0xbd5eb3e9f658bb10,
            0x6942bd126ed3e52e,
            0x1673786dd04ed6a,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( p ^ 5 - 1) / 3)
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
];

// z = u + 1
pub const FROBENIUS_COEFF_FQ6_C2: [Fq2; 6] = [
    // z ^ (( 2 * p ^ 0 - 2) / 3)
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
    // z ^ (( 2 * p ^ 1 - 2) / 3)
    Fq2 {
        c0: Fq([
            0xdacd106da5847973,
            0xd8fe2454bac2a79a,
            0x1ada4fd6fd832edc,
            0xfb9868449d150908,
            0xd63eb8aeea32285e,
            0x0167d6a36f873fd0,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( 2 * p ^ 2 - 2) / 3)
    Fq2 {
        c0: Fq([
            0x2c766f925a7b8727,
            0x03d7f6b0253d58b5,
            0x838ec0deec122131,
            0xbd5eb3e9f658bb10,
            0x6942bd126ed3e52e,
            0x01673786dd04ed6a,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( 2 * p ^ 3 - 2) / 3)
    Fq2 {
        c0: Fq([
            0x2cdffffffffff68,
            0x51409f837fffffb1,
            0x9f7db3a98a7d3ff2,
            0x7b4e97b76e7c6305,
            0x4cf495bf803c84e8,
            0x008d6661e2fdf49a,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( 2 * p ^ 4 - 2) / 3)
    Fq2 {
        c0: Fq([
            0xdacd106da5847973,
            0xd8fe2454bac2a79a,
            0x1ada4fd6fd832edc,
            0xfb9868449d150908,
            0xd63eb8aeea32285e,
            0x0167d6a36f873fd0,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
    // z ^ (( 2 * p ^ 5 - 2) / 3)
    Fq2 {
        c0: Fq([
            0x2c766f925a7b8727,
            0x03d7f6b0253d58b5,
            0x838ec0deec122131,
            0xbd5eb3e9f658bb10,
            0x6942bd126ed3e52e,
            0x01673786dd04ed6a,
        ]),
        c1: Fq([0, 0, 0, 0, 0, 0]),
    },
];

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fq6, "field_arithmetic");
    // extension field-specific
    crate::field_testing_suite!(Fq6, "cubic_sparse_mul", Fq2);
    crate::field_testing_suite!(
        Fq6,
        "frobenius",
        // Frobenius endomorphism power parameter for extension field
        //  ϕ: E → E
        //  (x, y) ↦ (x^p, y^p)
        // p: modulus of base field (Here, Fq::MODULUS)
        Fq::MODULUS_LIMBS
    );

    #[test]
    fn test_fq6_mul_nonresidue() {
        use ff::Field;
        let e = Fq6::random(rand_core::OsRng);
        let a0 = e.mul_by_nonresidue();
        let a1 = e * Fq6::NON_RESIDUE;
        assert_eq!(a0, a1);
    }
}
