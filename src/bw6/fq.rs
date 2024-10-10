use core::convert::TryInto;
use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

impl_field!(
    bw6_base,
    Fq,
    modulus = "122e824fb83ce0ad187c94004faff3eb926186a81d14688528275ef8087be41707ba638e584e91903cebaff25b423048689c8ed12f9fd9071dcd3dc73ebff2e98a116c25667a8f8160cf8aeeaf0a437e6913e6870000082f49d00000000008b",
    mul_gen = "2",
    zeta = "531dc16c6ecd27aa846c61024e4cca6c1f31e53bd9603c2d17be416c5e4426ee4a737f73b6f952ab5e57926fa701848e0a235a0a398300c65759fc45183151f2f082d4dcb5e37cb6290012d96f8819c547ba8a4000002f962140000000002a",
    from_uniform = [192],
    endian = "little",
);

crate::extend_field_legendre!(Fq);
crate::impl_binops_calls!(Fq);
crate::impl_binops_additive!(Fq, Fq);
crate::impl_binops_multiplicative!(Fq, Fq);
crate::field_bits!(Fq);
crate::serialize_deserialize_primefield!(Fq);
crate::impl_from_u64!(Fq);
use crate::ff_ext::ExtField;

impl ExtField for Fq {
    // -4
    const NON_RESIDUE: Self = Fq([
        0xe12e00000001e9c2,
        0x63c1e3faa001cd69,
        0xb1b4384fcbe29cf6,
        0xc79630bc713d5a1d,
        0x30127ac071851e2d,
        0x0979f350dcd36af1,
        0x6a66defed8b361f2,
        0x53abac78b24d4e23,
        0xb7ab89dede485a92,
        0x5c3a0745675e8452,
        0x446f17918c5f5700,
        0x00fdf24e3267fa1e,
    ]);
    fn mul_by_nonresidue(&self) -> Self {
        self.double().double().neg()
    }
    fn frobenius_map(&mut self, _: usize) {}
}

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fq, "field_arithmetic");
    crate::field_testing_suite!(Fq, "conversion");
    crate::field_testing_suite!(Fq, "serialization");
    crate::field_testing_suite!(Fq, "quadratic_residue");
    crate::field_testing_suite!(Fq, "bits");
    crate::field_testing_suite!(Fq, "serialization_check");
    crate::field_testing_suite!(Fq, "constants");
    crate::field_testing_suite!(Fq, "sqrt");
    crate::field_testing_suite!(Fq, "zeta");
    crate::field_testing_suite!(Fq, "from_uniform_bytes", 192);
    #[test]
    fn test_fq_mul_nonresidue() {
        let e = Fq::random(rand_core::OsRng);
        let a0 = e.mul_by_nonresidue();
        let a1 = e * Fq::NON_RESIDUE;
        assert_eq!(a0, a1);
    }
}
