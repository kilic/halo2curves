use core::convert::TryInto;
use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// zeta = 01ae3a4617c510eabc8756ba8f8c524eb8882a75cc9bc8e359064ee822fb5bffd1e945779fffffffffffffffffffffff
impl_field!(
    bls12377_base,
    Fq,
    modulus = "01ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001",
    mul_gen = "f",
    zeta = "9b3af05dd14f6ec619aaf7d34594aabc5ed1347970dec00452217cc900000008508c00000000001",
    from_uniform = [64, 96],
    endian = "big",
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
    // -5
    #[rustfmt::skip]
    const NON_RESIDUE: Self = Fq::from_raw([0x8508bffffffffffc, 0x170b5d4430000000, 0x1ef3622fba094800, 0x1a22d9f300f5138f, 0xc63b05c06ca1493b, 0x01ae3a4617c510ea]);
    fn mul_by_nonresidue(&self) -> Self {
        (self.double().double() + self).neg()
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
    crate::field_testing_suite!(Fq, "from_uniform_bytes", 64, 96);
    #[test]
    fn test_fq_mul_nonresidue() {
        let e = Fq::random(rand_core::OsRng);
        let a0 = e.mul_by_nonresidue();
        let a1 = e * Fq::NON_RESIDUE;
        assert_eq!(a0, a1);
    }
}
