use core::convert::TryInto;
use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

impl_field!(
    bls12377_scalar,
    Fr,
    modulus = "12ab655e9a2ca55660b44d1e5c37b00159aa76fed00000010a11800000000001",
    mul_gen = "16",
    zeta = "452217cc900000010a11800000000000",
    from_uniform = [64],
    endian = "little",
);

crate::extend_field_legendre!(Fr);
crate::impl_binops_calls!(Fr);
crate::impl_binops_additive!(Fr, Fr);
crate::impl_binops_multiplicative!(Fr, Fr);
crate::field_bits!(Fr);
crate::serialize_deserialize_primefield!(Fr);
crate::impl_from_u64!(Fr);

#[cfg(test)]
mod test {
    use super::Fr;
    use crate::{
        arith_test, constants_test, from_uniform_bytes_test, legendre_test, serde_test, test,
    };

    constants_test!(Fr);

    arith_test!(Fr);
    legendre_test!(Fr);
    test!(arith, Fr, sqrt_test, 1000);

    serde_test!(Fr PrimeFieldBits);
    from_uniform_bytes_test!(Fr, 1000, L 64);
}
