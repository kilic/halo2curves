use super::{Fq, Fq6, Fr, G1Affine, G2Affine, G1, G2};
use crate::ff_ext::quadratic::QuadSparseMul;
use crate::ff_ext::ExtField;
use core::borrow::Borrow;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{Field, PrimeField};
use group::prime::PrimeCurveAffine;
use group::{Curve, Group};
use pairing::{Engine, MillerLoopResult, MultiMillerLoop, PairingCurveAffine};
use pasta_curves::arithmetic::CurveExt;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

crate::impl_gt!(Gt, Fq6, Fr);
crate::impl_miller_loop_components!(BW6, G1, G1Affine, G2, G2Affine, Fq6, Gt, Fr);

fn ell(f: &mut Fq6, coeffs: &(Fq, Fq, Fq), p: &G1Affine) {
    let c0 = coeffs.0 * &p.y;
    let c1 = coeffs.1 * &p.x;
    Fq6::mul_by_014(f, &coeffs.2, &c1, &c0);
}

const X: u64 = 0x8508c00000000001;

impl MillerLoopResult for Fq6 {
    type Gt = Gt;

    fn final_exponentiation(&self) -> Gt {
        // Algorithm 4.4 from https://yelhousni.github.io/phd.pdf
        // Reference implementation 1:
        // https://github.com/Consensys/gnark-crypto/blob/ed1dc7cab3e2a0e1f6b9ac45680a8fdfb7dee47f/ecc/bw6-761/pairing.go#L62
        // Reference implementation 2:
        // https://github.com/arkworks-rs/algebra/blob/b33df5cce2d54cf4c9248e4b229c7d6708fa9375/ec/src/models/bw6/mod.rs#L247

        fn exp(f: &Fq6, x: u64) -> Fq6 {
            let bits = (0..64 - x.leading_zeros()).map(|i| (x >> i) as u8 & 1);
            let mut acc = Fq6::one();
            for (i, bit) in bits.rev().enumerate() {
                (i != 0).then(|| acc.cyclotomic_square());
                (bit == 1).then(|| acc *= f);
            }
            acc
        }

        // Easy part
        let t0 = self.invert().unwrap();
        let mut t1 = *self;
        t1.conjugate();
        let t1 = t1 * t0;
        let mut t0 = t1;
        t0.frobenius_map(1);
        let mut f = t0 * t1;

        // Hard part
        let mut t = f;
        t.frobenius_map(1);

        let a = exp(&f, X - 1);
        let a = exp(&a, X - 1);
        let a = a * t;
        f.conjugate();
        let b = exp(&a, X + 1) * f;
        let a = a.square() * a;
        let c = exp(&b, (X - 1) / 3);
        let mut d = exp(&c, X - 1);
        let e = exp(&d, X - 1);
        let e = exp(&e, X - 1) * d;
        d.conjugate();
        let mut fc = d * b;
        let g = exp(&e, X + 1) * fc;
        let h = g * c;
        fc.conjugate();
        let i = exp(&(g * d), X + 1) * fc;
        const D1: u64 = 11;
        const D2: u64 = 103;
        let j = exp(&h, D1) * e;
        let k = j.square() * j * b * exp(&i, D2);
        Gt(a * k)
    }
}

// x+1
// Alligned with LOOP_2_NAF
const LOOP_1_NAF: [i8; 190] = [
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

// x^3-x^2-x
// 0x23ed1347970dec008a442f991fffffffffffffffffffffff
const LOOP_2_NAF: [i8; 190] = [
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, -1,
    0, 1, 0, -1, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 1, 0, 1, 0, 0, 0, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0,
    1, 0, 0, 1,
];

pub fn multi_miller_loop(terms: &[(&G1Affine, &G2Affine)]) -> Fq6 {
    // Follows 6th equation at https://hackmd.io/@gnark/BW6-761-changes

    let terms = terms
        .iter()
        .filter_map(|&(p, q)| {
            if bool::from(p.is_identity()) || bool::from(q.is_identity()) {
                None
            } else {
                Some((*p, *q))
            }
        })
        .collect::<Vec<_>>();

    let mut f = Fq6::one();
    let mut r = terms
        .iter()
        .map(|(_, q)| q.to_curve().endo().neg())
        .collect::<Vec<_>>();

    for (x2, x1) in LOOP_2_NAF
        .iter()
        .rev()
        .skip(1)
        .zip(LOOP_1_NAF.iter().rev().skip(1))
    {
        f.square_assign();
        let x = x2 * 3 + x1;

        terms
            .iter()
            .zip(r.iter_mut())
            .for_each(|((p, _), r)| double(&mut f, r, p));

        match x {
            -3 => {
                // q1neg
                terms.iter().zip(r.iter_mut()).for_each(|((p, q), r)| {
                    add(&mut f, r, &q.to_curve().endo().to_affine(), p);
                });
            }
            -1 => {
                // q0neg
                terms.iter().zip(r.iter_mut()).for_each(|((p, q), r)| {
                    add(&mut f, r, &q.neg(), p);
                });
            }
            1 => {
                // q0
                terms.iter().zip(r.iter_mut()).for_each(|((p, q), r)| {
                    add(&mut f, r, &q, p);
                });
            }
            3 => {
                // q1
                terms.iter().zip(r.iter_mut()).for_each(|((p, q), r)| {
                    add(&mut f, r, &q.to_curve().endo().to_affine().neg(), p);
                });
            }
            _ => continue,
        }
    }

    f
}

#[cfg(test)]
mod test {
    use super::super::{Fr, BW6, G1, G2};
    use super::{multi_miller_loop, Fq6, G1Affine, G2Affine, Gt};
    use ff::Field;
    use group::{prime::PrimeCurveAffine, Curve, Group};
    use pairing::{Engine as _, MillerLoopResult, PairingCurveAffine};
    use rand_core::OsRng;
    crate::test_pairing!(BW6, G1, G1Affine, G2, G2Affine, Fq6, Gt, Fr);
}

// pub fn multi_miller_loop3(terms: &[(&G1Affine, &G2Affine)], witness: &Fq6) -> Fq6 {
//     let c = witness.clone();
//     let c_inv = c.invert().unwrap();

//     let mut cq = c.clone();
//     cq.frobenius_map(1);

//     let mut cq_inv0 = cq.clone();
//     cq_inv0 = cq_inv0.invert().unwrap();

//     let mut cq_inv1 = c_inv.clone();
//     cq_inv1.frobenius_map(1);

//     assert!(cq_inv0 == cq_inv1);
//     let cq_inv = cq_inv0;

//     let mut f0 = cq_inv;

//     let mut r = terms
//         .iter()
//         .map(|(_, q)| q.to_curve().endo().neg())
//         .collect::<Vec<_>>();

//     for (_, (x2, x1)) in LOOP_2_NAF
//         .iter()
//         .rev()
//         .skip(1)
//         .zip(LOOP_1_NAF.iter().rev().skip(1))
//         .enumerate()
//     {
//         f0.square_assign();

//         let x = x2 * 3 + x1;

//         terms.iter().zip(r.iter_mut()).for_each(|((p, _), r)| {
//             double(&mut f0, r, p);
//         });

//         match x {
//             -3 => {
//                 f0 = f0 * &cq;
//                 // q1neg
//                 terms
//                     .iter()
//                     .zip(r.iter_mut())
//                     .for_each(|((p, q), r)| add(&mut f0, r, &q.to_curve().endo().to_affine(), p));
//             }
//             -1 => {
//                 f0 = f0 * &c;
//                 // q0neg
//                 terms
//                     .iter()
//                     .zip(r.iter_mut())
//                     .for_each(|((p, q), r)| add(&mut f0, r, &q.neg(), p));
//             }
//             1 => {
//                 f0 = f0 * &c_inv;
//                 // q0
//                 terms
//                     .iter()
//                     .zip(r.iter_mut())
//                     .for_each(|((p, q), r)| add(&mut f0, r, &q, p));
//             }
//             3 => {
//                 f0 = f0 * &cq_inv;
//                 // q1
//                 terms.iter().zip(r.iter_mut()).for_each(|((p, q), r)| {
//                     add(&mut f0, r, &q.to_curve().endo().to_affine().neg(), p)
//                 });
//             }
//             _ => continue,
//         }
//     }

//     f0
// }

// fn eval(f: &mut Fq6, u0: Fq, u1: Fq, u2: Fq) {
//     let sparse = Fq6 {
//         c0: Fq3 {
//             c0: u0,
//             c1: u1,
//             c2: Fq::zero(),
//         },
//         c1: Fq3 {
//             c0: Fq::zero(),
//             c1: u2,
//             c2: Fq::zero(),
//         },
//     };
//     *f = *f * &sparse;
// }

// fn add_aff(f: &mut Fq6, r: &mut G2Affine, q: &G2Affine, p: &G1Affine) {
//     {
//         assert!(r != q);
//         assert!(&r.neg() != q);
//         let rx = r.x;
//         let ry = r.y;
//         // add
//         let lmd = {
//             let qx = q.x;
//             let qy = q.y;
//             let t0 = ry - qy;
//             let t1 = rx - qx;
//             let lmd = t1.invert().unwrap() * t0;
//             let x = lmd.square() - rx - qx;
//             let y = lmd * (rx - x) - ry;
//             let u = G2Affine::from_xy(x, y).unwrap();
//             assert_eq!(r.to_curve() + q.to_curve(), u.to_curve());
//             *r = u;
//             lmd
//         };
//         // eval
//         {
//             // let u0 = lmd * rx - ry;
//             // let u1 = -(lmd * p.x);
//             // let u2 = p.y;
//             // eval(f, u0, u1, u2);

//             let nu = lmd * rx - ry;
//             let u0 = Fq::one();
//             let u1 = lmd * p.x;
//             let u2 = nu * p.y;
//             eval(f, u2, u1, u0);
//         }
//     }
// }

// #[test]
// fn test_bw666() {
//     use num_integer::Integer;
//     use rand_xorshift::XorShiftRng;

//     let mut rng = XorShiftRng::from_seed([
//         0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
//         0xe5,
//     ]);

//     let p = BigUint::from_str_radix(
//         "122e824fb83ce0ad187c94004faff3eb926186a81d14688528275ef8087be41707ba638e584e91903cebaff25b423048689c8ed12f9fd9071dcd3dc73ebff2e98a116c25667a8f8160cf8aeeaf0a437e6913e6870000082f49d00000000008b",16
//     ).unwrap();

//     let r = BigUint::from_str_radix(
//         "01ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001",16
//     ).unwrap();

//     let (h, zero) = (p.pow(6) - 1usize).div_rem(&r);
//     assert_eq!(zero, BigUint::from(0usize));

//     let fin_exp = |f: Fq6| {
//         // let s = BigUint::from(0x8508c00000000001u64 + 1);
//         let s = BigUint::from(1usize);
//         let d = &s * &h;
//         let isone = f.pow(d.to_u64_digits());
//         println!("NAIVE: {}", isone == Fq6::one());
//     };

//     {
//         {
//             let n = 2;
//             let g1 = G1Affine::generator();
//             let g2 = G2Affine::generator();
//             let scalars = (0..n)
//                 .map(|_| (Fr::random(&mut rng), Fr::random(&mut rng)))
//                 .collect::<Vec<_>>();
//             let terms = scalars
//                 .iter()
//                 .map(|(a, b)| ((g1 * a).to_affine(), (g2 * b).to_affine()))
//                 .collect::<Vec<_>>();
//             let mut terms = terms.iter().map(|(a, b)| (a, b)).collect::<Vec<_>>();

//             // let u0 = multi_miller_loop(&terms[..]);

//             let last = scalars.iter().fold(Fr::ZERO, |acc, (u0, u1)| acc + u0 * u1);
//             let negg1 = -g1;
//             let accg2 = (g2 * last).into();
//             terms.push((&negg1, &accg2));

//             let p1 = terms.iter().map(|(a, _)| *a).cloned().collect::<Vec<_>>();
//             let p2 = terms.iter().map(|(_, b)| *b).cloned().collect::<Vec<_>>();

//             {
//                 // let u0 = multi_miller_loop3(&terms[..], &Fq6::one());
//                 let (u0, u1) = multi_miller_loop_aff(&p1, &p2);

//                 println!("PC: {}", u0.final_exponentiation().0 == Fq6::one());
//                 fin_exp(u0);

//                 let t = m_inv(&u0);
//                 let is_ml = m_exp(&t);
//                 println!("mi {}", u0 == is_ml);

//                 let t = m2_inv(&u0);
//                 let is_ml = m_exp(&t);
//                 println!("m2i {}", u0 == is_ml);

//                 let t = r_inv(&u0);
//                 let is_ml = r_exp(&t);
//                 println!("ri {}", u0 == is_ml);

//                 let w = r_inv(&u0);
//                 let w = m_inv(&w);
//                 let u1 = multi_miller_loop3(&terms[..], &w);
//                 println!("XXX ??? ?   {}", u1 == Fq6::one());

//                 let w = lambda_exp(&w);
//                 let w = w * &u0.invert().unwrap();
//                 println!("ISML {}", w == Fq6::one());
//             }

//             // let t = r_inv(&t);

//             // {
//             //     println!("XXX ??? ?   {:#?}", u0);
//             //     let t = r_inv(&u0);
//             //     let is_ml = r_exp(&t);
//             //     println!("{}", u0 == is_ml);

//             //     let t = m_inv(&u0);
//             //     let is_ml = m_exp(&t);
//             //     println!("{}", u0 == is_ml);

//             //     let t = m2_inv(&u0);
//             //     let is_ml = m_exp(&t);
//             //     println!("{}", u0 == is_ml);
//             // }

//             // {
//             //     let w = r_inv(&u0);
//             //     let w = m_inv(&w);
//             //     let is_ml0 = m_exp(&r_exp(&w));
//             //     let is_ml1 = lambda_exp(&w);
//             //     println!("ISML {}", u0 == is_ml0);
//             //     println!("ISML {}", u0 == is_ml1);
//             // }

//             // // // for i in 1..=190 {
//             // // let w = r_inv(&u0);
//             // // let w = m2_inv(&w);
//             // // let wi = w.invert().unwrap();
//             // // println!("gud inv?   {}", w * wi == Fq6::one());
//             // // let u1 = multi_miller_loop3(&terms[..], &wi, None);
//             // // println!("XXX ??? ?   {}", u1 == Fq6::one());

//             // let w = r_inv(&u0);
//             // let w = m_inv(&w);

//             // let u1 = multi_miller_loop3(&terms[..], &w);
//             // println!("XXX ??? ?   {}", u1 == Fq6::one());
//             // // println!("!!! {}", slambda_exp(&w) == u1);

//             // u1.frobenius_map(1);

//             // let w = r_inv(&u0);
//             // let w = m2_inv(&w);
//             // let mut u1 = multi_miller_loop3(&terms[..], &w);
//             // u1.frobenius_map(1);
//             // println!("XXX ??? ?   {}", u1 == Fq6::one());
//             // // w.frobenius_map(i);

//             // let w = r_inv(&u0);
//             // let w = m2_inv(&w);
//             // let mut w = w.invert().unwrap();
//             // // w.frobenius_map(i);
//             // let u1 = multi_miller_loop3(&terms[..], &w, None);
//             // println!("YYY ??? ?   {}", u1 == Fq6::one());

//             // let w = r_inv(&u0);
//             // let mut w = m_inv(&w);
//             // // w.frobenius_map(i);
//             // let w = w.invert().unwrap();
//             // let u1 = multi_miller_loop3(&terms[..], &w);
//             // println!("XXX ??? ?   {}", u1 == Fq6::one());

//             // let w = r_inv(&u0);
//             // let mut w = m2_inv(&w);
//             // // w.frobenius_map(i);
//             // let w = w.invert().unwrap();
//             // let u1 = multi_miller_loop3(&terms[..], &w);
//             // println!("XXX ??? ?   {}", u1 == Fq6::one());
//             // }
//         }
//     }
// }

// fn r_exp(e: &Fq6) -> Fq6 {
//     let r = "01ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001";
//     let r = BigUint::from_str_radix(&r, 16).unwrap();
//     e.pow(r.to_u64_digits())
// }

// fn r_inv(e: &Fq6) -> Fq6 {
//     let r_i = "dd165f72375dd203252eb746be13df990bfd5452b9ce1e01a99fe6ab47de540eef3e192658eee271d9f17174ae23a0680a7e25059aa190b3527f3b4c0dd073c49d1cf0f3922a12134d592226b8bcfd6e9c25dd434d7c7dc9a59f2149eb2d279c4b5b28e70170cbe10dd765feefe5efca4ff8f9ff95efb3af0a54e0dadf422f899df03066586d51310f88c1816d0696b508ca6ad4ba2a18c412b540ed9b68e7d6fd5e9a3bee952e60bc683d04a5efd292ff351bf525467ed4afc6b75d8be1b348ede2b05edb62bbc12648d9634e82ff1c1b14dd384bca77927b701b3774a9f764bb51bdcded75c0b39babf0e68ca0c774970a5550263fd5216431e4292710e4d2d37a489dba4c4cdb36f93d523691204a72c79e1088f77f697df50f5505bd4cc620470e974730a3eaaae6fde70572bf14d98b7c99a03bab03b7a88bd88f94391c54c81b2082bd7695fb4830c7216241a950b7d5359d86fce5ab0ed93c32d8e1d39c16741db5123841a0f87fc256e1bbf887fdbe6335a18717eba7084af7b9fb7a126c5ca8ad6ed44f1505fb8162152f738717a0ad969341043ab138a16f57c287cb810a8b94bdc0631c0d7b80f970d40d5eeb745fafbd16e138afa04ac943668002e132e4e32a4106abeacdf626533f75a52eb5fe5c0fb494a17f1f644c80745fa6493d64f0230e061b5bc3dd0085b251f358c2a051d17cf6c3a571e4517617784e21ef26547136cd9426d9";
//     let r_i = BigUint::from_str_radix(&r_i, 16).unwrap();
//     e.pow(r_i.to_u64_digits())
// }

// fn m_exp(e: &Fq6) -> Fq6 {
//     let r = "184ac835122297cbd3c5cd47279693ddb3614d475afc330a22307fde5ad0866a048921a6aa874ce909e7fb6f3a17c02dd6e827cb0784603e4ab81c51bfffffc4c21a7fffffffff77";
//     let r = BigUint::from_str_radix(&r, 16).unwrap();
//     e.pow(r.to_u64_digits())
// }

// fn m_inv(e: &Fq6) -> Fq6 {
//     let r_i = "536699dbd98c57931276d1bdc2d5d4515cb10370d6c21cecdb6b351c9a1e40157d00d1b05f69e9d8cb6318edb26fe948ab63764a65f2b285a4fe4757b86ef3215b9966c92781bdbe20101ac0c2d719a16190e6c459167e22d5bc43024c279f1c888df10336c8ab28e3dc39d34027f2fb0f09ef5453c5377a1eaf8620ad036c5d821746ab754f9ba0985ecfb6ac65904ea31971bf28bf0bccde72adb0343abb1a1728fefaab909339bb0763930158bc62245870e3c7e234301cf98db672c498b16fecf0485fc5ba2d9e5dbcfc95c9f487114b4ffcb0a7485f6308a2f37aa1fda35b48d9acec42aa651830a0ccc861240e91691b35688a7fcf1ce479a7f49f330d0022323fbc3a67d90c826d462b14f1b44c674fa5ed7cc2acbf9862d3688b907a0da98668273c296c0211e907ca3966092222d8f73f0a3972b9dec1360cfe4f5a20433e4999c229ca39113118115df5ab5a77359e9d1df6c1098e4ee31ac1a45500b7ec78bfadc6dbf5fd15aa34758c8f9c85652c01e3d5c958aa9e4e0230b6bee9dd9c9a541eb7c6057fdb93950f81562ca5215539f563e1f01a19ffdc9b57f46886eef25b28f94de9e48a5c5cd8a98441da20a6b4db888960ca4779e721bd3c5ae14274b4bd12909f749781791f723b7b03b6e4d0d8ea77bf3fbc4fcc0e6b96cd6047131d7c01a1530cadda765e9ecb03009a890b6992bb724cd1aef7fc54a06c0e98cf728fc9c0884ec7";
//     let r_i = BigUint::from_str_radix(&r_i, 16).unwrap();
//     e.pow(r_i.to_u64_digits())
// }

// fn m2_inv(e: &Fq6) -> Fq6 {
//     let r_i = "22d75c0470adc82631ed3c8d6cffbb2e390c3abc3661db50a15cb267f96af927284fec35c865a11367f041fe98671c90838ec3271d62a6e4767e2e3195d6454a68e71b286bd2ecb90422e8210f0c3d7dbf996a7eb17d5c1b3d5e346f06974f30566e5e9b80ecaa3e1f9b866de056b7098871396a82abcd5c5d6be0f4ef905c0fd96211d4e5f35798929f790c22739ddc66a15e0661ef7228b7fa54ca57ff66f79741c48c6af5fde91ca5145ee4ce7415cd641201243bef0bd3143355c998b14351d1a848b0a642a8b0cbc4b5a2d629b6918a62704e9cd16cc2ec8a3997e2f4fae505cd85c35674d20ecdd2ef9c2f9a950876eac9cbd6c5d1ded1f05b5ea50a34358bb5e053d7a680e24c78baeffe3faf3a15b82eea0b0d7709b5f4fd6b7fcf57eb996c100890e559b81093e14d34a6110461734e1f5b76c72c78d73370a510260899ad6451ff975a05eaee285f19476019cffa8502a1f6cdaaae8b6f3a28bb1d36ed8489fef1823b8bbb55034eec257a7a6b4f64e181ef98fbd669b0fd5acc8688a2a24a33bef8e27ae8b6432417eef51d7ee00c385a095227c8e32a61bc20295a84e29c01151df4477e2192cab1b20e1ba68d9115801b65efb099a6efe6e92e021ac0514a28a6ed7edea65df7dbebe6f4b2883645033dd33907f5a6a33321f00e020d4f4fc487a61f98fa55a2ba757ba71864cbdc1b3a7628d3048e7fd9ff4102563f90375dab4986a995df7bacd522d5af838948c3b90bedb1675258481a397ff70883ab68c47fc8b6b8ea710de7d823b8d96b026732202657";
//     let r_i = BigUint::from_str_radix(&r_i, 16).unwrap();
//     e.pow(r_i.to_u64_digits())
// }

// fn lambda_exp(e: &Fq6) -> Fq6 {
//     let r_i = "28d323e134a5dd8d5d387d4ba8f487bf84d6e12a68cc2d52135e42d0ff3f010dd5d8b7279b08997e56ec93f17b1bc533d554658aff2a0f2bf3422352e59f8ab88cbda88b2f316fd54b4eef0b3d16c3969da45c667e271cb976b11ac1be78a3e06a74217d0720132c7499bbefffff7d906bbfffffffff77";
//     let r_i = BigUint::from_str_radix(&r_i, 16).unwrap();
//     e.pow(r_i.to_u64_digits())
// }

// fn slambda_exp(e: &Fq6) -> Fq6 {
//     let r_i = "23ed1347970dec008a442f99200000008508c00000000001";
//     let r_i = BigUint::from_str_radix(&r_i, 16).unwrap();
//     e.pow(r_i.to_u64_digits())
// }

// #[test]
// fn test_pairing_check_xxx() {
//     let n = 10;
//     let g1 = G1::generator().to_affine();
//     let g2 = G2::generator().to_affine();
//     let scalars = (0..n)
//         .map(|_| (Fr::random(OsRng), Fr::random(OsRng)))
//         .collect::<Vec<_>>();
//     let terms = scalars
//         .iter()
//         .map(|(a, b)| ((g1 * a).to_affine(), (g2 * b).to_affine()))
//         .collect::<Vec<_>>();
//     let mut terms = terms.iter().map(|(a, b)| (a, b)).collect::<Vec<_>>();
//     let gt = BW6::pairing(&g1, &g2);
//     let u0 = scalars
//         .iter()
//         .fold(Gt::identity(), |acc, (a, b)| acc + gt * a * b);
//     let u1 = multi_miller_loop(&terms[..]).final_exponentiation();
//     assert_eq!(u1, u0);

//     let last = scalars.iter().fold(Fr::ZERO, |acc, (u0, u1)| acc + u0 * u1);
//     let negg1 = -g1;
//     let accg2 = (g2 * last).into();
//     terms.push((&negg1, &accg2));
//     let mmlr = multi_miller_loop(&terms[..]);
//     let one = mmlr.final_exponentiation();
//     assert_eq!(one, Gt::identity());
// }

// pub fn multi_miller_loop(terms: &[(&G1Affine, &G2Affine)]) -> Fq6 {
//     let terms = terms
//         .iter()
//         .filter_map(|&(p, q)| {
//             (bool::from(p.is_identity()) || bool::from(q.is_identity())).then_some((*p, *q))
//         })
//         .collect::<Vec<_>>();

//     let loop_1 = (0..64)
//         .map(|i| ((0x8508c00000000001u64 + 1) >> i) & 1 == 1)
//         .rev()
//         .skip(1);
//     let mut f1 = Fq6::one();
//     let mut r = terms.iter().map(|(_, q)| q.to_curve()).collect::<Vec<_>>();

//     for (i, x) in loop_1.enumerate() {
//         (i != 0).then(|| f1.square_assign());

//         terms.iter().zip(r.iter_mut()).for_each(|((p, _), r)| {
//             double(&mut f1, r, p);
//         });

//         terms.iter().zip(r.iter_mut()).for_each(|((p, q), r)| {
//             x.then(|| add(&mut f1, r, q, p));
//         });
//     }

//     let mut f2 = Fq6::one();
//     let mut r = terms.iter().map(|(_, q)| q.to_curve()).collect::<Vec<_>>();

//     for (i, x) in LOOP_2_NAF.iter().rev().skip(1).enumerate() {
//         (i != 0).then(|| f2.square_assign());

//         terms.iter().zip(r.iter_mut()).for_each(|((p, _), r)| {
//             double(&mut f2, r, p);
//         });

//         match x {
//             &val @ (1 | -1) => {
//                 for ((p, q), r) in terms.iter().zip(r.iter_mut()) {
//                     if val == 1 {
//                         add(&mut f2, r, q, p);
//                     } else {
//                         add(&mut f2, r, &q.neg(), p);
//                     }
//                 }
//             }
//             _ => continue,
//         }
//     }

//     f2.frobenius_map(1);
//     f1 * f2
// }

// let p = BigUint::from_str_radix(
//     "122e824fb83ce0ad187c94004faff3eb926186a81d14688528275ef8087be41707ba638e584e91903cebaff25b423048689c8ed12f9fd9071dcd3dc73ebff2e98a116c25667a8f8160cf8aeeaf0a437e6913e6870000082f49d00000000008b",16
// ).unwrap();

// let r = BigUint::from_str_radix(
//     "01ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001",16
// ).unwrap();

// let (h, zero) = (p.pow(6) - 1usize).div_rem(&r);
// assert_eq!(zero, BigUint::from(0usize));
// let s = BigUint::from(0x8508c00000000001u64 + 1);
// let d = s * h;
// Gt(self.pow(d.to_u64_digits()))
// #[cfg(test)]
// mod test {
//     use super::super::{Fr, BW6, G1, G2};
//     use super::{multi_miller_loop, Fq6, G1Affine, G2Affine, Gt};
//     use ff::Field;
//     use group::{prime::PrimeCurveAffine, Curve, Group};
//     use pairing::{Engine as _, MillerLoopResult, PairingCurveAffine};
//     use rand_core::OsRng;
//     crate::test_pairing!(BW6, G1, G1Affine, G2, G2Affine, Fq6, Gt, Fr);
// }

// fn dbl_aff(f: &mut Fq6, r: &mut G2Affine, p: &G1Affine) {
//     let rx = r.x;
//     let ry = r.y;

//     let lmd = {
//         let xx = rx.square();
//         let xx3 = xx + xx + xx;
//         let yy2 = ry + ry;
//         let iyy2 = yy2.invert().unwrap();
//         let lmd = xx3 * iyy2;
//         let x = lmd.square() - rx - rx;
//         let y = lmd * (rx - x) - ry;
//         let u = G2Affine::from_xy(x, y).unwrap();
//         assert_eq!(r.to_curve().double(), u.to_curve());
//         *r = u;
//         lmd
//     };

//     // eval
//     {
//         // let u0 = lmd * rx - ry;
//         // let u1 = -(lmd * p.x);
//         // let u2 = p.y;
//         // eval(f, u0, u1, u2);

//         let nu = lmd * rx - ry;
//         let u0 = Fq::one();
//         let u1 = lmd * p.x;
//         let u2 = nu * p.y;
//         eval(f, u2, u1, u0);
//     }
// }

// pub fn multi_miller_loop_aff(p1: &[G1Affine], p2: &[G2Affine]) -> (Fq6, Fq6) {
//     // Follows 6th equation at https://hackmd.io/@gnark/BW6-761-changes
//     let mut r_aff = p2
//         .iter()
//         .map(|q| q.to_curve().endo().neg().to_affine())
//         .collect::<Vec<_>>();

//     let mut r_pro = p2
//         .iter()
//         .map(|q| q.to_curve().endo().neg())
//         .collect::<Vec<_>>();

//     let p1 = p1
//         .iter()
//         .map(|p| {
//             let y = p.y.invert().unwrap();
//             let x = (p.x * y).neg();
//             G1Affine { x, y }
//         })
//         .collect::<Vec<_>>();

//     let f_aff = &mut Fq6::one();
//     let f_pro = &mut Fq6::one();

//     let p2_endo = p2
//         .iter()
//         .map(|q| q.to_curve().endo().to_affine())
//         .collect::<Vec<_>>();

//     let mut c = 0;
//     for (x2, x1) in LOOP_2_NAF
//         .iter()
//         .skip(1)
//         .rev()
//         .skip(1)
//         .zip(LOOP_1_NAF.iter().skip(1).rev().skip(1))
//     {
//         let x = x2 * 3 + x1;
//         f_aff.square_assign();
//         f_pro.square_assign();

//         // println!("CAFF {} {}", c, x);
//         p1.iter()
//             .zip(r_aff.iter_mut())
//             .zip(r_pro.iter_mut())
//             .for_each(|((p, r_aff), r_pro)| {
//                 dbl_aff(f_aff, r_aff, p);
//                 double(f_pro, r_pro, p);
//             });
//         c += 1;

//         match x {
//             -3 => {
//                 p1.iter()
//                     .zip(p2_endo.iter())
//                     .zip(r_aff.iter_mut())
//                     .zip(r_pro.iter_mut())
//                     .for_each(|(((p1, p2), r_aff), r_pro)| {
//                         add_aff(f_aff, r_aff, &p2, p1);
//                         add(f_pro, r_pro, &p2, p1);
//                     });
//             }
//             -1 => {
//                 p1.iter()
//                     .zip(p2)
//                     .zip(r_aff.iter_mut())
//                     .zip(r_pro.iter_mut())
//                     .for_each(|(((p1, p2), r_aff), r_pro)| {
//                         add_aff(f_aff, r_aff, &p2.neg(), p1);
//                         add(f_pro, r_pro, &p2.neg(), p1);
//                     });
//             }
//             1 => {
//                 p1.iter()
//                     .zip(p2)
//                     .zip(r_aff.iter_mut())
//                     .zip(r_pro.iter_mut())
//                     .for_each(|(((p1, p2), r_aff), r_pro)| {
//                         let qx = *p2;
//                         add_aff(f_aff, r_aff, &qx, p1);
//                         add(f_pro, r_pro, &qx, p1);
//                     });
//             }
//             3 => {
//                 p1.iter()
//                     .zip(p2_endo.iter())
//                     .zip(r_aff.iter_mut())
//                     .zip(r_pro.iter_mut())
//                     .for_each(|(((p1, p2), r_aff), r_pro)| {
//                         add_aff(f_aff, r_aff, &p2.neg(), p1);
//                         add(f_pro, r_pro, &p2.neg(), p1);
//                     });
//             }
//             _ => continue,
//         }
//     }
//     f_aff.square_assign();
//     f_pro.square_assign();

//     p1.iter()
//         .zip(r_aff.iter_mut())
//         .zip(r_pro.iter_mut())
//         .for_each(|((p, r_aff), r_pro)| {
//             dbl_aff(f_aff, r_aff, p);
//             double(f_pro, r_pro, p);
//         });

//     p1.iter()
//         .zip(p2_endo.iter())
//         .zip(r_aff.iter_mut())
//         .zip(r_pro.iter_mut())
//         .for_each(|(((p1, p2), r_aff), r_pro)| {
//             add(f_pro, r_pro, &p2, p1);
//         });

//     (*f_aff, *f_pro)
// }
