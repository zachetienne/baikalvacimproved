
#ifndef __FD_FUNCTIONS_H__
#define __FD_FUNCTIONS_H__
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#define _UNUSED   CCTK_ATTRIBUTE_UNUSED
#define _NOINLINE CCTK_ATTRIBUTE_NOINLINE
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dKOD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_16 = 1.0/16.0;
  const REAL_SIMD_ARRAY _Rational_1_16 = ConstSIMD(tmp_Rational_1_16);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_8 = 3.0/8.0;
  const REAL_SIMD_ARRAY _Rational_3_8 = ConstSIMD(tmp_Rational_3_8);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_1_4, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), NegFusedMulSubSIMD(_Rational_1_16, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2), MulSIMD(_Rational_3_8, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dKOD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_16 = 1.0/16.0;
  const REAL_SIMD_ARRAY _Rational_1_16 = ConstSIMD(tmp_Rational_1_16);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_8 = 3.0/8.0;
  const REAL_SIMD_ARRAY _Rational_3_8 = ConstSIMD(tmp_Rational_3_8);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_4, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), NegFusedMulSubSIMD(_Rational_1_16, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2), MulSIMD(_Rational_3_8, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dKOD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_16 = 1.0/16.0;
  const REAL_SIMD_ARRAY _Rational_1_16 = ConstSIMD(tmp_Rational_1_16);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_8 = 3.0/8.0;
  const REAL_SIMD_ARRAY _Rational_3_8 = ConstSIMD(tmp_Rational_3_8);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_4, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), NegFusedMulSubSIMD(_Rational_1_16, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2), MulSIMD(_Rational_3_8, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_ddnD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_3_2, f, FusedMulSubSIMD(_Rational_1_2, f_i0m2_i1_i2, MulSIMD(_Integer_2, f_i0m1_i1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_ddnD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_3_2, f, FusedMulSubSIMD(_Rational_1_2, f_i0_i1m2_i2, MulSIMD(_Integer_2, f_i0_i1m1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_ddnD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_3_2, f, FusedMulSubSIMD(_Rational_1_2, f_i0_i1_i2m2, MulSIMD(_Integer_2, f_i0_i1_i2m1))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dupD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  return MulSIMD(invdx0, SubSIMD(FusedMulSubSIMD(_Integer_2, f_i0p1_i1_i2, MulSIMD(_Rational_1_2, f_i0p2_i1_i2)), MulSIMD(_Rational_3_2, f)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dupD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  return MulSIMD(invdx1, SubSIMD(FusedMulSubSIMD(_Integer_2, f_i0_i1p1_i2, MulSIMD(_Rational_1_2, f_i0_i1p2_i2)), MulSIMD(_Rational_3_2, f)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dupD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  return MulSIMD(invdx2, SubSIMD(FusedMulSubSIMD(_Integer_2, f_i0_i1_i2p1, MulSIMD(_Rational_1_2, f_i0_i1_i2p2)), MulSIMD(_Rational_3_2, f)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f_i0p1_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  return MulSIMD(_Rational_1_2, MulSIMD(invdx0, SubSIMD(f_i0p1_i1_i2, f_i0m1_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f_i0_i1p1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  return MulSIMD(_Rational_1_2, MulSIMD(invdx1, SubSIMD(f_i0_i1p1_i2, f_i0_i1m1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f_i0_i1_i2p1) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  return MulSIMD(_Rational_1_2, MulSIMD(invdx2, SubSIMD(f_i0_i1_i2p1, f_i0_i1_i2m1)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dDD00(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  return MulSIMD(MulSIMD(invdx0, invdx0), AddSIMD(f_i0p1_i1_i2, NegFusedMulAddSIMD(_Integer_2, f, f_i0m1_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dDD01(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0m1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p1_i1m1_i2,const REAL_SIMD_ARRAY f_i0m1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p1_i1p1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  return MulSIMD(MulSIMD(_Rational_1_4, invdx0), MulSIMD(invdx1, AddSIMD(f_i0p1_i1p1_i2, SubSIMD(SubSIMD(f_i0m1_i1m1_i2, f_i0m1_i1p1_i2), f_i0p1_i1m1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dDD02(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0m1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p1_i1_i2m1,const REAL_SIMD_ARRAY f_i0m1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p1_i1_i2p1) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  return MulSIMD(MulSIMD(_Rational_1_4, invdx0), MulSIMD(invdx2, AddSIMD(f_i0p1_i1_i2p1, SubSIMD(SubSIMD(f_i0m1_i1_i2m1, f_i0m1_i1_i2p1), f_i0p1_i1_i2m1))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dDD11(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  return MulSIMD(MulSIMD(invdx1, invdx1), AddSIMD(f_i0_i1p1_i2, NegFusedMulAddSIMD(_Integer_2, f, f_i0_i1m1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dDD12(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1m1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p1_i2m1,const REAL_SIMD_ARRAY f_i0_i1m1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p1_i2p1) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  return MulSIMD(MulSIMD(_Rational_1_4, invdx1), MulSIMD(invdx2, AddSIMD(f_i0_i1p1_i2p1, SubSIMD(SubSIMD(f_i0_i1m1_i2m1, f_i0_i1m1_i2p1), f_i0_i1p1_i2m1))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_2_f_dDD22(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1) {

  const double tmp_Integer_2 = 2.0;
  const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  return MulSIMD(MulSIMD(invdx2, invdx2), AddSIMD(f_i0_i1_i2p1, NegFusedMulAddSIMD(_Integer_2, f, f_i0_i1_i2m1)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dKOD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_15_64 = 15.0/64.0;
  const REAL_SIMD_ARRAY _Rational_15_64 = ConstSIMD(tmp_Rational_15_64);

  const double tmp_Rational_1_64 = 1.0/64.0;
  const REAL_SIMD_ARRAY _Rational_1_64 = ConstSIMD(tmp_Rational_1_64);

  const double tmp_Rational_3_32 = 3.0/32.0;
  const REAL_SIMD_ARRAY _Rational_3_32 = ConstSIMD(tmp_Rational_3_32);

  const double tmp_Rational_5_16 = 5.0/16.0;
  const REAL_SIMD_ARRAY _Rational_5_16 = ConstSIMD(tmp_Rational_5_16);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_1_64, AddSIMD(f_i0m3_i1_i2, f_i0p3_i1_i2), FusedMulAddSIMD(_Rational_3_32, MulSIMD(_NegativeOne_, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2)), FusedMulSubSIMD(_Rational_15_64, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), MulSIMD(_Rational_5_16, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dKOD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_15_64 = 15.0/64.0;
  const REAL_SIMD_ARRAY _Rational_15_64 = ConstSIMD(tmp_Rational_15_64);

  const double tmp_Rational_1_64 = 1.0/64.0;
  const REAL_SIMD_ARRAY _Rational_1_64 = ConstSIMD(tmp_Rational_1_64);

  const double tmp_Rational_3_32 = 3.0/32.0;
  const REAL_SIMD_ARRAY _Rational_3_32 = ConstSIMD(tmp_Rational_3_32);

  const double tmp_Rational_5_16 = 5.0/16.0;
  const REAL_SIMD_ARRAY _Rational_5_16 = ConstSIMD(tmp_Rational_5_16);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_64, AddSIMD(f_i0_i1m3_i2, f_i0_i1p3_i2), FusedMulAddSIMD(_Rational_3_32, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2)), FusedMulSubSIMD(_Rational_15_64, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), MulSIMD(_Rational_5_16, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dKOD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_15_64 = 15.0/64.0;
  const REAL_SIMD_ARRAY _Rational_15_64 = ConstSIMD(tmp_Rational_15_64);

  const double tmp_Rational_1_64 = 1.0/64.0;
  const REAL_SIMD_ARRAY _Rational_1_64 = ConstSIMD(tmp_Rational_1_64);

  const double tmp_Rational_3_32 = 3.0/32.0;
  const REAL_SIMD_ARRAY _Rational_3_32 = ConstSIMD(tmp_Rational_3_32);

  const double tmp_Rational_5_16 = 5.0/16.0;
  const REAL_SIMD_ARRAY _Rational_5_16 = ConstSIMD(tmp_Rational_5_16);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_64, AddSIMD(f_i0_i1_i2m3, f_i0_i1_i2p3), FusedMulAddSIMD(_Rational_3_32, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2)), FusedMulSubSIMD(_Rational_15_64, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), MulSIMD(_Rational_5_16, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_ddnD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_5_6 = 5.0/6.0;
  const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_1_4, f_i0p1_i1_i2, FusedMulAddSIMD(_Rational_5_6, f, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0m2_i1_i2, MulSIMD(_Rational_1_12, f_i0m3_i1_i2)), MulSIMD(_Rational_3_2, f_i0m1_i1_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_ddnD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_5_6 = 5.0/6.0;
  const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_4, f_i0_i1p1_i2, FusedMulAddSIMD(_Rational_5_6, f, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1m2_i2, MulSIMD(_Rational_1_12, f_i0_i1m3_i2)), MulSIMD(_Rational_3_2, f_i0_i1m1_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_ddnD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_5_6 = 5.0/6.0;
  const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_4, f_i0_i1_i2p1, FusedMulAddSIMD(_Rational_5_6, f, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1_i2m2, MulSIMD(_Rational_1_12, f_i0_i1_i2m3)), MulSIMD(_Rational_3_2, f_i0_i1_i2m1)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dupD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_5_6 = 5.0/6.0;
  const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_3_2, f_i0p1_i1_i2, SubSIMD(SubSIMD(FusedMulSubSIMD(_Rational_1_12, f_i0p3_i1_i2, MulSIMD(_Rational_1_2, f_i0p2_i1_i2)), MulSIMD(_Rational_5_6, f)), MulSIMD(_Rational_1_4, f_i0m1_i1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dupD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_5_6 = 5.0/6.0;
  const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_3_2, f_i0_i1p1_i2, SubSIMD(SubSIMD(FusedMulSubSIMD(_Rational_1_12, f_i0_i1p3_i2, MulSIMD(_Rational_1_2, f_i0_i1p2_i2)), MulSIMD(_Rational_5_6, f)), MulSIMD(_Rational_1_4, f_i0_i1m1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dupD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_4 = 1.0/4.0;
  const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_5_6 = 5.0/6.0;
  const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_3_2, f_i0_i1_i2p1, SubSIMD(SubSIMD(FusedMulSubSIMD(_Rational_1_12, f_i0_i1_i2p3, MulSIMD(_Rational_1_2, f_i0_i1_i2p2)), MulSIMD(_Rational_5_6, f)), MulSIMD(_Rational_1_4, f_i0_i1_i2m1))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_2_3 = 2.0/3.0;
  const REAL_SIMD_ARRAY _Rational_2_3 = ConstSIMD(tmp_Rational_2_3);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_1_12, SubSIMD(f_i0m2_i1_i2, f_i0p2_i1_i2), MulSIMD(_Rational_2_3, SubSIMD(f_i0p1_i1_i2, f_i0m1_i1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_2_3 = 2.0/3.0;
  const REAL_SIMD_ARRAY _Rational_2_3 = ConstSIMD(tmp_Rational_2_3);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_12, SubSIMD(f_i0_i1m2_i2, f_i0_i1p2_i2), MulSIMD(_Rational_2_3, SubSIMD(f_i0_i1p1_i2, f_i0_i1m1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_2_3 = 2.0/3.0;
  const REAL_SIMD_ARRAY _Rational_2_3 = ConstSIMD(tmp_Rational_2_3);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_12, SubSIMD(f_i0_i1_i2m2, f_i0_i1_i2p2), MulSIMD(_Rational_2_3, SubSIMD(f_i0_i1_i2p1, f_i0_i1_i2m1))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD00(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_5_2 = 5.0/2.0;
  const REAL_SIMD_ARRAY _Rational_5_2 = ConstSIMD(tmp_Rational_5_2);

  return MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(_Rational_4_3, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), NegFusedMulSubSIMD(_Rational_1_12, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2), MulSIMD(_Rational_5_2, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD01(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0m2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p2_i2,const REAL_SIMD_ARRAY f_i0m1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p2_i1p2_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_144 = 1.0/144.0;
  const REAL_SIMD_ARRAY _Rational_1_144 = ConstSIMD(tmp_Rational_1_144);

  const double tmp_Rational_1_18 = 1.0/18.0;
  const REAL_SIMD_ARRAY _Rational_1_18 = ConstSIMD(tmp_Rational_1_18);

  const double tmp_Rational_4_9 = 4.0/9.0;
  const REAL_SIMD_ARRAY _Rational_4_9 = ConstSIMD(tmp_Rational_4_9);

  return MulSIMD(invdx0, MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_18, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1m1_i2, f_i0m2_i1m1_i2), AddSIMD(AddSIMD(f_i0m2_i1p1_i2, f_i0p1_i1m2_i2), SubSIMD(SubSIMD(f_i0m1_i1p2_i2, f_i0m1_i1m2_i2), f_i0p2_i1p1_i2))), f_i0p1_i1p2_i2), FusedMulAddSIMD(_Rational_4_9, AddSIMD(f_i0p1_i1p1_i2, SubSIMD(SubSIMD(f_i0m1_i1m1_i2, f_i0m1_i1p1_i2), f_i0p1_i1m1_i2)), MulSIMD(_Rational_1_144, AddSIMD(f_i0p2_i1p2_i2, SubSIMD(SubSIMD(f_i0m2_i1m2_i2, f_i0m2_i1p2_i2), f_i0p2_i1m2_i2)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD02(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p2,const REAL_SIMD_ARRAY f_i0m1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p2_i1_i2p2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_144 = 1.0/144.0;
  const REAL_SIMD_ARRAY _Rational_1_144 = ConstSIMD(tmp_Rational_1_144);

  const double tmp_Rational_1_18 = 1.0/18.0;
  const REAL_SIMD_ARRAY _Rational_1_18 = ConstSIMD(tmp_Rational_1_18);

  const double tmp_Rational_4_9 = 4.0/9.0;
  const REAL_SIMD_ARRAY _Rational_4_9 = ConstSIMD(tmp_Rational_4_9);

  return MulSIMD(invdx0, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_18, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1_i2m1, f_i0m2_i1_i2m1), AddSIMD(AddSIMD(f_i0m2_i1_i2p1, f_i0p1_i1_i2m2), SubSIMD(SubSIMD(f_i0m1_i1_i2p2, f_i0m1_i1_i2m2), f_i0p2_i1_i2p1))), f_i0p1_i1_i2p2), FusedMulAddSIMD(_Rational_4_9, AddSIMD(f_i0p1_i1_i2p1, SubSIMD(SubSIMD(f_i0m1_i1_i2m1, f_i0m1_i1_i2p1), f_i0p1_i1_i2m1)), MulSIMD(_Rational_1_144, AddSIMD(f_i0p2_i1_i2p2, SubSIMD(SubSIMD(f_i0m2_i1_i2m2, f_i0m2_i1_i2p2), f_i0p2_i1_i2m2)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD11(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_5_2 = 5.0/2.0;
  const REAL_SIMD_ARRAY _Rational_5_2 = ConstSIMD(tmp_Rational_5_2);

  return MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(_Rational_4_3, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), NegFusedMulSubSIMD(_Rational_1_12, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2), MulSIMD(_Rational_5_2, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD12(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p2,const REAL_SIMD_ARRAY f_i0_i1m1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p2_i2p2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_144 = 1.0/144.0;
  const REAL_SIMD_ARRAY _Rational_1_144 = ConstSIMD(tmp_Rational_1_144);

  const double tmp_Rational_1_18 = 1.0/18.0;
  const REAL_SIMD_ARRAY _Rational_1_18 = ConstSIMD(tmp_Rational_1_18);

  const double tmp_Rational_4_9 = 4.0/9.0;
  const REAL_SIMD_ARRAY _Rational_4_9 = ConstSIMD(tmp_Rational_4_9);

  return MulSIMD(invdx1, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_18, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p2_i2m1, f_i0_i1m2_i2m1), AddSIMD(AddSIMD(f_i0_i1m2_i2p1, f_i0_i1p1_i2m2), SubSIMD(SubSIMD(f_i0_i1m1_i2p2, f_i0_i1m1_i2m2), f_i0_i1p2_i2p1))), f_i0_i1p1_i2p2), FusedMulAddSIMD(_Rational_4_9, AddSIMD(f_i0_i1p1_i2p1, SubSIMD(SubSIMD(f_i0_i1m1_i2m1, f_i0_i1m1_i2p1), f_i0_i1p1_i2m1)), MulSIMD(_Rational_1_144, AddSIMD(f_i0_i1p2_i2p2, SubSIMD(SubSIMD(f_i0_i1m2_i2m2, f_i0_i1m2_i2p2), f_i0_i1p2_i2m2)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD22(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_12 = 1.0/12.0;
  const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_5_2 = 5.0/2.0;
  const REAL_SIMD_ARRAY _Rational_5_2 = ConstSIMD(tmp_Rational_5_2);

  return MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(_Rational_4_3, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), NegFusedMulSubSIMD(_Rational_1_12, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2), MulSIMD(_Rational_5_2, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dKOD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_256 = 1.0/256.0;
  const REAL_SIMD_ARRAY _Rational_1_256 = ConstSIMD(tmp_Rational_1_256);

  const double tmp_Rational_1_32 = 1.0/32.0;
  const REAL_SIMD_ARRAY _Rational_1_32 = ConstSIMD(tmp_Rational_1_32);

  const double tmp_Rational_35_128 = 35.0/128.0;
  const REAL_SIMD_ARRAY _Rational_35_128 = ConstSIMD(tmp_Rational_35_128);

  const double tmp_Rational_7_32 = 7.0/32.0;
  const REAL_SIMD_ARRAY _Rational_7_32 = ConstSIMD(tmp_Rational_7_32);

  const double tmp_Rational_7_64 = 7.0/64.0;
  const REAL_SIMD_ARRAY _Rational_7_64 = ConstSIMD(tmp_Rational_7_64);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_1_32, AddSIMD(f_i0m3_i1_i2, f_i0p3_i1_i2), FusedMulAddSIMD(_Rational_7_32, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), FusedMulAddSIMD(_Rational_7_64, MulSIMD(_NegativeOne_, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2)), NegFusedMulSubSIMD(_Rational_1_256, AddSIMD(f_i0p4_i1_i2, f_i0m4_i1_i2), MulSIMD(_Rational_35_128, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dKOD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_256 = 1.0/256.0;
  const REAL_SIMD_ARRAY _Rational_1_256 = ConstSIMD(tmp_Rational_1_256);

  const double tmp_Rational_1_32 = 1.0/32.0;
  const REAL_SIMD_ARRAY _Rational_1_32 = ConstSIMD(tmp_Rational_1_32);

  const double tmp_Rational_35_128 = 35.0/128.0;
  const REAL_SIMD_ARRAY _Rational_35_128 = ConstSIMD(tmp_Rational_35_128);

  const double tmp_Rational_7_32 = 7.0/32.0;
  const REAL_SIMD_ARRAY _Rational_7_32 = ConstSIMD(tmp_Rational_7_32);

  const double tmp_Rational_7_64 = 7.0/64.0;
  const REAL_SIMD_ARRAY _Rational_7_64 = ConstSIMD(tmp_Rational_7_64);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_32, AddSIMD(f_i0_i1m3_i2, f_i0_i1p3_i2), FusedMulAddSIMD(_Rational_7_32, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), FusedMulAddSIMD(_Rational_7_64, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2)), NegFusedMulSubSIMD(_Rational_1_256, AddSIMD(f_i0_i1p4_i2, f_i0_i1m4_i2), MulSIMD(_Rational_35_128, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dKOD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_256 = 1.0/256.0;
  const REAL_SIMD_ARRAY _Rational_1_256 = ConstSIMD(tmp_Rational_1_256);

  const double tmp_Rational_1_32 = 1.0/32.0;
  const REAL_SIMD_ARRAY _Rational_1_32 = ConstSIMD(tmp_Rational_1_32);

  const double tmp_Rational_35_128 = 35.0/128.0;
  const REAL_SIMD_ARRAY _Rational_35_128 = ConstSIMD(tmp_Rational_35_128);

  const double tmp_Rational_7_32 = 7.0/32.0;
  const REAL_SIMD_ARRAY _Rational_7_32 = ConstSIMD(tmp_Rational_7_32);

  const double tmp_Rational_7_64 = 7.0/64.0;
  const REAL_SIMD_ARRAY _Rational_7_64 = ConstSIMD(tmp_Rational_7_64);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_32, AddSIMD(f_i0_i1_i2m3, f_i0_i1_i2p3), FusedMulAddSIMD(_Rational_7_32, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), FusedMulAddSIMD(_Rational_7_64, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2)), NegFusedMulSubSIMD(_Rational_1_256, AddSIMD(f_i0_i1_i2p4, f_i0_i1_i2m4), MulSIMD(_Rational_35_128, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_ddnD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_30 = 1.0/30.0;
  const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_2_15 = 2.0/15.0;
  const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

  const double tmp_Rational_2_5 = 2.0/5.0;
  const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_7_12 = 7.0/12.0;
  const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_7_12, f, SubSIMD(FusedMulAddSIMD(_Rational_1_60, f_i0m4_i1_i2, FusedMulAddSIMD(_Rational_2_5, f_i0p1_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0m2_i1_i2, MulSIMD(_Rational_1_30, f_i0p2_i1_i2)), MulSIMD(_Rational_4_3, f_i0m1_i1_i2)))), MulSIMD(_Rational_2_15, f_i0m3_i1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_ddnD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_30 = 1.0/30.0;
  const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_2_15 = 2.0/15.0;
  const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

  const double tmp_Rational_2_5 = 2.0/5.0;
  const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_7_12 = 7.0/12.0;
  const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_7_12, f, SubSIMD(FusedMulAddSIMD(_Rational_1_60, f_i0_i1m4_i2, FusedMulAddSIMD(_Rational_2_5, f_i0_i1p1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1m2_i2, MulSIMD(_Rational_1_30, f_i0_i1p2_i2)), MulSIMD(_Rational_4_3, f_i0_i1m1_i2)))), MulSIMD(_Rational_2_15, f_i0_i1m3_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_ddnD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_30 = 1.0/30.0;
  const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_2_15 = 2.0/15.0;
  const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

  const double tmp_Rational_2_5 = 2.0/5.0;
  const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_7_12 = 7.0/12.0;
  const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_7_12, f, SubSIMD(FusedMulAddSIMD(_Rational_1_60, f_i0_i1_i2m4, FusedMulAddSIMD(_Rational_2_5, f_i0_i1_i2p1, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1_i2m2, MulSIMD(_Rational_1_30, f_i0_i1_i2p2)), MulSIMD(_Rational_4_3, f_i0_i1_i2m1)))), MulSIMD(_Rational_2_15, f_i0_i1_i2m3))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dupD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_30 = 1.0/30.0;
  const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_2_15 = 2.0/15.0;
  const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

  const double tmp_Rational_2_5 = 2.0/5.0;
  const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_7_12 = 7.0/12.0;
  const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

  return MulSIMD(invdx0, SubSIMD(SubSIMD(FusedMulAddSIMD(_Rational_2_15, f_i0p3_i1_i2, FusedMulAddSIMD(_Rational_4_3, f_i0p1_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_30, f_i0m2_i1_i2, MulSIMD(_Rational_1_2, f_i0p2_i1_i2)), MulSIMD(_Rational_7_12, f)))), MulSIMD(_Rational_2_5, f_i0m1_i1_i2)), MulSIMD(_Rational_1_60, f_i0p4_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dupD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_30 = 1.0/30.0;
  const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_2_15 = 2.0/15.0;
  const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

  const double tmp_Rational_2_5 = 2.0/5.0;
  const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_7_12 = 7.0/12.0;
  const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

  return MulSIMD(invdx1, SubSIMD(SubSIMD(FusedMulAddSIMD(_Rational_2_15, f_i0_i1p3_i2, FusedMulAddSIMD(_Rational_4_3, f_i0_i1p1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_30, f_i0_i1m2_i2, MulSIMD(_Rational_1_2, f_i0_i1p2_i2)), MulSIMD(_Rational_7_12, f)))), MulSIMD(_Rational_2_5, f_i0_i1m1_i2)), MulSIMD(_Rational_1_60, f_i0_i1p4_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dupD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_30 = 1.0/30.0;
  const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_2_15 = 2.0/15.0;
  const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

  const double tmp_Rational_2_5 = 2.0/5.0;
  const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

  const double tmp_Rational_4_3 = 4.0/3.0;
  const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

  const double tmp_Rational_7_12 = 7.0/12.0;
  const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

  return MulSIMD(invdx2, SubSIMD(SubSIMD(FusedMulAddSIMD(_Rational_2_15, f_i0_i1_i2p3, FusedMulAddSIMD(_Rational_4_3, f_i0_i1_i2p1, SubSIMD(FusedMulSubSIMD(_Rational_1_30, f_i0_i1_i2m2, MulSIMD(_Rational_1_2, f_i0_i1_i2p2)), MulSIMD(_Rational_7_12, f)))), MulSIMD(_Rational_2_5, f_i0_i1_i2m1)), MulSIMD(_Rational_1_60, f_i0_i1_i2p4)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_3_20 = 3.0/20.0;
  const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

  const double tmp_Rational_3_4 = 3.0/4.0;
  const REAL_SIMD_ARRAY _Rational_3_4 = ConstSIMD(tmp_Rational_3_4);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_3_20, SubSIMD(f_i0m2_i1_i2, f_i0p2_i1_i2), FusedMulAddSIMD(_Rational_3_4, SubSIMD(f_i0p1_i1_i2, f_i0m1_i1_i2), MulSIMD(_Rational_1_60, SubSIMD(f_i0p3_i1_i2, f_i0m3_i1_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_3_20 = 3.0/20.0;
  const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

  const double tmp_Rational_3_4 = 3.0/4.0;
  const REAL_SIMD_ARRAY _Rational_3_4 = ConstSIMD(tmp_Rational_3_4);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_3_20, SubSIMD(f_i0_i1m2_i2, f_i0_i1p2_i2), FusedMulAddSIMD(_Rational_3_4, SubSIMD(f_i0_i1p1_i2, f_i0_i1m1_i2), MulSIMD(_Rational_1_60, SubSIMD(f_i0_i1p3_i2, f_i0_i1m3_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_60 = 1.0/60.0;
  const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

  const double tmp_Rational_3_20 = 3.0/20.0;
  const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

  const double tmp_Rational_3_4 = 3.0/4.0;
  const REAL_SIMD_ARRAY _Rational_3_4 = ConstSIMD(tmp_Rational_3_4);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_3_20, SubSIMD(f_i0_i1_i2m2, f_i0_i1_i2p2), FusedMulAddSIMD(_Rational_3_4, SubSIMD(f_i0_i1_i2p1, f_i0_i1_i2m1), MulSIMD(_Rational_1_60, SubSIMD(f_i0_i1_i2p3, f_i0_i1_i2m3)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD00(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_90 = 1.0/90.0;
  const REAL_SIMD_ARRAY _Rational_1_90 = ConstSIMD(tmp_Rational_1_90);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_3_20 = 3.0/20.0;
  const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

  const double tmp_Rational_49_18 = 49.0/18.0;
  const REAL_SIMD_ARRAY _Rational_49_18 = ConstSIMD(tmp_Rational_49_18);

  return MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(_Rational_3_2, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), FusedMulAddSIMD(_Rational_3_20, MulSIMD(_NegativeOne_, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2)), FusedMulSubSIMD(_Rational_1_90, AddSIMD(f_i0m3_i1_i2, f_i0p3_i1_i2), MulSIMD(_Rational_49_18, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD01(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0m3_i1m3_i2,const REAL_SIMD_ARRAY f_i0m2_i1m3_i2,const REAL_SIMD_ARRAY f_i0m1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p2_i1m3_i2,const REAL_SIMD_ARRAY f_i0p3_i1m3_i2,const REAL_SIMD_ARRAY f_i0m3_i1m2_i2,const REAL_SIMD_ARRAY f_i0m2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p2_i1m2_i2,const REAL_SIMD_ARRAY f_i0p3_i1m2_i2,const REAL_SIMD_ARRAY f_i0m3_i1m1_i2,const REAL_SIMD_ARRAY f_i0m2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p2_i1m1_i2,const REAL_SIMD_ARRAY f_i0p3_i1m1_i2,const REAL_SIMD_ARRAY f_i0m3_i1p1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p2_i1p1_i2,const REAL_SIMD_ARRAY f_i0p3_i1p1_i2,const REAL_SIMD_ARRAY f_i0m3_i1p2_i2,const REAL_SIMD_ARRAY f_i0m2_i1p2_i2,const REAL_SIMD_ARRAY f_i0m1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p2_i1p2_i2,const REAL_SIMD_ARRAY f_i0p3_i1p2_i2,const REAL_SIMD_ARRAY f_i0m3_i1p3_i2,const REAL_SIMD_ARRAY f_i0m2_i1p3_i2,const REAL_SIMD_ARRAY f_i0m1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p2_i1p3_i2,const REAL_SIMD_ARRAY f_i0p3_i1p3_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_3600 = 1.0/3600.0;
  const REAL_SIMD_ARRAY _Rational_1_3600 = ConstSIMD(tmp_Rational_1_3600);

  const double tmp_Rational_1_400 = 1.0/400.0;
  const REAL_SIMD_ARRAY _Rational_1_400 = ConstSIMD(tmp_Rational_1_400);

  const double tmp_Rational_1_80 = 1.0/80.0;
  const REAL_SIMD_ARRAY _Rational_1_80 = ConstSIMD(tmp_Rational_1_80);

  const double tmp_Rational_9_16 = 9.0/16.0;
  const REAL_SIMD_ARRAY _Rational_9_16 = ConstSIMD(tmp_Rational_9_16);

  const double tmp_Rational_9_400 = 9.0/400.0;
  const REAL_SIMD_ARRAY _Rational_9_400 = ConstSIMD(tmp_Rational_9_400);

  const double tmp_Rational_9_80 = 9.0/80.0;
  const REAL_SIMD_ARRAY _Rational_9_80 = ConstSIMD(tmp_Rational_9_80);

  return MulSIMD(invdx0, MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_80, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1p1_i2, f_i0m3_i1p1_i2), AddSIMD(AddSIMD(f_i0m3_i1m1_i2, f_i0p1_i1p3_i2), SubSIMD(SubSIMD(f_i0m1_i1m3_i2, f_i0m1_i1p3_i2), f_i0p3_i1m1_i2))), f_i0p1_i1m3_i2), FusedMulAddSIMD(_Rational_9_16, AddSIMD(f_i0p1_i1p1_i2, SubSIMD(SubSIMD(f_i0m1_i1m1_i2, f_i0m1_i1p1_i2), f_i0p1_i1m1_i2)), FusedMulAddSIMD(_Rational_9_400, AddSIMD(f_i0p2_i1p2_i2, SubSIMD(SubSIMD(f_i0m2_i1m2_i2, f_i0m2_i1p2_i2), f_i0p2_i1m2_i2)), FusedMulAddSIMD(_Rational_9_80, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1m1_i2, f_i0m2_i1m1_i2), AddSIMD(AddSIMD(f_i0m2_i1p1_i2, f_i0p1_i1m2_i2), SubSIMD(SubSIMD(f_i0m1_i1p2_i2, f_i0m1_i1m2_i2), f_i0p2_i1p1_i2))), f_i0p1_i1p2_i2), FusedMulAddSIMD(_Rational_1_3600, AddSIMD(f_i0p3_i1p3_i2, SubSIMD(SubSIMD(f_i0m3_i1m3_i2, f_i0m3_i1p3_i2), f_i0p3_i1m3_i2)), MulSIMD(_Rational_1_400, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1m2_i2, f_i0m3_i1m2_i2), AddSIMD(AddSIMD(f_i0m3_i1p2_i2, f_i0p2_i1m3_i2), SubSIMD(SubSIMD(f_i0m2_i1p3_i2, f_i0m2_i1m3_i2), f_i0p3_i1p2_i2))), f_i0p2_i1p3_i2)))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD02(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0m3_i1_i2m3,const REAL_SIMD_ARRAY f_i0m2_i1_i2m3,const REAL_SIMD_ARRAY f_i0m1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p2_i1_i2m3,const REAL_SIMD_ARRAY f_i0p3_i1_i2m3,const REAL_SIMD_ARRAY f_i0m3_i1_i2m2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p2_i1_i2m2,const REAL_SIMD_ARRAY f_i0p3_i1_i2m2,const REAL_SIMD_ARRAY f_i0m3_i1_i2m1,const REAL_SIMD_ARRAY f_i0m2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p2_i1_i2m1,const REAL_SIMD_ARRAY f_i0p3_i1_i2m1,const REAL_SIMD_ARRAY f_i0m3_i1_i2p1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p2_i1_i2p1,const REAL_SIMD_ARRAY f_i0p3_i1_i2p1,const REAL_SIMD_ARRAY f_i0m3_i1_i2p2,const REAL_SIMD_ARRAY f_i0m2_i1_i2p2,const REAL_SIMD_ARRAY f_i0m1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p2_i1_i2p2,const REAL_SIMD_ARRAY f_i0p3_i1_i2p2,const REAL_SIMD_ARRAY f_i0m3_i1_i2p3,const REAL_SIMD_ARRAY f_i0m2_i1_i2p3,const REAL_SIMD_ARRAY f_i0m1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p2_i1_i2p3,const REAL_SIMD_ARRAY f_i0p3_i1_i2p3) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_3600 = 1.0/3600.0;
  const REAL_SIMD_ARRAY _Rational_1_3600 = ConstSIMD(tmp_Rational_1_3600);

  const double tmp_Rational_1_400 = 1.0/400.0;
  const REAL_SIMD_ARRAY _Rational_1_400 = ConstSIMD(tmp_Rational_1_400);

  const double tmp_Rational_1_80 = 1.0/80.0;
  const REAL_SIMD_ARRAY _Rational_1_80 = ConstSIMD(tmp_Rational_1_80);

  const double tmp_Rational_9_16 = 9.0/16.0;
  const REAL_SIMD_ARRAY _Rational_9_16 = ConstSIMD(tmp_Rational_9_16);

  const double tmp_Rational_9_400 = 9.0/400.0;
  const REAL_SIMD_ARRAY _Rational_9_400 = ConstSIMD(tmp_Rational_9_400);

  const double tmp_Rational_9_80 = 9.0/80.0;
  const REAL_SIMD_ARRAY _Rational_9_80 = ConstSIMD(tmp_Rational_9_80);

  return MulSIMD(invdx0, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_80, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2p1, f_i0m3_i1_i2p1), AddSIMD(AddSIMD(f_i0m3_i1_i2m1, f_i0p1_i1_i2p3), SubSIMD(SubSIMD(f_i0m1_i1_i2m3, f_i0m1_i1_i2p3), f_i0p3_i1_i2m1))), f_i0p1_i1_i2m3), FusedMulAddSIMD(_Rational_9_16, AddSIMD(f_i0p1_i1_i2p1, SubSIMD(SubSIMD(f_i0m1_i1_i2m1, f_i0m1_i1_i2p1), f_i0p1_i1_i2m1)), FusedMulAddSIMD(_Rational_9_400, AddSIMD(f_i0p2_i1_i2p2, SubSIMD(SubSIMD(f_i0m2_i1_i2m2, f_i0m2_i1_i2p2), f_i0p2_i1_i2m2)), FusedMulAddSIMD(_Rational_9_80, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1_i2m1, f_i0m2_i1_i2m1), AddSIMD(AddSIMD(f_i0m2_i1_i2p1, f_i0p1_i1_i2m2), SubSIMD(SubSIMD(f_i0m1_i1_i2p2, f_i0m1_i1_i2m2), f_i0p2_i1_i2p1))), f_i0p1_i1_i2p2), FusedMulAddSIMD(_Rational_1_3600, AddSIMD(f_i0p3_i1_i2p3, SubSIMD(SubSIMD(f_i0m3_i1_i2m3, f_i0m3_i1_i2p3), f_i0p3_i1_i2m3)), MulSIMD(_Rational_1_400, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2m2, f_i0m3_i1_i2m2), AddSIMD(AddSIMD(f_i0m3_i1_i2p2, f_i0p2_i1_i2m3), SubSIMD(SubSIMD(f_i0m2_i1_i2p3, f_i0m2_i1_i2m3), f_i0p3_i1_i2p2))), f_i0p2_i1_i2p3)))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD11(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_90 = 1.0/90.0;
  const REAL_SIMD_ARRAY _Rational_1_90 = ConstSIMD(tmp_Rational_1_90);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_3_20 = 3.0/20.0;
  const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

  const double tmp_Rational_49_18 = 49.0/18.0;
  const REAL_SIMD_ARRAY _Rational_49_18 = ConstSIMD(tmp_Rational_49_18);

  return MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(_Rational_3_2, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), FusedMulAddSIMD(_Rational_3_20, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2)), FusedMulSubSIMD(_Rational_1_90, AddSIMD(f_i0_i1m3_i2, f_i0_i1p3_i2), MulSIMD(_Rational_49_18, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD12(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1m3_i2m3,const REAL_SIMD_ARRAY f_i0_i1m2_i2m3,const REAL_SIMD_ARRAY f_i0_i1m1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p2_i2m3,const REAL_SIMD_ARRAY f_i0_i1p3_i2m3,const REAL_SIMD_ARRAY f_i0_i1m3_i2m2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p2_i2m2,const REAL_SIMD_ARRAY f_i0_i1p3_i2m2,const REAL_SIMD_ARRAY f_i0_i1m3_i2m1,const REAL_SIMD_ARRAY f_i0_i1m2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p2_i2m1,const REAL_SIMD_ARRAY f_i0_i1p3_i2m1,const REAL_SIMD_ARRAY f_i0_i1m3_i2p1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p2_i2p1,const REAL_SIMD_ARRAY f_i0_i1p3_i2p1,const REAL_SIMD_ARRAY f_i0_i1m3_i2p2,const REAL_SIMD_ARRAY f_i0_i1m2_i2p2,const REAL_SIMD_ARRAY f_i0_i1m1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p2_i2p2,const REAL_SIMD_ARRAY f_i0_i1p3_i2p2,const REAL_SIMD_ARRAY f_i0_i1m3_i2p3,const REAL_SIMD_ARRAY f_i0_i1m2_i2p3,const REAL_SIMD_ARRAY f_i0_i1m1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p2_i2p3,const REAL_SIMD_ARRAY f_i0_i1p3_i2p3) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_3600 = 1.0/3600.0;
  const REAL_SIMD_ARRAY _Rational_1_3600 = ConstSIMD(tmp_Rational_1_3600);

  const double tmp_Rational_1_400 = 1.0/400.0;
  const REAL_SIMD_ARRAY _Rational_1_400 = ConstSIMD(tmp_Rational_1_400);

  const double tmp_Rational_1_80 = 1.0/80.0;
  const REAL_SIMD_ARRAY _Rational_1_80 = ConstSIMD(tmp_Rational_1_80);

  const double tmp_Rational_9_16 = 9.0/16.0;
  const REAL_SIMD_ARRAY _Rational_9_16 = ConstSIMD(tmp_Rational_9_16);

  const double tmp_Rational_9_400 = 9.0/400.0;
  const REAL_SIMD_ARRAY _Rational_9_400 = ConstSIMD(tmp_Rational_9_400);

  const double tmp_Rational_9_80 = 9.0/80.0;
  const REAL_SIMD_ARRAY _Rational_9_80 = ConstSIMD(tmp_Rational_9_80);

  return MulSIMD(invdx1, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_80, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2p1, f_i0_i1m3_i2p1), AddSIMD(AddSIMD(f_i0_i1m3_i2m1, f_i0_i1p1_i2p3), SubSIMD(SubSIMD(f_i0_i1m1_i2m3, f_i0_i1m1_i2p3), f_i0_i1p3_i2m1))), f_i0_i1p1_i2m3), FusedMulAddSIMD(_Rational_9_16, AddSIMD(f_i0_i1p1_i2p1, SubSIMD(SubSIMD(f_i0_i1m1_i2m1, f_i0_i1m1_i2p1), f_i0_i1p1_i2m1)), FusedMulAddSIMD(_Rational_9_400, AddSIMD(f_i0_i1p2_i2p2, SubSIMD(SubSIMD(f_i0_i1m2_i2m2, f_i0_i1m2_i2p2), f_i0_i1p2_i2m2)), FusedMulAddSIMD(_Rational_9_80, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p2_i2m1, f_i0_i1m2_i2m1), AddSIMD(AddSIMD(f_i0_i1m2_i2p1, f_i0_i1p1_i2m2), SubSIMD(SubSIMD(f_i0_i1m1_i2p2, f_i0_i1m1_i2m2), f_i0_i1p2_i2p1))), f_i0_i1p1_i2p2), FusedMulAddSIMD(_Rational_1_3600, AddSIMD(f_i0_i1p3_i2p3, SubSIMD(SubSIMD(f_i0_i1m3_i2m3, f_i0_i1m3_i2p3), f_i0_i1p3_i2m3)), MulSIMD(_Rational_1_400, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2m2, f_i0_i1m3_i2m2), AddSIMD(AddSIMD(f_i0_i1m3_i2p2, f_i0_i1p2_i2m3), SubSIMD(SubSIMD(f_i0_i1m2_i2p3, f_i0_i1m2_i2m3), f_i0_i1p3_i2p2))), f_i0_i1p2_i2p3)))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD22(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_90 = 1.0/90.0;
  const REAL_SIMD_ARRAY _Rational_1_90 = ConstSIMD(tmp_Rational_1_90);

  const double tmp_Rational_3_2 = 3.0/2.0;
  const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

  const double tmp_Rational_3_20 = 3.0/20.0;
  const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

  const double tmp_Rational_49_18 = 49.0/18.0;
  const REAL_SIMD_ARRAY _Rational_49_18 = ConstSIMD(tmp_Rational_49_18);

  return MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(_Rational_3_2, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), FusedMulAddSIMD(_Rational_3_20, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2)), FusedMulSubSIMD(_Rational_1_90, AddSIMD(f_i0_i1_i2m3, f_i0_i1_i2p3), MulSIMD(_Rational_49_18, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dKOD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m5_i1_i2,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2,const REAL_SIMD_ARRAY f_i0p5_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_105_512 = 105.0/512.0;
  const REAL_SIMD_ARRAY _Rational_105_512 = ConstSIMD(tmp_Rational_105_512);

  const double tmp_Rational_15_128 = 15.0/128.0;
  const REAL_SIMD_ARRAY _Rational_15_128 = ConstSIMD(tmp_Rational_15_128);

  const double tmp_Rational_1_1024 = 1.0/1024.0;
  const REAL_SIMD_ARRAY _Rational_1_1024 = ConstSIMD(tmp_Rational_1_1024);

  const double tmp_Rational_45_1024 = 45.0/1024.0;
  const REAL_SIMD_ARRAY _Rational_45_1024 = ConstSIMD(tmp_Rational_45_1024);

  const double tmp_Rational_5_512 = 5.0/512.0;
  const REAL_SIMD_ARRAY _Rational_5_512 = ConstSIMD(tmp_Rational_5_512);

  const double tmp_Rational_63_256 = 63.0/256.0;
  const REAL_SIMD_ARRAY _Rational_63_256 = ConstSIMD(tmp_Rational_63_256);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_45_1024, AddSIMD(f_i0m3_i1_i2, f_i0p3_i1_i2), FusedMulAddSIMD(_Rational_15_128, MulSIMD(_NegativeOne_, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2)), FusedMulAddSIMD(_Rational_1_1024, AddSIMD(f_i0m5_i1_i2, f_i0p5_i1_i2), FusedMulAddSIMD(_Rational_5_512, MulSIMD(_NegativeOne_, AddSIMD(f_i0p4_i1_i2, f_i0m4_i1_i2)), FusedMulSubSIMD(_Rational_105_512, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), MulSIMD(_Rational_63_256, f)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dKOD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m5_i2,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2,const REAL_SIMD_ARRAY f_i0_i1p5_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_105_512 = 105.0/512.0;
  const REAL_SIMD_ARRAY _Rational_105_512 = ConstSIMD(tmp_Rational_105_512);

  const double tmp_Rational_15_128 = 15.0/128.0;
  const REAL_SIMD_ARRAY _Rational_15_128 = ConstSIMD(tmp_Rational_15_128);

  const double tmp_Rational_1_1024 = 1.0/1024.0;
  const REAL_SIMD_ARRAY _Rational_1_1024 = ConstSIMD(tmp_Rational_1_1024);

  const double tmp_Rational_45_1024 = 45.0/1024.0;
  const REAL_SIMD_ARRAY _Rational_45_1024 = ConstSIMD(tmp_Rational_45_1024);

  const double tmp_Rational_5_512 = 5.0/512.0;
  const REAL_SIMD_ARRAY _Rational_5_512 = ConstSIMD(tmp_Rational_5_512);

  const double tmp_Rational_63_256 = 63.0/256.0;
  const REAL_SIMD_ARRAY _Rational_63_256 = ConstSIMD(tmp_Rational_63_256);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_45_1024, AddSIMD(f_i0_i1m3_i2, f_i0_i1p3_i2), FusedMulAddSIMD(_Rational_15_128, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2)), FusedMulAddSIMD(_Rational_1_1024, AddSIMD(f_i0_i1m5_i2, f_i0_i1p5_i2), FusedMulAddSIMD(_Rational_5_512, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p4_i2, f_i0_i1m4_i2)), FusedMulSubSIMD(_Rational_105_512, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), MulSIMD(_Rational_63_256, f)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for Kreiss-Oliger derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dKOD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m5,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4,const REAL_SIMD_ARRAY f_i0_i1_i2p5) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_105_512 = 105.0/512.0;
  const REAL_SIMD_ARRAY _Rational_105_512 = ConstSIMD(tmp_Rational_105_512);

  const double tmp_Rational_15_128 = 15.0/128.0;
  const REAL_SIMD_ARRAY _Rational_15_128 = ConstSIMD(tmp_Rational_15_128);

  const double tmp_Rational_1_1024 = 1.0/1024.0;
  const REAL_SIMD_ARRAY _Rational_1_1024 = ConstSIMD(tmp_Rational_1_1024);

  const double tmp_Rational_45_1024 = 45.0/1024.0;
  const REAL_SIMD_ARRAY _Rational_45_1024 = ConstSIMD(tmp_Rational_45_1024);

  const double tmp_Rational_5_512 = 5.0/512.0;
  const REAL_SIMD_ARRAY _Rational_5_512 = ConstSIMD(tmp_Rational_5_512);

  const double tmp_Rational_63_256 = 63.0/256.0;
  const REAL_SIMD_ARRAY _Rational_63_256 = ConstSIMD(tmp_Rational_63_256);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_45_1024, AddSIMD(f_i0_i1_i2m3, f_i0_i1_i2p3), FusedMulAddSIMD(_Rational_15_128, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2)), FusedMulAddSIMD(_Rational_1_1024, AddSIMD(f_i0_i1_i2m5, f_i0_i1_i2p5), FusedMulAddSIMD(_Rational_5_512, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p4, f_i0_i1_i2m4)), FusedMulSubSIMD(_Rational_105_512, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), MulSIMD(_Rational_63_256, f)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_ddnD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m5_i1_i2,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_14 = 1.0/14.0;
  const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

  const double tmp_Rational_1_168 = 1.0/168.0;
  const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_28 = 1.0/28.0;
  const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_6 = 1.0/6.0;
  const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

  const double tmp_Rational_5_4 = 5.0/4.0;
  const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

  const double tmp_Rational_9_20 = 9.0/20.0;
  const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

  return MulSIMD(invdx0, SubSIMD(FusedMulAddSIMD(_Rational_9_20, f, SubSIMD(FusedMulAddSIMD(_Rational_1_2, AddSIMD(f_i0m2_i1_i2, f_i0p1_i1_i2), FusedMulAddSIMD(_Rational_1_28, f_i0m4_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_168, f_i0p3_i1_i2, MulSIMD(_Rational_1_14, f_i0p2_i1_i2)), MulSIMD(_Rational_5_4, f_i0m1_i1_i2)))), MulSIMD(_Rational_1_280, f_i0m5_i1_i2))), MulSIMD(_Rational_1_6, f_i0m3_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_ddnD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m5_i2,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_14 = 1.0/14.0;
  const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

  const double tmp_Rational_1_168 = 1.0/168.0;
  const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_28 = 1.0/28.0;
  const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_6 = 1.0/6.0;
  const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

  const double tmp_Rational_5_4 = 5.0/4.0;
  const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

  const double tmp_Rational_9_20 = 9.0/20.0;
  const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

  return MulSIMD(invdx1, SubSIMD(FusedMulAddSIMD(_Rational_9_20, f, SubSIMD(FusedMulAddSIMD(_Rational_1_2, AddSIMD(f_i0_i1m2_i2, f_i0_i1p1_i2), FusedMulAddSIMD(_Rational_1_28, f_i0_i1m4_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_168, f_i0_i1p3_i2, MulSIMD(_Rational_1_14, f_i0_i1p2_i2)), MulSIMD(_Rational_5_4, f_i0_i1m1_i2)))), MulSIMD(_Rational_1_280, f_i0_i1m5_i2))), MulSIMD(_Rational_1_6, f_i0_i1m3_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_ddnD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m5,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_14 = 1.0/14.0;
  const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

  const double tmp_Rational_1_168 = 1.0/168.0;
  const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_28 = 1.0/28.0;
  const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_6 = 1.0/6.0;
  const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

  const double tmp_Rational_5_4 = 5.0/4.0;
  const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

  const double tmp_Rational_9_20 = 9.0/20.0;
  const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

  return MulSIMD(invdx2, SubSIMD(FusedMulAddSIMD(_Rational_9_20, f, SubSIMD(FusedMulAddSIMD(_Rational_1_2, AddSIMD(f_i0_i1_i2m2, f_i0_i1_i2p1), FusedMulAddSIMD(_Rational_1_28, f_i0_i1_i2m4, SubSIMD(FusedMulSubSIMD(_Rational_1_168, f_i0_i1_i2p3, MulSIMD(_Rational_1_14, f_i0_i1_i2p2)), MulSIMD(_Rational_5_4, f_i0_i1_i2m1)))), MulSIMD(_Rational_1_280, f_i0_i1_i2m5))), MulSIMD(_Rational_1_6, f_i0_i1_i2m3)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dupD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2,const REAL_SIMD_ARRAY f_i0p5_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_14 = 1.0/14.0;
  const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

  const double tmp_Rational_1_168 = 1.0/168.0;
  const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_28 = 1.0/28.0;
  const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_6 = 1.0/6.0;
  const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

  const double tmp_Rational_5_4 = 5.0/4.0;
  const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

  const double tmp_Rational_9_20 = 9.0/20.0;
  const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

  return MulSIMD(invdx0, SubSIMD(FusedMulAddSIMD(_Rational_1_6, f_i0p3_i1_i2, FusedMulAddSIMD(_Rational_5_4, f_i0p1_i1_i2, FusedMulAddSIMD(_Rational_1_2, MulSIMD(_NegativeOne_, AddSIMD(f_i0p2_i1_i2, f_i0m1_i1_i2)), FusedMulAddSIMD(_Rational_1_280, f_i0p5_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_14, f_i0m2_i1_i2, MulSIMD(_Rational_1_168, f_i0m3_i1_i2)), MulSIMD(_Rational_9_20, f)))))), MulSIMD(_Rational_1_28, f_i0p4_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dupD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2,const REAL_SIMD_ARRAY f_i0_i1p5_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_14 = 1.0/14.0;
  const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

  const double tmp_Rational_1_168 = 1.0/168.0;
  const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_28 = 1.0/28.0;
  const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_6 = 1.0/6.0;
  const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

  const double tmp_Rational_5_4 = 5.0/4.0;
  const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

  const double tmp_Rational_9_20 = 9.0/20.0;
  const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

  return MulSIMD(invdx1, SubSIMD(FusedMulAddSIMD(_Rational_1_6, f_i0_i1p3_i2, FusedMulAddSIMD(_Rational_5_4, f_i0_i1p1_i2, FusedMulAddSIMD(_Rational_1_2, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p2_i2, f_i0_i1m1_i2)), FusedMulAddSIMD(_Rational_1_280, f_i0_i1p5_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_14, f_i0_i1m2_i2, MulSIMD(_Rational_1_168, f_i0_i1m3_i2)), MulSIMD(_Rational_9_20, f)))))), MulSIMD(_Rational_1_28, f_i0_i1p4_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dupD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4,const REAL_SIMD_ARRAY f_i0_i1_i2p5) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_14 = 1.0/14.0;
  const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

  const double tmp_Rational_1_168 = 1.0/168.0;
  const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

  const double tmp_Rational_1_2 = 1.0/2.0;
  const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

  const double tmp_Rational_1_28 = 1.0/28.0;
  const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_6 = 1.0/6.0;
  const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

  const double tmp_Rational_5_4 = 5.0/4.0;
  const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

  const double tmp_Rational_9_20 = 9.0/20.0;
  const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

  return MulSIMD(invdx2, SubSIMD(FusedMulAddSIMD(_Rational_1_6, f_i0_i1_i2p3, FusedMulAddSIMD(_Rational_5_4, f_i0_i1_i2p1, FusedMulAddSIMD(_Rational_1_2, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m1)), FusedMulAddSIMD(_Rational_1_280, f_i0_i1_i2p5, SubSIMD(FusedMulSubSIMD(_Rational_1_14, f_i0_i1_i2m2, MulSIMD(_Rational_1_168, f_i0_i1_i2m3)), MulSIMD(_Rational_9_20, f)))))), MulSIMD(_Rational_1_28, f_i0_i1_i2p4)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_5 = 1.0/5.0;
  const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

  const double tmp_Rational_4_105 = 4.0/105.0;
  const REAL_SIMD_ARRAY _Rational_4_105 = ConstSIMD(tmp_Rational_4_105);

  const double tmp_Rational_4_5 = 4.0/5.0;
  const REAL_SIMD_ARRAY _Rational_4_5 = ConstSIMD(tmp_Rational_4_5);

  return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_4_105, SubSIMD(f_i0p3_i1_i2, f_i0m3_i1_i2), FusedMulAddSIMD(_Rational_4_5, SubSIMD(f_i0p1_i1_i2, f_i0m1_i1_i2), FusedMulAddSIMD(_Rational_1_280, SubSIMD(f_i0m4_i1_i2, f_i0p4_i1_i2), MulSIMD(_Rational_1_5, SubSIMD(f_i0m2_i1_i2, f_i0p2_i1_i2))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_5 = 1.0/5.0;
  const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

  const double tmp_Rational_4_105 = 4.0/105.0;
  const REAL_SIMD_ARRAY _Rational_4_105 = ConstSIMD(tmp_Rational_4_105);

  const double tmp_Rational_4_5 = 4.0/5.0;
  const REAL_SIMD_ARRAY _Rational_4_5 = ConstSIMD(tmp_Rational_4_5);

  return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_4_105, SubSIMD(f_i0_i1p3_i2, f_i0_i1m3_i2), FusedMulAddSIMD(_Rational_4_5, SubSIMD(f_i0_i1p1_i2, f_i0_i1m1_i2), FusedMulAddSIMD(_Rational_1_280, SubSIMD(f_i0_i1m4_i2, f_i0_i1p4_i2), MulSIMD(_Rational_1_5, SubSIMD(f_i0_i1m2_i2, f_i0_i1p2_i2))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_280 = 1.0/280.0;
  const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

  const double tmp_Rational_1_5 = 1.0/5.0;
  const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

  const double tmp_Rational_4_105 = 4.0/105.0;
  const REAL_SIMD_ARRAY _Rational_4_105 = ConstSIMD(tmp_Rational_4_105);

  const double tmp_Rational_4_5 = 4.0/5.0;
  const REAL_SIMD_ARRAY _Rational_4_5 = ConstSIMD(tmp_Rational_4_5);

  return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_4_105, SubSIMD(f_i0_i1_i2p3, f_i0_i1_i2m3), FusedMulAddSIMD(_Rational_4_5, SubSIMD(f_i0_i1_i2p1, f_i0_i1_i2m1), FusedMulAddSIMD(_Rational_1_280, SubSIMD(f_i0_i1_i2m4, f_i0_i1_i2p4), MulSIMD(_Rational_1_5, SubSIMD(f_i0_i1_i2m2, f_i0_i1_i2p2))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD00(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_5 = 1.0/5.0;
  const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

  const double tmp_Rational_1_560 = 1.0/560.0;
  const REAL_SIMD_ARRAY _Rational_1_560 = ConstSIMD(tmp_Rational_1_560);

  const double tmp_Rational_205_72 = 205.0/72.0;
  const REAL_SIMD_ARRAY _Rational_205_72 = ConstSIMD(tmp_Rational_205_72);

  const double tmp_Rational_8_315 = 8.0/315.0;
  const REAL_SIMD_ARRAY _Rational_8_315 = ConstSIMD(tmp_Rational_8_315);

  const double tmp_Rational_8_5 = 8.0/5.0;
  const REAL_SIMD_ARRAY _Rational_8_5 = ConstSIMD(tmp_Rational_8_5);

  return MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(_Rational_1_560, MulSIMD(_NegativeOne_, AddSIMD(f_i0p4_i1_i2, f_i0m4_i1_i2)), FusedMulAddSIMD(_Rational_8_315, AddSIMD(f_i0m3_i1_i2, f_i0p3_i1_i2), FusedMulAddSIMD(_Rational_8_5, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), NegFusedMulSubSIMD(_Rational_1_5, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2), MulSIMD(_Rational_205_72, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD01(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0m4_i1m4_i2,const REAL_SIMD_ARRAY f_i0m3_i1m4_i2,const REAL_SIMD_ARRAY f_i0m2_i1m4_i2,const REAL_SIMD_ARRAY f_i0m1_i1m4_i2,const REAL_SIMD_ARRAY f_i0p1_i1m4_i2,const REAL_SIMD_ARRAY f_i0p2_i1m4_i2,const REAL_SIMD_ARRAY f_i0p3_i1m4_i2,const REAL_SIMD_ARRAY f_i0p4_i1m4_i2,const REAL_SIMD_ARRAY f_i0m4_i1m3_i2,const REAL_SIMD_ARRAY f_i0m3_i1m3_i2,const REAL_SIMD_ARRAY f_i0m2_i1m3_i2,const REAL_SIMD_ARRAY f_i0m1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p2_i1m3_i2,const REAL_SIMD_ARRAY f_i0p3_i1m3_i2,const REAL_SIMD_ARRAY f_i0p4_i1m3_i2,const REAL_SIMD_ARRAY f_i0m4_i1m2_i2,const REAL_SIMD_ARRAY f_i0m3_i1m2_i2,const REAL_SIMD_ARRAY f_i0m2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p2_i1m2_i2,const REAL_SIMD_ARRAY f_i0p3_i1m2_i2,const REAL_SIMD_ARRAY f_i0p4_i1m2_i2,const REAL_SIMD_ARRAY f_i0m4_i1m1_i2,const REAL_SIMD_ARRAY f_i0m3_i1m1_i2,const REAL_SIMD_ARRAY f_i0m2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p2_i1m1_i2,const REAL_SIMD_ARRAY f_i0p3_i1m1_i2,const REAL_SIMD_ARRAY f_i0p4_i1m1_i2,const REAL_SIMD_ARRAY f_i0m4_i1p1_i2,const REAL_SIMD_ARRAY f_i0m3_i1p1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p2_i1p1_i2,const REAL_SIMD_ARRAY f_i0p3_i1p1_i2,const REAL_SIMD_ARRAY f_i0p4_i1p1_i2,const REAL_SIMD_ARRAY f_i0m4_i1p2_i2,const REAL_SIMD_ARRAY f_i0m3_i1p2_i2,const REAL_SIMD_ARRAY f_i0m2_i1p2_i2,const REAL_SIMD_ARRAY f_i0m1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p2_i1p2_i2,const REAL_SIMD_ARRAY f_i0p3_i1p2_i2,const REAL_SIMD_ARRAY f_i0p4_i1p2_i2,const REAL_SIMD_ARRAY f_i0m4_i1p3_i2,const REAL_SIMD_ARRAY f_i0m3_i1p3_i2,const REAL_SIMD_ARRAY f_i0m2_i1p3_i2,const REAL_SIMD_ARRAY f_i0m1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p2_i1p3_i2,const REAL_SIMD_ARRAY f_i0p3_i1p3_i2,const REAL_SIMD_ARRAY f_i0p4_i1p3_i2,const REAL_SIMD_ARRAY f_i0m4_i1p4_i2,const REAL_SIMD_ARRAY f_i0m3_i1p4_i2,const REAL_SIMD_ARRAY f_i0m2_i1p4_i2,const REAL_SIMD_ARRAY f_i0m1_i1p4_i2,const REAL_SIMD_ARRAY f_i0p1_i1p4_i2,const REAL_SIMD_ARRAY f_i0p2_i1p4_i2,const REAL_SIMD_ARRAY f_i0p3_i1p4_i2,const REAL_SIMD_ARRAY f_i0p4_i1p4_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_16_11025 = 16.0/11025.0;
  const REAL_SIMD_ARRAY _Rational_16_11025 = ConstSIMD(tmp_Rational_16_11025);

  const double tmp_Rational_16_25 = 16.0/25.0;
  const REAL_SIMD_ARRAY _Rational_16_25 = ConstSIMD(tmp_Rational_16_25);

  const double tmp_Rational_16_525 = 16.0/525.0;
  const REAL_SIMD_ARRAY _Rational_16_525 = ConstSIMD(tmp_Rational_16_525);

  const double tmp_Rational_1_1400 = 1.0/1400.0;
  const REAL_SIMD_ARRAY _Rational_1_1400 = ConstSIMD(tmp_Rational_1_1400);

  const double tmp_Rational_1_25 = 1.0/25.0;
  const REAL_SIMD_ARRAY _Rational_1_25 = ConstSIMD(tmp_Rational_1_25);

  const double tmp_Rational_1_350 = 1.0/350.0;
  const REAL_SIMD_ARRAY _Rational_1_350 = ConstSIMD(tmp_Rational_1_350);

  const double tmp_Rational_1_7350 = 1.0/7350.0;
  const REAL_SIMD_ARRAY _Rational_1_7350 = ConstSIMD(tmp_Rational_1_7350);

  const double tmp_Rational_1_78400 = 1.0/78400.0;
  const REAL_SIMD_ARRAY _Rational_1_78400 = ConstSIMD(tmp_Rational_1_78400);

  const double tmp_Rational_4_25 = 4.0/25.0;
  const REAL_SIMD_ARRAY _Rational_4_25 = ConstSIMD(tmp_Rational_4_25);

  const double tmp_Rational_4_525 = 4.0/525.0;
  const REAL_SIMD_ARRAY _Rational_4_525 = ConstSIMD(tmp_Rational_4_525);

  return MulSIMD(invdx0, MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_7350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1m3_i2, f_i0m4_i1m3_i2), AddSIMD(AddSIMD(f_i0m4_i1p3_i2, f_i0p3_i1m4_i2), SubSIMD(SubSIMD(f_i0m3_i1p4_i2, f_i0m3_i1m4_i2), f_i0p4_i1p3_i2))), f_i0p3_i1p4_i2), FusedMulAddSIMD(_Rational_1_78400, AddSIMD(f_i0p4_i1p4_i2, SubSIMD(SubSIMD(f_i0m4_i1m4_i2, f_i0m4_i1p4_i2), f_i0p4_i1m4_i2)), FusedMulAddSIMD(_Rational_1_25, AddSIMD(f_i0p2_i1p2_i2, SubSIMD(SubSIMD(f_i0m2_i1m2_i2, f_i0m2_i1p2_i2), f_i0p2_i1m2_i2)), FusedMulAddSIMD(_Rational_1_350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1m1_i2, f_i0m4_i1m1_i2), AddSIMD(AddSIMD(f_i0m4_i1p1_i2, f_i0p1_i1m4_i2), SubSIMD(SubSIMD(f_i0m1_i1p4_i2, f_i0m1_i1m4_i2), f_i0p4_i1p1_i2))), f_i0p1_i1p4_i2), FusedMulAddSIMD(_Rational_16_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1p1_i2, f_i0m3_i1p1_i2), AddSIMD(AddSIMD(f_i0m3_i1m1_i2, f_i0p1_i1p3_i2), SubSIMD(SubSIMD(f_i0m1_i1m3_i2, f_i0m1_i1p3_i2), f_i0p3_i1m1_i2))), f_i0p1_i1m3_i2), FusedMulAddSIMD(_Rational_1_1400, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1p2_i2, f_i0m4_i1p2_i2), AddSIMD(AddSIMD(f_i0m4_i1m2_i2, f_i0p2_i1p4_i2), SubSIMD(SubSIMD(f_i0m2_i1m4_i2, f_i0m2_i1p4_i2), f_i0p4_i1m2_i2))), f_i0p2_i1m4_i2), FusedMulAddSIMD(_Rational_4_25, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1m1_i2, f_i0m2_i1m1_i2), AddSIMD(AddSIMD(f_i0m2_i1p1_i2, f_i0p1_i1m2_i2), SubSIMD(SubSIMD(f_i0m1_i1p2_i2, f_i0m1_i1m2_i2), f_i0p2_i1p1_i2))), f_i0p1_i1p2_i2), FusedMulAddSIMD(_Rational_4_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1m2_i2, f_i0m3_i1m2_i2), AddSIMD(AddSIMD(f_i0m3_i1p2_i2, f_i0p2_i1m3_i2), SubSIMD(SubSIMD(f_i0m2_i1p3_i2, f_i0m2_i1m3_i2), f_i0p3_i1p2_i2))), f_i0p2_i1p3_i2), FusedMulAddSIMD(_Rational_16_11025, AddSIMD(f_i0p3_i1p3_i2, SubSIMD(SubSIMD(f_i0m3_i1m3_i2, f_i0m3_i1p3_i2), f_i0p3_i1m3_i2)), MulSIMD(_Rational_16_25, AddSIMD(f_i0p1_i1p1_i2, SubSIMD(SubSIMD(f_i0m1_i1m1_i2, f_i0m1_i1p1_i2), f_i0p1_i1m1_i2))))))))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD02(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0m4_i1_i2m4,const REAL_SIMD_ARRAY f_i0m3_i1_i2m4,const REAL_SIMD_ARRAY f_i0m2_i1_i2m4,const REAL_SIMD_ARRAY f_i0m1_i1_i2m4,const REAL_SIMD_ARRAY f_i0p1_i1_i2m4,const REAL_SIMD_ARRAY f_i0p2_i1_i2m4,const REAL_SIMD_ARRAY f_i0p3_i1_i2m4,const REAL_SIMD_ARRAY f_i0p4_i1_i2m4,const REAL_SIMD_ARRAY f_i0m4_i1_i2m3,const REAL_SIMD_ARRAY f_i0m3_i1_i2m3,const REAL_SIMD_ARRAY f_i0m2_i1_i2m3,const REAL_SIMD_ARRAY f_i0m1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p2_i1_i2m3,const REAL_SIMD_ARRAY f_i0p3_i1_i2m3,const REAL_SIMD_ARRAY f_i0p4_i1_i2m3,const REAL_SIMD_ARRAY f_i0m4_i1_i2m2,const REAL_SIMD_ARRAY f_i0m3_i1_i2m2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p2_i1_i2m2,const REAL_SIMD_ARRAY f_i0p3_i1_i2m2,const REAL_SIMD_ARRAY f_i0p4_i1_i2m2,const REAL_SIMD_ARRAY f_i0m4_i1_i2m1,const REAL_SIMD_ARRAY f_i0m3_i1_i2m1,const REAL_SIMD_ARRAY f_i0m2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p2_i1_i2m1,const REAL_SIMD_ARRAY f_i0p3_i1_i2m1,const REAL_SIMD_ARRAY f_i0p4_i1_i2m1,const REAL_SIMD_ARRAY f_i0m4_i1_i2p1,const REAL_SIMD_ARRAY f_i0m3_i1_i2p1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p2_i1_i2p1,const REAL_SIMD_ARRAY f_i0p3_i1_i2p1,const REAL_SIMD_ARRAY f_i0p4_i1_i2p1,const REAL_SIMD_ARRAY f_i0m4_i1_i2p2,const REAL_SIMD_ARRAY f_i0m3_i1_i2p2,const REAL_SIMD_ARRAY f_i0m2_i1_i2p2,const REAL_SIMD_ARRAY f_i0m1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p2_i1_i2p2,const REAL_SIMD_ARRAY f_i0p3_i1_i2p2,const REAL_SIMD_ARRAY f_i0p4_i1_i2p2,const REAL_SIMD_ARRAY f_i0m4_i1_i2p3,const REAL_SIMD_ARRAY f_i0m3_i1_i2p3,const REAL_SIMD_ARRAY f_i0m2_i1_i2p3,const REAL_SIMD_ARRAY f_i0m1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p2_i1_i2p3,const REAL_SIMD_ARRAY f_i0p3_i1_i2p3,const REAL_SIMD_ARRAY f_i0p4_i1_i2p3,const REAL_SIMD_ARRAY f_i0m4_i1_i2p4,const REAL_SIMD_ARRAY f_i0m3_i1_i2p4,const REAL_SIMD_ARRAY f_i0m2_i1_i2p4,const REAL_SIMD_ARRAY f_i0m1_i1_i2p4,const REAL_SIMD_ARRAY f_i0p1_i1_i2p4,const REAL_SIMD_ARRAY f_i0p2_i1_i2p4,const REAL_SIMD_ARRAY f_i0p3_i1_i2p4,const REAL_SIMD_ARRAY f_i0p4_i1_i2p4) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_16_11025 = 16.0/11025.0;
  const REAL_SIMD_ARRAY _Rational_16_11025 = ConstSIMD(tmp_Rational_16_11025);

  const double tmp_Rational_16_25 = 16.0/25.0;
  const REAL_SIMD_ARRAY _Rational_16_25 = ConstSIMD(tmp_Rational_16_25);

  const double tmp_Rational_16_525 = 16.0/525.0;
  const REAL_SIMD_ARRAY _Rational_16_525 = ConstSIMD(tmp_Rational_16_525);

  const double tmp_Rational_1_1400 = 1.0/1400.0;
  const REAL_SIMD_ARRAY _Rational_1_1400 = ConstSIMD(tmp_Rational_1_1400);

  const double tmp_Rational_1_25 = 1.0/25.0;
  const REAL_SIMD_ARRAY _Rational_1_25 = ConstSIMD(tmp_Rational_1_25);

  const double tmp_Rational_1_350 = 1.0/350.0;
  const REAL_SIMD_ARRAY _Rational_1_350 = ConstSIMD(tmp_Rational_1_350);

  const double tmp_Rational_1_7350 = 1.0/7350.0;
  const REAL_SIMD_ARRAY _Rational_1_7350 = ConstSIMD(tmp_Rational_1_7350);

  const double tmp_Rational_1_78400 = 1.0/78400.0;
  const REAL_SIMD_ARRAY _Rational_1_78400 = ConstSIMD(tmp_Rational_1_78400);

  const double tmp_Rational_4_25 = 4.0/25.0;
  const REAL_SIMD_ARRAY _Rational_4_25 = ConstSIMD(tmp_Rational_4_25);

  const double tmp_Rational_4_525 = 4.0/525.0;
  const REAL_SIMD_ARRAY _Rational_4_525 = ConstSIMD(tmp_Rational_4_525);

  return MulSIMD(invdx0, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_7350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1_i2m3, f_i0m4_i1_i2m3), AddSIMD(AddSIMD(f_i0m4_i1_i2p3, f_i0p3_i1_i2m4), SubSIMD(SubSIMD(f_i0m3_i1_i2p4, f_i0m3_i1_i2m4), f_i0p4_i1_i2p3))), f_i0p3_i1_i2p4), FusedMulAddSIMD(_Rational_1_78400, AddSIMD(f_i0p4_i1_i2p4, SubSIMD(SubSIMD(f_i0m4_i1_i2m4, f_i0m4_i1_i2p4), f_i0p4_i1_i2m4)), FusedMulAddSIMD(_Rational_1_25, AddSIMD(f_i0p2_i1_i2p2, SubSIMD(SubSIMD(f_i0m2_i1_i2m2, f_i0m2_i1_i2p2), f_i0p2_i1_i2m2)), FusedMulAddSIMD(_Rational_1_350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1_i2m1, f_i0m4_i1_i2m1), AddSIMD(AddSIMD(f_i0m4_i1_i2p1, f_i0p1_i1_i2m4), SubSIMD(SubSIMD(f_i0m1_i1_i2p4, f_i0m1_i1_i2m4), f_i0p4_i1_i2p1))), f_i0p1_i1_i2p4), FusedMulAddSIMD(_Rational_16_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2p1, f_i0m3_i1_i2p1), AddSIMD(AddSIMD(f_i0m3_i1_i2m1, f_i0p1_i1_i2p3), SubSIMD(SubSIMD(f_i0m1_i1_i2m3, f_i0m1_i1_i2p3), f_i0p3_i1_i2m1))), f_i0p1_i1_i2m3), FusedMulAddSIMD(_Rational_1_1400, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1_i2p2, f_i0m4_i1_i2p2), AddSIMD(AddSIMD(f_i0m4_i1_i2m2, f_i0p2_i1_i2p4), SubSIMD(SubSIMD(f_i0m2_i1_i2m4, f_i0m2_i1_i2p4), f_i0p4_i1_i2m2))), f_i0p2_i1_i2m4), FusedMulAddSIMD(_Rational_4_25, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1_i2m1, f_i0m2_i1_i2m1), AddSIMD(AddSIMD(f_i0m2_i1_i2p1, f_i0p1_i1_i2m2), SubSIMD(SubSIMD(f_i0m1_i1_i2p2, f_i0m1_i1_i2m2), f_i0p2_i1_i2p1))), f_i0p1_i1_i2p2), FusedMulAddSIMD(_Rational_4_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2m2, f_i0m3_i1_i2m2), AddSIMD(AddSIMD(f_i0m3_i1_i2p2, f_i0p2_i1_i2m3), SubSIMD(SubSIMD(f_i0m2_i1_i2p3, f_i0m2_i1_i2m3), f_i0p3_i1_i2p2))), f_i0p2_i1_i2p3), FusedMulAddSIMD(_Rational_16_11025, AddSIMD(f_i0p3_i1_i2p3, SubSIMD(SubSIMD(f_i0m3_i1_i2m3, f_i0m3_i1_i2p3), f_i0p3_i1_i2m3)), MulSIMD(_Rational_16_25, AddSIMD(f_i0p1_i1_i2p1, SubSIMD(SubSIMD(f_i0m1_i1_i2m1, f_i0m1_i1_i2p1), f_i0p1_i1_i2m1))))))))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD11(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_5 = 1.0/5.0;
  const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

  const double tmp_Rational_1_560 = 1.0/560.0;
  const REAL_SIMD_ARRAY _Rational_1_560 = ConstSIMD(tmp_Rational_1_560);

  const double tmp_Rational_205_72 = 205.0/72.0;
  const REAL_SIMD_ARRAY _Rational_205_72 = ConstSIMD(tmp_Rational_205_72);

  const double tmp_Rational_8_315 = 8.0/315.0;
  const REAL_SIMD_ARRAY _Rational_8_315 = ConstSIMD(tmp_Rational_8_315);

  const double tmp_Rational_8_5 = 8.0/5.0;
  const REAL_SIMD_ARRAY _Rational_8_5 = ConstSIMD(tmp_Rational_8_5);

  return MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(_Rational_1_560, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p4_i2, f_i0_i1m4_i2)), FusedMulAddSIMD(_Rational_8_315, AddSIMD(f_i0_i1m3_i2, f_i0_i1p3_i2), FusedMulAddSIMD(_Rational_8_5, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), NegFusedMulSubSIMD(_Rational_1_5, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2), MulSIMD(_Rational_205_72, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD12(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1m4_i2m4,const REAL_SIMD_ARRAY f_i0_i1m3_i2m4,const REAL_SIMD_ARRAY f_i0_i1m2_i2m4,const REAL_SIMD_ARRAY f_i0_i1m1_i2m4,const REAL_SIMD_ARRAY f_i0_i1p1_i2m4,const REAL_SIMD_ARRAY f_i0_i1p2_i2m4,const REAL_SIMD_ARRAY f_i0_i1p3_i2m4,const REAL_SIMD_ARRAY f_i0_i1p4_i2m4,const REAL_SIMD_ARRAY f_i0_i1m4_i2m3,const REAL_SIMD_ARRAY f_i0_i1m3_i2m3,const REAL_SIMD_ARRAY f_i0_i1m2_i2m3,const REAL_SIMD_ARRAY f_i0_i1m1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p2_i2m3,const REAL_SIMD_ARRAY f_i0_i1p3_i2m3,const REAL_SIMD_ARRAY f_i0_i1p4_i2m3,const REAL_SIMD_ARRAY f_i0_i1m4_i2m2,const REAL_SIMD_ARRAY f_i0_i1m3_i2m2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p2_i2m2,const REAL_SIMD_ARRAY f_i0_i1p3_i2m2,const REAL_SIMD_ARRAY f_i0_i1p4_i2m2,const REAL_SIMD_ARRAY f_i0_i1m4_i2m1,const REAL_SIMD_ARRAY f_i0_i1m3_i2m1,const REAL_SIMD_ARRAY f_i0_i1m2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p2_i2m1,const REAL_SIMD_ARRAY f_i0_i1p3_i2m1,const REAL_SIMD_ARRAY f_i0_i1p4_i2m1,const REAL_SIMD_ARRAY f_i0_i1m4_i2p1,const REAL_SIMD_ARRAY f_i0_i1m3_i2p1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p2_i2p1,const REAL_SIMD_ARRAY f_i0_i1p3_i2p1,const REAL_SIMD_ARRAY f_i0_i1p4_i2p1,const REAL_SIMD_ARRAY f_i0_i1m4_i2p2,const REAL_SIMD_ARRAY f_i0_i1m3_i2p2,const REAL_SIMD_ARRAY f_i0_i1m2_i2p2,const REAL_SIMD_ARRAY f_i0_i1m1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p2_i2p2,const REAL_SIMD_ARRAY f_i0_i1p3_i2p2,const REAL_SIMD_ARRAY f_i0_i1p4_i2p2,const REAL_SIMD_ARRAY f_i0_i1m4_i2p3,const REAL_SIMD_ARRAY f_i0_i1m3_i2p3,const REAL_SIMD_ARRAY f_i0_i1m2_i2p3,const REAL_SIMD_ARRAY f_i0_i1m1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p2_i2p3,const REAL_SIMD_ARRAY f_i0_i1p3_i2p3,const REAL_SIMD_ARRAY f_i0_i1p4_i2p3,const REAL_SIMD_ARRAY f_i0_i1m4_i2p4,const REAL_SIMD_ARRAY f_i0_i1m3_i2p4,const REAL_SIMD_ARRAY f_i0_i1m2_i2p4,const REAL_SIMD_ARRAY f_i0_i1m1_i2p4,const REAL_SIMD_ARRAY f_i0_i1p1_i2p4,const REAL_SIMD_ARRAY f_i0_i1p2_i2p4,const REAL_SIMD_ARRAY f_i0_i1p3_i2p4,const REAL_SIMD_ARRAY f_i0_i1p4_i2p4) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_16_11025 = 16.0/11025.0;
  const REAL_SIMD_ARRAY _Rational_16_11025 = ConstSIMD(tmp_Rational_16_11025);

  const double tmp_Rational_16_25 = 16.0/25.0;
  const REAL_SIMD_ARRAY _Rational_16_25 = ConstSIMD(tmp_Rational_16_25);

  const double tmp_Rational_16_525 = 16.0/525.0;
  const REAL_SIMD_ARRAY _Rational_16_525 = ConstSIMD(tmp_Rational_16_525);

  const double tmp_Rational_1_1400 = 1.0/1400.0;
  const REAL_SIMD_ARRAY _Rational_1_1400 = ConstSIMD(tmp_Rational_1_1400);

  const double tmp_Rational_1_25 = 1.0/25.0;
  const REAL_SIMD_ARRAY _Rational_1_25 = ConstSIMD(tmp_Rational_1_25);

  const double tmp_Rational_1_350 = 1.0/350.0;
  const REAL_SIMD_ARRAY _Rational_1_350 = ConstSIMD(tmp_Rational_1_350);

  const double tmp_Rational_1_7350 = 1.0/7350.0;
  const REAL_SIMD_ARRAY _Rational_1_7350 = ConstSIMD(tmp_Rational_1_7350);

  const double tmp_Rational_1_78400 = 1.0/78400.0;
  const REAL_SIMD_ARRAY _Rational_1_78400 = ConstSIMD(tmp_Rational_1_78400);

  const double tmp_Rational_4_25 = 4.0/25.0;
  const REAL_SIMD_ARRAY _Rational_4_25 = ConstSIMD(tmp_Rational_4_25);

  const double tmp_Rational_4_525 = 4.0/525.0;
  const REAL_SIMD_ARRAY _Rational_4_525 = ConstSIMD(tmp_Rational_4_525);

  return MulSIMD(invdx1, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_7350, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p4_i2m3, f_i0_i1m4_i2m3), AddSIMD(AddSIMD(f_i0_i1m4_i2p3, f_i0_i1p3_i2m4), SubSIMD(SubSIMD(f_i0_i1m3_i2p4, f_i0_i1m3_i2m4), f_i0_i1p4_i2p3))), f_i0_i1p3_i2p4), FusedMulAddSIMD(_Rational_1_78400, AddSIMD(f_i0_i1p4_i2p4, SubSIMD(SubSIMD(f_i0_i1m4_i2m4, f_i0_i1m4_i2p4), f_i0_i1p4_i2m4)), FusedMulAddSIMD(_Rational_1_25, AddSIMD(f_i0_i1p2_i2p2, SubSIMD(SubSIMD(f_i0_i1m2_i2m2, f_i0_i1m2_i2p2), f_i0_i1p2_i2m2)), FusedMulAddSIMD(_Rational_1_350, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p4_i2m1, f_i0_i1m4_i2m1), AddSIMD(AddSIMD(f_i0_i1m4_i2p1, f_i0_i1p1_i2m4), SubSIMD(SubSIMD(f_i0_i1m1_i2p4, f_i0_i1m1_i2m4), f_i0_i1p4_i2p1))), f_i0_i1p1_i2p4), FusedMulAddSIMD(_Rational_16_525, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2p1, f_i0_i1m3_i2p1), AddSIMD(AddSIMD(f_i0_i1m3_i2m1, f_i0_i1p1_i2p3), SubSIMD(SubSIMD(f_i0_i1m1_i2m3, f_i0_i1m1_i2p3), f_i0_i1p3_i2m1))), f_i0_i1p1_i2m3), FusedMulAddSIMD(_Rational_1_1400, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p4_i2p2, f_i0_i1m4_i2p2), AddSIMD(AddSIMD(f_i0_i1m4_i2m2, f_i0_i1p2_i2p4), SubSIMD(SubSIMD(f_i0_i1m2_i2m4, f_i0_i1m2_i2p4), f_i0_i1p4_i2m2))), f_i0_i1p2_i2m4), FusedMulAddSIMD(_Rational_4_25, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p2_i2m1, f_i0_i1m2_i2m1), AddSIMD(AddSIMD(f_i0_i1m2_i2p1, f_i0_i1p1_i2m2), SubSIMD(SubSIMD(f_i0_i1m1_i2p2, f_i0_i1m1_i2m2), f_i0_i1p2_i2p1))), f_i0_i1p1_i2p2), FusedMulAddSIMD(_Rational_4_525, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2m2, f_i0_i1m3_i2m2), AddSIMD(AddSIMD(f_i0_i1m3_i2p2, f_i0_i1p2_i2m3), SubSIMD(SubSIMD(f_i0_i1m2_i2p3, f_i0_i1m2_i2m3), f_i0_i1p3_i2p2))), f_i0_i1p2_i2p3), FusedMulAddSIMD(_Rational_16_11025, AddSIMD(f_i0_i1p3_i2p3, SubSIMD(SubSIMD(f_i0_i1m3_i2m3, f_i0_i1m3_i2p3), f_i0_i1p3_i2m3)), MulSIMD(_Rational_16_25, AddSIMD(f_i0_i1p1_i2p1, SubSIMD(SubSIMD(f_i0_i1m1_i2m1, f_i0_i1m1_i2p1), f_i0_i1p1_i2m1))))))))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22 direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD22(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4) {

  const double tmp_NegativeOne_ = -1.0;
  const REAL_SIMD_ARRAY CCTK_ATTRIBUTE_UNUSED _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

  const double tmp_Rational_1_5 = 1.0/5.0;
  const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

  const double tmp_Rational_1_560 = 1.0/560.0;
  const REAL_SIMD_ARRAY _Rational_1_560 = ConstSIMD(tmp_Rational_1_560);

  const double tmp_Rational_205_72 = 205.0/72.0;
  const REAL_SIMD_ARRAY _Rational_205_72 = ConstSIMD(tmp_Rational_205_72);

  const double tmp_Rational_8_315 = 8.0/315.0;
  const REAL_SIMD_ARRAY _Rational_8_315 = ConstSIMD(tmp_Rational_8_315);

  const double tmp_Rational_8_5 = 8.0/5.0;
  const REAL_SIMD_ARRAY _Rational_8_5 = ConstSIMD(tmp_Rational_8_5);

  return MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(_Rational_1_560, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p4, f_i0_i1_i2m4)), FusedMulAddSIMD(_Rational_8_315, AddSIMD(f_i0_i1_i2m3, f_i0_i1_i2p3), FusedMulAddSIMD(_Rational_8_5, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), NegFusedMulSubSIMD(_Rational_1_5, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2), MulSIMD(_Rational_205_72, f))))));
}
#endif // #ifndef __FD_FUNCTIONS_H__
