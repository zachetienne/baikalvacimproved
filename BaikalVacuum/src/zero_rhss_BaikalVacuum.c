#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * Zero RHSs for NRPy+-generated thorn BaikalVacuum
 */
void zero_rhss_BaikalVacuum(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  #pragma omp parallel for
  for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
    for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
      for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
        aDD00_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        aDD01_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        aDD02_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        aDD11_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        aDD12_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        aDD22_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        alpha_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        betU0_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        betU1_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        betU2_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        cf_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        hDD00_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        hDD01_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        hDD02_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        hDD11_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        hDD12_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        hDD22_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        lambdaU0_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        lambdaU1_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        lambdaU2_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        trK_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        vetU0_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        vetU1_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
        vetU2_rhsGF[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;
      } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
    } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
  } // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)
}
