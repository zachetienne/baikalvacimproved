#include "math.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * Converting from ADM to BSSN quantities is required in the Einstein Toolkit,
 * as initial data are given in terms of ADM quantities, and Baikal evolves the BSSN quantities.
 */
void BSSN_to_ADM_BaikalVacuum(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  #pragma omp parallel for
  for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
    for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
      for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
        /*
         * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
         */
        const double hDD00 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD01 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD02 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD11 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD12 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double hDD22 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double aDD00 = aDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double aDD01 = aDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double aDD02 = aDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double aDD11 = aDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double aDD12 = aDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double aDD22 = aDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double vetU0 = vetU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double vetU1 = vetU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double vetU2 = vetU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double betU0 = betU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double betU1 = betU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double betU2 = betU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double trK = trKGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double cf = cfGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const double alpha = alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        /*
         * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
         */
        const double FDPart3_0 = (1.0/((cf)*(cf)));
        const double FDPart3_1 = FDPart3_0*(hDD00 + 1);
        const double FDPart3_4 = FDPart3_0*(hDD11 + 1);
        const double FDPart3_6 = FDPart3_0*(hDD22 + 1);
        const double FDPart3_7 = (1.0/3.0)*trK;
        gxx[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_1;
        gxy[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*hDD01;
        gxz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*hDD02;
        gyy[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_4;
        gyz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*hDD12;
        gzz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_6;
        kxx[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*aDD00 + FDPart3_1*FDPart3_7;
        kxy[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*FDPart3_7*hDD01 + FDPart3_0*aDD01;
        kxz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*FDPart3_7*hDD02 + FDPart3_0*aDD02;
        kyy[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*aDD11 + FDPart3_4*FDPart3_7;
        kyz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*FDPart3_7*hDD12 + FDPart3_0*aDD12;
        kzz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = FDPart3_0*aDD22 + FDPart3_6*FDPart3_7;
        alp[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = alpha;
        betax[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = vetU0;
        betay[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = vetU1;
        betaz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = vetU2;
        dtbetax[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = betU0;
        dtbetay[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = betU1;
        dtbetaz[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = betU2;
      } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
    } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
  } // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)
}
