#include "math.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * 
 * Set up NewRad boundary conditions.
 *    As explained in lean_public/LeanBSSNMoL/src/calc_bssn_rhs.F90,
 *    the function NewRad_Apply takes the following arguments:
 *    NewRad_Apply(cctkGH, var, rhs, var0, v0, radpower),
 *      which implement the boundary condition:
 *        var  =  var_at_infinite_r + u(r-var_char_speed*t)/r^var_radpower
 *   Obviously for var_radpower>0, var_at_infinite_r is the value of
 *     the variable at r->infinity. var_char_speed is the propagation
 *     speed at the outer boundary, and var_radpower is the radial
 *     falloff rate.
 */
void specify_NewRad_BoundaryConditions_parameters_BaikalVacuum(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  NewRad_Apply(cctkGH, aDD00GF, aDD00_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD01GF, aDD01_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD02GF, aDD02_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD11GF, aDD11_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD12GF, aDD12_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD22GF, aDD22_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, alphaGF, alpha_rhsGF, 1.0, 1.41421356237310, 1.0);
  NewRad_Apply(cctkGH, betU0GF, betU0_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, betU1GF, betU1_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, betU2GF, betU2_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, cfGF, cf_rhsGF, 1.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD00GF, hDD00_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD01GF, hDD01_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD02GF, hDD02_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD11GF, hDD11_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD12GF, hDD12_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD22GF, hDD22_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, lambdaU0GF, lambdaU0_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, lambdaU1GF, lambdaU1_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, lambdaU2GF, lambdaU2_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, trKGF, trK_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, vetU0GF, vetU0_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, vetU1GF, vetU1_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, vetU2GF, vetU2_rhsGF, 0.0, 1.0, 1.0);
}
