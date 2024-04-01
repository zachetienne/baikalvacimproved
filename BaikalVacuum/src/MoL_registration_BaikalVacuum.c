#include "stdio.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * Register evolved gridfunctions & RHSs
 * with the Method of Lines timestepper
 * MoL (the Einstein Toolkit Method of Lines thorn)
 * (MoL thorn, found in arrangements/CactusBase/MoL).
 * MoL documentation located in arrangements/CactusBase/MoL/doc
 */
void MoL_registration_BaikalVacuum(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  // Register evolution & RHS gridfunction groups with MoL, so it knows
  //  how to perform the appropriate timestepping

  group = CCTK_GroupIndex("BaikalVacuum::evol_variables");
  rhs = CCTK_GroupIndex("BaikalVacuum::evol_variables_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
}
