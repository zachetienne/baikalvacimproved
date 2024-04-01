#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
/*
 * Register symmetries for NRPy+-generated thorn BaikalVacuum
 */
void Symmetry_registration_oldCartGrid3D_BaikalVacuum(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Stores gridfunction parity across x=0, y=0, and z=0 planes, respectively
  int sym[3];

  // Next register parities for each gridfunction based on its name
  //    (to ensure this algorithm is robust, gridfunctions with integers
  //     in their base names are forbidden in NRPy+).

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[0] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::aDD00GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[1] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::aDD01GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::aDD02GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] *= -1;
  sym[1] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::aDD11GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::aDD12GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[2] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::aDD22GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  // (this gridfunction is a scalar -- no need to change default sym[]'s!)
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::alphaGF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::betU0GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::betU1GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::betU2GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  // (this gridfunction is a scalar -- no need to change default sym[]'s!)
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::cfGF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[0] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::hDD00GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[1] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::hDD01GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::hDD02GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] *= -1;
  sym[1] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::hDD11GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::hDD12GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[2] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::hDD22GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::lambdaU0GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::lambdaU1GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::lambdaU2GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  // (this gridfunction is a scalar -- no need to change default sym[]'s!)
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::trKGF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::vetU0GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::vetU1GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::vetU2GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[0] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::RbarDD00GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[1] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::RbarDD01GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::RbarDD02GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] *= -1;
  sym[1] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::RbarDD11GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::RbarDD12GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[2] *= -1;
  sym[2] *= -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::RbarDD22GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  // (this gridfunction is a scalar -- no need to change default sym[]'s!)
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::HGF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  // (this gridfunction is a scalar -- no need to change default sym[]'s!)
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::MSQUAREDGF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[0] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::MU0GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[1] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::MU1GF");

  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed
  //    to account for gridfunction parity across
  //    x=0, y=0, and/or z=0 planes, respectively
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "BaikalVacuum::MU2GF");
}
