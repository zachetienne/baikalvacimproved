#include "cctk.h"
#include "Slicing.h"
/*
 * Register slicing condition for NRPy+-generated thorn BaikalVacuum
 */
int RegisterSlicing_BaikalVacuum() {

  Einstein_RegisterSlicing ("BaikalVacuum");
  return 0;
}
