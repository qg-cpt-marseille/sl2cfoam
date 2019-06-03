#include <math.h>
#include <stdio.h>

#include "sl2cfoam.h"
#include "recouplingSL2C.h"

int main() {

  dspin two_j1, two_j2, two_j3, two_j4, two_j5, two_j6;
  dspin two_k1, two_k2, two_k3, two_k4, two_k5, two_k6;
  float two_rho1, two_rho2, two_rho3, two_rho4, two_rho5, two_rho6;

  dspin two_Dl;
  double Immirzi = 1.2;

  two_j1 = 2; two_j2 = 2; two_j3 = 2;
  two_j4 = 2; two_j5 = 2; two_j6 = 2;

  two_k1 = 2; two_k2 = 2; two_k3 = 2;
  two_k4 = 2; two_k5 = 2; two_k6 = 2;

  two_rho1 = 2*Immirzi; two_rho2 = 2*Immirzi; two_rho3 = 2*Immirzi;
  two_rho4 = 2*Immirzi; two_rho5 = 2*Immirzi; two_rho6 = 2*Immirzi;

  sl2cfoam_init();

  double value = wigner_3j_sl2c(two_j1, two_j2, two_j3,
                                two_k1, two_k2, two_k3,
                                two_rho1, two_rho2, two_rho3);

  printf("\nSL(2,C) {3j}:  %17g \n\n", value);

  dspin two_Dj = two_k1;

  double value1 = wigner_6j_sl2c(two_k1, two_k2, two_k3,
                                 two_k4, two_k5, two_k6,
                                 two_rho1, two_rho2, two_rho3,
                                 two_rho4, two_rho5, two_rho6,
                                 two_Dj);

  printf("SL(2,C) {6j}: %17g for two_Dj = %i\n", value1 , two_Dj);

  sl2cfoam_free();

  return 0;
}
