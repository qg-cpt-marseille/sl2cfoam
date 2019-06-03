#include <math.h>
#include <stdio.h>

#include "sl2cfoam.h"
#include "b3function.h"

int main() {

  dspin two_j1, two_j2, two_j3;
  dspin two_l1, two_l2, two_l3;

  two_j1 = 2; two_j2 = 2; two_j3 = 2;
  two_l1 = 2; two_l2 = 2; two_l3 = 2;

  double Immirzi = 1.2;

  sl2cfoam_init();

  double Boost3 = B3(two_j1, two_j2, two_j3,
                     two_l1, two_l2, two_l3,
                     Immirzi);

  printf("B3: %i %i %i %i %i %i | %17g \n",
         two_j1, two_j2, two_j3, two_l1, two_l2, two_l3,
         Boost3);

  sl2cfoam_free();

  return 0;
}
