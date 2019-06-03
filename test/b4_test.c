#include <math.h>
#include <stdio.h>

#include "sl2cfoam.h"
#include "b4function.h"
#include "common.h"

int main() {

  //////////////////////// Spins //////////////////////

  dspin two_j1 = 2, two_j2 = 2, two_j3 = 2, two_j4 = 2;
  dspin two_l1 = 2, two_l2 = 2, two_l3 = 2, two_l4 = 2;

  //////////////////////// Check values and triangular inequalities//////////////////////

  if (two_l1 < two_j1 || two_l2 < two_j2 ||
      two_l3 < two_j3 || two_l4 < two_j4 ||
      two_j4 > two_j1+two_j2+two_j3      ||
      two_l4 > two_l1+two_l2+two_l3)
  {
      return 0;
  }

  //////////////////////// Immirzi parameter //////////////////////

  float Immirzi = 1.2;

  //////////////////////// Booster's parameters //////////////////////

  dspin two_k1 = two_j1, two_k2 = two_j2, two_k3 = two_j3, two_k4 = two_j4;
  float two_rho1 = Immirzi*(float)two_j1, two_rho2 = Immirzi*(float)two_j2;
  float two_rho3 = Immirzi*(float)two_j3, two_rho4 = Immirzi*(float)two_j4;

  dspin two_Ji1_min = (dspin)max(abs(two_j1-two_j2), abs(two_j3-two_j4));
  dspin two_Ji1_max = (dspin)min(two_j1+two_j2, two_j3+two_j4);
  dspin two_Ji1_bound = two_Ji1_max-two_Ji1_min;

  dspin two_Ji2_min = (dspin)max(abs(two_l1-two_l2),abs(two_l3-two_l4));
  dspin two_Ji2_max = (dspin)min(two_l1+two_l2, two_l3+two_l4);
  dspin two_Ji2_bound = two_Ji2_max-two_Ji2_min;

  //////////////////////// Initialize  the library //////////////////////

  sl2cfoam_init();

  //////////////////////// Compute Booster function //////////////////////

  // Double pointer to contain all possible intertwiner combinations
  long double  **B4_moy = B4Function( two_k1,  two_k2,  two_k3,  two_k4,
                                      two_rho1,  two_rho2,  two_rho3,  two_rho4,
                                      two_j1,  two_j2,  two_j3,  two_j4,
                                      two_l1,  two_l2,  two_l3,  two_l4
                                     );

  //////////////////////// Print Results //////////////////////

  printf ("\n");
  for (dspin two_Ji1=0; two_Ji1<=two_Ji1_bound; two_Ji1 += 2) {
      for (dspin two_Ji2=0; two_Ji2<=two_Ji2_bound; two_Ji2 += 2) {

          dspin Ji1 = two_Ji1/2;
          dspin Ji2 = two_Ji2/2;

          // NB To get the actual value for a couple of intertwiners
          // the following conversion is needed. From bound to actual intertwiner
          float Int1 = trunc(100*(float)(two_Ji1+two_Ji1_min)/2) / 100;
          float Int2 = trunc(100*(float)(two_Ji2+two_Ji2_min)/2) / 100;

          printf("B4: %.1f %.1f | %17Lg \n",Int1, Int2, B4_moy[Ji2][Ji1]);

      }
  }

  //////////////////////// Deallocate memory //////////////////////

  sl2cfoam_free();

}
