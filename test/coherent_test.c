#include <math.h>
#include <stdio.h>

#include "sl2cfoam.h"
#include "coherentstates.h"
#include "common.h"

int main(){

  //////////////////////// Spins //////////////////////

  int two_j1 = 2, two_j2 = 2, two_j3 = 2, two_j4 = 2;

  //////////////////////// Check values and triangular inequalities //////////////////////

  if (two_j4 > two_j1+two_j2+two_j3) {
      return 0;
  }

  dspin limi1 = max(abs(two_j1-two_j2), abs(two_j3-two_j4));
  dspin limi2 = min(two_j1+two_j2, two_j3+ two_j4);

  //////////////////////// Define Normals //////////////////////

  double f1, f2, f3, f4, t1, t2, t3, t4;

  f1 = 0;
  t1 = 0;
  f2 = 0;
  t2 = 1.9106332362490186;
  f3 = 2.0943951023931953;
  t3 = 1.9106332362490186;
  f4 = -2.0943951023931953;
  t4 = 1.9106332362490186;

  //////////////////////// Initialize the library //////////////////////

  sl2cfoam_init();

  /////////////////////////////////////////
  // Coherent State - Creation and printing
  /////////////////////////////////////////

  for (dspin two_i = limi1; two_i <= limi2; two_i += 2) {

    double complex coherentStatei = 0.0 + 0.0*I;

    coherentStatei = CoherentStateI(two_i, two_j1, two_j2, two_j3, two_j4,
                                    M_PI+(f1), M_PI-(t1),
                                    M_PI+(f2), M_PI-(t2),
                                    M_PI+(f3), M_PI-(t3),
                                    M_PI+(f4), M_PI-(t4),
                                    -1, -1, -1, -1);

    printf("%g | %17g %17g*I \n", (spin)(two_i)/2, creal(coherentStatei), cimag(coherentStatei));

  }


  //////////////////////// Deallocate memory //////////////////////

  sl2cfoam_free();

}
