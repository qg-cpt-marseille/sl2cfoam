#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "sl2cfoam.h"
#include "error.h"


int main() {

  //////////////////// Boundary Spins and Dl ////////////////////

  dspin two_j1, two_j2, two_j3, two_j4, two_j5, two_j6;

  two_j1 = 4; two_j2 = 4;
  two_j3 = 4; two_j4 = 4;
  two_j5 = 4; two_j6 = 4;

  dspin two_Dl = 10;

  //////////////////// Immirzi Definition ////////////////////

  float Immirzi = 1.2;

	//////////////////// Initialize wigxjpf and hash table////////////////////

	struct sl2cfoam_config libconf;
	libconf.data_folder = "data_sl2cfoam/";
	libconf.store_ampls = 0;
	libconf.verbosity = 0;
	libconf.coupling = SL2CFOAM_COUPLING_REDUCIBLE;

	sl2cfoam_init_conf(&libconf);

  //////////////////// Start Computation ////////////////////

  double ampl = sl2cfoam_three_simplex( two_j1, two_j2, two_j3,
                                        two_j4, two_j5, two_j6,
                                        two_Dl, Immirzi);
  if (fabs(ampl) > 1.0) {
    warning("Non-sense amplitude value obtained.");
  }

  printf("EPRL 3-Simplex: %i %i %i %i %i %i | dl = %i | -> %17g \n",
        two_j1, two_j2, two_j3, two_j4, two_j5, two_j6, two_Dl, ampl);

  //////////////////// Free Library ////////////////////

  sl2cfoam_free();

  return 0;
}
