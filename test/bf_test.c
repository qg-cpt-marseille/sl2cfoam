#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "sl2cfoam.h"
#include "error.h"

int main() {

	//////////////////// Initialize wigxjpf and hash table////////////////////

	struct sl2cfoam_config libconf;
	libconf.data_folder = "data_sl2cfoam/";
	libconf.store_ampls = 0;
	libconf.verbosity = 0;
	libconf.coupling = SL2CFOAM_COUPLING_IRREDUCIBLE;

	sl2cfoam_init_conf(&libconf);

	dspin two_j = 2;

	dspin two_i1 = 0; 
	dspin two_i2 = 0;
	dspin two_i3 = 0;
	dspin two_i4 = 0;
	dspin two_i5 = 0;


	//////////////////////////////////////////////////////
	// Test NON-PARALLEL version
	//////////////////////////////////////////////////////

	//////////////////// Start Computation ////////////////////

    printf("Testing 4-SIMPLEX BF NON-PARALLEL version...\n");

	double ampl = sl2cfoam_four_simplex_BF(two_j, two_j, two_j, two_j, two_j,
										   two_j, two_j, two_j, two_j, two_j,
										   two_i1, two_i2, two_i3, two_i4, two_i5);

	if (fabs(ampl) > 1.0) {
		error("Non-sense amplitude value obtained.");
	}

	printf("Amplitude: [all j = %d] %d %d %d %d %d | -> %17g \n\n",
				two_j, two_i1, two_i2, two_i3, two_i4, two_i5, ampl);


	//////////////////// Free Library ////////////////////

	sl2cfoam_free();

	return 0;
	
}
