/* Copyright 2019 Giorgio Sarno, Pietro Donà and Francesco Gozzini */

/* sl2cfoam is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   sl2cfoam is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.*/

///////////////////////////////////////////////////////////////
// Initialization functions.
///////////////////////////////////////////////////////////////

#include <mpfr.h>

#include "cgamma.h"
#include "jsymbols.h"
#include "sl2cfoam.h"
#include "config.h"

void sl2cfoam_init_conf(struct sl2cfoam_config* conf) {

    set_data_root(conf->data_folder);
    set_coupling(conf->coupling);
    set_verbosity(conf->verbosity);

	init_wigxjpf_global();
	init_complex_gamma();

    // set default precision for MPFR
    // gamma functions precision is set in config.h
    // with the MPBITS value
    mpfr_set_default_prec(128);
    mpfr_set_default_rounding_mode(MPFR_RNDN);

}

void sl2cfoam_init() {
	
    struct sl2cfoam_config def;
    def.data_folder = "../data_sl2cfoam";
    def.coupling = SL2CFOAM_COUPLING_REDUCIBLE;
    def.verbosity = SL2CFOAM_VERBOSE_OFF;

    sl2cfoam_init_conf(&def);

}

void sl2cfoam_init_thread() {
	init_wigxjpf_thread();
}

void sl2cfoam_free() {
	clear_wigxjpf_global();
	clear_complex_gamma();
};

void sl2cfoam_free_thread() {
	clear_wigxjpf_thread();
};