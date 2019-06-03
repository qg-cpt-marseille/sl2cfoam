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

#include "common.h"
#include "hashing.h"
#include "khash.h"
#include "utilities.h"
#include "config.h"
#include "sl2cfoam.h"


double sl2cfoam_three_simplex(dspin two_j1, dspin two_j2, dspin two_j3, 
                              dspin two_j4, dspin two_j5, dspin two_j6,
                              dspin two_Dl, float Immirzi) {

    sl2cfoam_hash_three_ampl(two_j1, two_j1,
                             two_j2, two_j2,
                             two_j3, two_j3,
                             two_j4, two_j4,
                             two_j5, two_j5,
                             two_j6, two_j6,                                                                                                        
                             two_Dl, Immirzi);

    sl2cfoam_three_ampl_load(two_j1, two_j2, two_j3, two_Dl, Immirzi);

    double ampl;
    ampl = sl2cfoam_three_ampl_get(two_j1, two_j2, two_j3,
                                   two_j4, two_j5, two_j6,
                                   two_Dl, Immirzi);

    sl2cfoam_three_ampl_unload();

    return ampl;

}


// pointer to static table loaded into memory
// initialized to NULL
static kh_HashTableEPRL3_t* tableEPRL;

// properties of loaded table
// for safety check
static dspin two_js_loaded[3];
static dspin two_Dl_loaded;
static float Immirzi_loaded;

void sl2cfoam_three_ampl_load(dspin two_j1, dspin two_j2, dspin two_j3,
                              dspin two_Dl, float Immirzi) {

    if (tableEPRL != NULL) {
        error("one EPRL table is already loaded");
    }

    // setup data folder
    check_data_3simplex(Immirzi);

    struct data_folders* fd = get_data_folders(Immirzi);
    char filename[1024];
    char path_ampl[1024];

    // build path for amplitude values
    sprintf(filename, "%i.%i.%i_%i.eprl", two_j1, two_j2, two_j3, two_Dl);
    strcpy(path_ampl, fd->threesimp_imm_ampl);
    strcat(path_ampl, filename);


    // load hash table of amplitudes from disk
    if (!file_exist(path_ampl)) {
         error("table for given spins is not precomputed");
    }
    tableEPRL = kh_load(HashTableEPRL3, path_ampl);

    two_Dl_loaded = two_Dl;
    two_js_loaded[0] = two_j1;
    two_js_loaded[1] = two_j2;
    two_js_loaded[2] = two_j3;
    Immirzi_loaded = Immirzi;

}

void sl2cfoam_three_ampl_unload() {

    if (tableEPRL == NULL) {
         error("no previously loaded EPRL table");
    }

    kh_destroy(HashTableEPRL3, tableEPRL);
    tableEPRL = NULL;
    two_Dl_loaded = 1; // half-integer, invalid
    two_js_loaded[0] = 0;
    two_js_loaded[1] = 0;
    two_js_loaded[2] = 0;
    Immirzi_loaded = 0.0;

}

double sl2cfoam_three_ampl_get(dspin two_j1, dspin two_j2, dspin two_j3,
                               dspin two_j4, dspin two_j5, dspin two_j6,
                               dspin two_Dl, float Immirzi) {

    // check if table is loaded and correct
    if (tableEPRL == NULL) {
        error("EPRL table must be loaded before calling this function");
    }

    if (two_Dl_loaded != two_Dl) {
        error("loaded EPRL table has different Dl")
    }

    if (two_js_loaded[0] != two_j1 || two_js_loaded[1] != two_j2 || two_js_loaded[2] != two_j3) {
        error("loaded EPRL table has different boundary spins")
    }

    if (Immirzi_loaded != Immirzi) {
        error("loaded EPRL table has different Immirzi")
    }

    // if table loaded is empty return 0
    if (kh_size(tableEPRL) == 0) {
        return 0.0; 
    }

    // key for disk cache of current amplitude
    HashTableEPRL3_key_t keyEPRL = {two_j4, two_j5, two_j6};

    if (kh_get(HashTableEPRL3, tableEPRL, keyEPRL) != kh_end(tableEPRL)) {

        // value found!
        khint_t s = kh_get(HashTableEPRL3, tableEPRL, keyEPRL);
        double ampl = kh_val(tableEPRL, s);

        return ampl;

    }

    // value not found, error
    error("value not found in loaded EPRL table");

}