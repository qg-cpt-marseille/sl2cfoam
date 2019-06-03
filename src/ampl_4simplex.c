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


double sl2cfoam_four_simplex(dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4,
                             dspin two_j5, dspin two_j6, dspin two_j7, dspin two_j8,
                             dspin two_j9, dspin two_j10,
                             dspin two_i1, dspin two_i2, dspin two_i3,
                             dspin two_i4, dspin two_i5,
                             dspin two_Dl, float Immirzi) {

    dspin two_js[10] = {two_j1, two_j2, two_j3, two_j4, two_j5,
                        two_j6, two_j7, two_j8, two_j9, two_j10};

    sl2cfoam_hash_symbols(two_js, two_i4, two_i4,
                          two_Dl, Immirzi);

    sl2cfoam_hash_four_ampl(two_js,
                            two_i1, two_i1,
                            two_i2, two_i2,
                            two_i3, two_i3,
                            two_i4, two_i4,
                            two_i5, two_i5,
                            two_Dl, Immirzi);

    sl2cfoam_four_ampl_load(two_js, two_i4, two_Dl, Immirzi);

    double ampl;
    ampl = sl2cfoam_four_ampl_get(two_js,
                                  two_i1, two_i2, two_i3, two_i4, two_i5,
                                  two_Dl, Immirzi);

    sl2cfoam_four_ampl_unload();

    return ampl;

}

double sl2cfoam_four_simplex_BF(dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4,
                                dspin two_j5, dspin two_j6, dspin two_j7, dspin two_j8,
                                dspin two_j9, dspin two_j10,
                                dspin two_i1, dspin two_i2, dspin two_i3,
                                dspin two_i4, dspin two_i5) {

    dspin two_js[10] = {two_j1, two_j2, two_j3, two_j4, two_j5,
                        two_j6, two_j7, two_j8, two_j9, two_j10};

    sl2cfoam_hash_symbols_BF(two_js, two_i4, two_i4);

    sl2cfoam_hash_four_ampl_BF(two_js,
                               two_i1, two_i1,
                               two_i2, two_i2,
                               two_i3, two_i3,
                               two_i4, two_i4,
                               two_i5, two_i5);

    sl2cfoam_four_ampl_load_BF(two_js, two_i4);

    double ampl;
    ampl = sl2cfoam_four_ampl_get_BF(two_js,
                                     two_i1, two_i2, two_i3, two_i4, two_i5);

    sl2cfoam_four_ampl_unload();

    return ampl;

}

// pointer to static table loaded into memory
// initialized to NULL
static kh_HashTableEPRL_t* tableEPRL;

// properties of loaded table
// for safety check
static dspin two_i4_loaded;
static dspin two_js_loaded[4];
static dspin two_Dl_loaded;
static float Immirzi_loaded;

void sl2cfoam_four_ampl_load(dspin two_js[10], dspin two_i4,
                             dspin two_Dl, float Immirzi) {

    if (tableEPRL != NULL) {
        error("one EPRL table is already loaded");
    }

    dspin two_j5, two_j6, two_j9, two_j10;

    two_j5 = two_js[4];
    two_j6 = two_js[5];
    two_j9 = two_js[8];
    two_j10 = two_js[9];

    // setup data folder
    check_data_4simplex(Immirzi);

    struct data_folders* fd = get_data_folders(Immirzi);
    char filename[1024];
    char path_ampl[1024];

    // build path for amplitude values
    sprintf(filename, "%i.%i.%i.%i_%i_%i.eprl", two_j5, two_j6, two_j9, two_j10, two_i4, two_Dl);
    strcpy(path_ampl, fd->foursimp_imm_hashtableampl_ampl);
    strcat(path_ampl, filename);

    // load hash table of amplitudes from disk
    if (!file_exist(path_ampl)) {
         error("table for given i4 not precomputed");
    }
    tableEPRL = kh_load(HashTableEPRL, path_ampl);

    two_i4_loaded = two_i4;
    two_Dl_loaded = two_Dl;
    two_js_loaded[0] = two_j5;
    two_js_loaded[1] = two_j6;
    two_js_loaded[2] = two_j9;
    two_js_loaded[3] = two_j10;
    Immirzi_loaded = Immirzi;

}

void sl2cfoam_four_ampl_load_BF(dspin two_js[10], dspin two_i4) {

    if (tableEPRL != NULL) {
        error("one EPRL table is already loaded");
    }

    dspin two_j5, two_j6, two_j9, two_j10;

    two_j5 = two_js[4];
    two_j6 = two_js[5];
    two_j9 = two_js[8];
    two_j10 = two_js[9];

    // setup data folder
    check_data_4simplex(0);

    struct data_folders* fd = get_data_folders(0);
    char filename[1024];
    char path_ampl[1024];

    // build path for amplitude values
    sprintf(filename, "%i.%i.%i.%i_%i.bf", two_j5, two_j6, two_j9, two_j10, two_i4);
    strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
    strcat(path_ampl, filename);

    // load hash table of amplitudes from disk
    if (!file_exist(path_ampl)) {
         error("table for given i4 not precomputed");
    }
    tableEPRL = kh_load(HashTableEPRL, path_ampl);

    two_i4_loaded = two_i4;
    two_Dl_loaded = 0;
    two_js_loaded[0] = two_j5;
    two_js_loaded[1] = two_j6;
    two_js_loaded[2] = two_j9;
    two_js_loaded[3] = two_j10;
    Immirzi_loaded = 0.0;

}

void sl2cfoam_four_ampl_unload() {

    if (tableEPRL == NULL) {
         error("no previously loaded EPRL table");
    }

    kh_destroy(HashTableEPRL, tableEPRL);
    tableEPRL = NULL;
    two_i4_loaded = 0;
    two_Dl_loaded = 1; // half-integer, invalid
    two_js_loaded[0] = 0;
    two_js_loaded[1] = 0;
    two_js_loaded[2] = 0;
    two_js_loaded[3] = 0;
    Immirzi_loaded = 0.0;

}

double sl2cfoam_four_ampl_get(dspin two_js[10],
                              dspin two_i1, dspin two_i2, dspin two_i3,
                              dspin two_i4, dspin two_i5,
                              dspin two_Dl, float Immirzi) {

    dspin two_j1, two_j2, two_j3, two_j4, two_j5,
          two_j6, two_j7, two_j8, two_j9, two_j10;

    two_j1 = two_js[0];
    two_j2 = two_js[1];
    two_j3 = two_js[2];
    two_j4 = two_js[3];
    two_j5 = two_js[4];
    two_j6 = two_js[5];
    two_j7 = two_js[6];
    two_j8 = two_js[7];
    two_j9 = two_js[8];
    two_j10 = two_js[9];

    // check if table is loaded and correct
    if (tableEPRL == NULL) {
        error("EPRL table must be loaded before calling this function");
    }

    if (two_i4_loaded != two_i4) {
        error("loaded EPRL table has different i4")
    }

    if (two_Dl_loaded != two_Dl) {
        error("loaded EPRL table has different Dl")
    }

    if (two_js_loaded[0] != two_j5 || two_js_loaded[1] != two_j6 ||
        two_js_loaded[2] != two_j9 || two_js_loaded[1] != two_j10) {
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
    HashTableEPRL_key_t keyEPRL = {two_j1, two_j2, two_j3, two_j4, two_j7, two_j8,
                                   two_i1, two_i2, two_i3, two_i5};

    if (kh_get(HashTableEPRL, tableEPRL, keyEPRL) != kh_end(tableEPRL)) {

        // value found!
        khint_t s = kh_get(HashTableEPRL, tableEPRL, keyEPRL);
        double ampl = kh_val(tableEPRL, s);

        return ampl;

    }

    // value not found, error
    error("value not found in loaded EPRL table");

}

double sl2cfoam_four_ampl_get_BF(dspin two_js[10],
                                 dspin two_i1, dspin two_i2, dspin two_i3,
                                 dspin two_i4, dspin two_i5) {

    return sl2cfoam_four_ampl_get(two_js,
                                  two_i1, two_i2, two_i3, two_i4, two_i5,
                                  0, 0);

}

