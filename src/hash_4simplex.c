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

#include <omp.h>

#include "common.h"
#include "hashing.h"
#include "khash.h"
#include "config.h"
#include "utilities.h"
#include "jsymbols.h"
#include "b4function.h"
#include "sl2cfoam.h"


void sl2cfoam_hash_symbols(dspin two_js[10],
                           dspin two_i4_min, dspin two_i4_max,
                           dspin two_Dl, float Immirzi) {

    check_data_4simplex(Immirzi);

    for (dspin two_d = 0; two_d <= two_Dl; two_d += 2) {

        B4_Hash(two_js[0], two_js[1], two_js[2], two_js[3], two_js[4],
                two_js[5], two_js[6], two_js[7], two_js[8], two_js[9],
                two_d, Immirzi);

    }

    for (dspin two_i4 = two_i4_min; two_i4 <= two_i4_max; two_i4 += 2) {

        for (dspin two_d = 0; two_d <= two_Dl; two_d += 2) {

            J15Symbol_Hash(two_js[0], two_js[1], two_js[2], two_js[3], two_js[4],
                           two_js[5], two_js[6], two_js[7], two_js[8], two_js[9],
                           two_i4, two_d);

        }

    }

}

void sl2cfoam_hash_symbols_BF(dspin two_js[10],
                              dspin two_i4_min, dspin two_i4_max) {

    check_data_4simplex(0);

    for (dspin two_i4 = two_i4_min; two_i4 <= two_i4_max; two_i4 += 2) {

        J15Symbol_Hash(two_js[0], two_js[1], two_js[2], two_js[3], two_js[4],
                        two_js[5], two_js[6], two_js[7], two_js[8], two_js[9],
                        two_i4, 0);

    }

}

void sl2cfoam_hash_four_ampl(dspin two_js[10],
                             dspin two_i1_min, dspin two_i1_max, 
                             dspin two_i2_min, dspin two_i2_max, 
                             dspin two_i3_min, dspin two_i3_max, 
                             dspin two_i4_min, dspin two_i4_max, 
                             dspin two_i5_min, dspin two_i5_max, 
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

    // setup data folder
    check_data_4simplex(Immirzi);
    struct data_folders* fd = get_data_folders(Immirzi);

    //////////////////////////////////////////////////////////////////////////////
    // loop over i4 (safe to be parallelized)
    //////////////////////////////////////////////////////////////////////////////

    // optimize OMP for multi-shells computations
    if (two_Dl > 0) {
        omp_set_nested(1);      // enable nested parallelization
    } else {
        omp_set_nested(0);      // disable nested parallelization
    }
    
    #pragma omp parallel
    {

    dspin two_i4;

    // per-thread initialization
    sl2cfoam_init_thread();
        
    #pragma omp for
    for (two_i4 = two_i4_min; two_i4 <= two_i4_max; two_i4 += 2) {

        char filename[1024];
        char path_ampl[1024];
        char path_6j[1024];
        char path_boost[1024];

        // build path for amplitude values
        sprintf(filename, "%i.%i.%i.%i_%i_%i.eprl", two_j5, two_j6, two_j9, two_j10, two_i4, two_Dl);
        strcpy(path_ampl, fd->foursimp_imm_hashtableampl_ampl);
        strcat(path_ampl, filename);

        // build path for 6js
        sprintf(filename, "%i.%i.%i.%i_%i_%i.6j", two_j5, two_j6, two_j9, two_j10, two_i4, two_Dl);
        strcpy(path_6j, fd->foursimp_hashtable6j_6j);
        strcat(path_6j, filename);

        // build path for boost
        sprintf(filename, "%i.%i.%i.%i_%i.boost", two_j5, two_j6, two_j9, two_j10, two_Dl);
        strcpy(path_boost, fd->foursimp_imm_boost);
        strcat(path_boost, filename);

        // load hash tables for BOOSTERS and 6Js
        khash_t(HashTableJ6) *h = kh_load(HashTableJ6, path_6j);
        khash_t(HashTableBooster) *h2 = kh_load(HashTableBooster, path_boost);

        // initialize or load hash table of amplitudes on disk
        kh_HashTableEPRL_t *h3 = NULL;
        if (file_exist(path_ampl)) {

            h3 = kh_load(HashTableEPRL, path_ampl);

        } else {

            // table not found, initialize it
            h3 = kh_init(HashTableEPRL);
            
        }

        // check 6js table: if all 15js zero for this i4
        // write empty EPRL table and continue;
        if (kh_size(h) == 0) {
            
            // write EMPTY hash table for this i4 to disk
            kh_write(HashTableEPRL, h3, path_ampl);

            // free memory
            kh_destroy(HashTableJ6, h);
            kh_destroy(HashTableBooster, h2);
            kh_destroy(HashTableEPRL, h3);

            continue;
            
        }


        dspin two_i1, two_i2, two_i3, two_i5;

        int ret; // return codes from khash

        for (two_i5 = two_i5_min; two_i5 <= two_i5_max; two_i5 += 2) {
        for (two_i3 = two_i3_min; two_i3 <= two_i3_max; two_i3 += 2) {
        for (two_i2 = two_i2_min; two_i2 <= two_i2_max; two_i2 += 2) {
        for (two_i1 = two_i1_min; two_i1 <= two_i1_max; two_i1 += 2) {


            // key for disk cache of current amplitude
            HashTableEPRL_key_t keyEPRL = {two_j1, two_j2, two_j3, two_j4, two_j7, two_j8,
                                           two_i1, two_i2, two_i3, two_i5};

            // check if value has been already computed and stored
            if (kh_get(HashTableEPRL, h3, keyEPRL) != kh_end(h3)) {

                // value found!
                // go to next intertwiner set
                continue;

            }

            //////////////////////////////////////////////////////////////////////
            // value was not found or table was initialized
            // compute it then
            //////////////////////////////////////////////////////////////////////

            if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

                // with reducible coupling amplitude there are
                // additional triangular inequalities on i1

                dspin two_i1_min_red = (dspin) max(abs(two_i2-two_i3), abs(two_i5-two_i4));
                dspin two_i1_max_red = (dspin) min(two_i2+two_i3, two_i5+two_i4);

                if (!(max(two_i1_min, two_i1_min_red) <= two_i1 && two_i1 <= min(two_i1_max, two_i1_max_red))) {

                    // save 0 and pass to next

                    khint_t s = kh_put(HashTableEPRL, h3, keyEPRL, &ret);

                    if (ret == -1) {
                        error("error inserting key into EPRL hash table");
                    }

                    if (ret == 1) {
                        kh_val(h3, s) = 0.0;
                    }

                    continue;

                }

            } 

            double ampl = 0.0;
            double err = 0.0;
            double cutoff = 1e-40;

            //////////////////////////////////////////////////////////////////////
            // start parallel summations over l
            //////////////////////////////////////////////////////////////////////

            #pragma omp parallel reduction (+:ampl, err)
            {

            init_wigxjpf_thread();  

            #pragma omp for collapse(6)
            for (dspin two_l1 = two_j1; two_l1 <= two_j1 + two_Dl; two_l1 += 2) {
            for (dspin two_l2 = two_j2; two_l2 <= two_j2 + two_Dl; two_l2 += 2) {
            for (dspin two_l3 = two_j3; two_l3 <= two_j3 + two_Dl; two_l3 += 2) {
            for (dspin two_l4 = two_j4; two_l4 <= two_j4 + two_Dl; two_l4 += 2) {
            for (dspin two_l7 = two_j7; two_l7 <= two_j7 + two_Dl; two_l7 += 2) {
            for (dspin two_l8 = two_j8; two_l8 <= two_j8 + two_Dl; two_l8 += 2) {

                //////////////////// Range for internal K intertwiners ////////////////////

                dspin two_limk1_1, two_limk1_2, two_limk2_1, two_limk2_2, 
                        two_limk3_1, two_limk3_2, two_limk5_1, two_limk5_2;


                if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

                    two_limk2_1 = max(abs(two_j10-two_l7), abs(two_l1-two_l2));
                    two_limk2_2 = min(two_j10+two_l7, two_l1+two_l2);
                    two_limk3_1 = max(abs(two_j9-two_l8), abs(two_l3-two_l1));
                    two_limk3_2 = min(two_j9+two_l8, two_l3+two_l1);
                    two_limk5_1 = max(abs(two_l4-two_j6), abs(two_l7-two_l8));
                    two_limk5_2 = min(two_l4+two_j6, two_l7+two_l8);

                } else {

                    two_limk1_1 = max(abs(two_l2-two_l3), abs(two_j5-two_l4));
                    two_limk1_2 = min(two_l2+two_l3, two_j5+two_l4);
                    two_limk2_1 = max(abs(two_l1-two_j10), abs(two_l7-two_l2));
                    two_limk2_2 = min(two_l1+two_j10, two_l7+two_l2);
                    two_limk3_1 = max(abs(two_j9-two_l8), abs(two_l3-two_l1));
                    two_limk3_2 = min(two_j9+two_l8, two_l3+two_l1);
                    two_limk5_1 = max(abs(two_l4-two_l7), abs(two_l8-two_j6));
                    two_limk5_2 = min(two_l4+two_l7, two_l8+two_j6);

                }

                ///////////////// Initialize pointers for 6j keys and boosters /////////////////

                dspin *A = malloc(6 * sizeof(dspin));

                khint_t hintBooster;
                double boosterB, boosterC, boosterE, boosterA;

                //////////////////// Start K summations ////////////////////

                //////////////////// K2 for loop ////////////////////

                for (dspin  two_k2 = two_limk2_1; two_k2 <= two_limk2_2; two_k2 += 2) {

                    HashTableBooster_key_t keyBoosterB;

                    if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

                        keyBoosterB = (HashTableBooster_key_t) {two_j10, two_j7, two_j1, two_j2,
                                                                two_j10, two_l7, two_l1, two_l2,
                                                                two_i2, two_k2};

                    } else {

                        keyBoosterB = (HashTableBooster_key_t) {two_j1, two_j10, two_j7, two_j2,
                                                                two_l1, two_j10, two_l7, two_l2,
                                                                two_i2, two_k2};
                    }

                    hintBooster = kh_get(HashTableBooster, h2, keyBoosterB);
                    boosterB = d(two_k2) * kh_val(h2, hintBooster);

                    // Apply K cutoffs
                    if (fabs(boosterB) < cutoff) continue;

                //////////////////// K3 for loop ////////////////////

                for (dspin two_k3 = two_limk3_1; two_k3 <= two_limk3_2; two_k3 += 2) {

                    HashTableBooster_key_t keyBoosterC = {two_j9, two_j8, two_j3, two_j1,
                                                          two_j9, two_l8, two_l3, two_l1,
                                                          two_i3, two_k3};

                    hintBooster = kh_get(HashTableBooster, h2, keyBoosterC);
                    boosterC = d(two_k3) * kh_val(h2, hintBooster);

                    if (fabs(boosterC) < cutoff) continue;

                //////////////////// K5 for loop ////////////////////

                for (dspin two_k5 = two_limk5_1; two_k5 <= two_limk5_2; two_k5 += 2) {

                    HashTableBooster_key_t keyBoosterE;

                    if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

                        keyBoosterE = (HashTableBooster_key_t) {two_j4, two_j6, two_j7, two_j8,
                                                                two_l4, two_j6, two_l7, two_l8,
                                                                two_i5, two_k5};

                    } else {

                        keyBoosterE = (HashTableBooster_key_t) {two_j4, two_j7, two_j8, two_j6,
                                                                two_l4, two_l7, two_l8, two_j6,
                                                                two_i5, two_k5};

                    }

                    hintBooster = kh_get(HashTableBooster, h2, keyBoosterE);
                    boosterE = d(two_k5) * kh_val(h2, hintBooster);

                    if ( fabs(boosterE) < cutoff ) continue;

                    if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {
                    
                        two_limk1_1 = max(max(abs(two_l3-two_l2),abs(two_j5-two_l4)), max(abs(two_k3-two_k2),abs(two_k5-two_i4)));
                        two_limk1_2 = min(min(two_l3+two_l2,two_j5+two_l4), min(two_k2+two_k3,two_k5+two_i4));

                    }

                //////////////////// K1 for loop ////////////////////

                for (dspin two_k1 = two_limk1_1; two_k1 <= two_limk1_2; two_k1 += 2) {

                    HashTableBooster_key_t keyBoosterA = {two_j2, two_j3, two_j5, two_j4,
                                                          two_l2, two_l3, two_j5, two_l4,
                                                          two_i1, two_k1};

                    hintBooster = kh_get(HashTableBooster, h2, keyBoosterA);
                    boosterA = d(two_k1) * kh_val(h2, hintBooster);

                    if (fabs(boosterA) < cutoff) continue;

                    // Convert Boosters to MPFR variables
                    mpfr_t  boosterMPFR_A, boosterMPFR_B, boosterMPFR_C, boosterMPFR_E, boostersMPFR;

                    mpfr_init_set_d(boosterMPFR_A, boosterA, MPFR_RNDN);
                    mpfr_init_set_d(boosterMPFR_B, boosterB, MPFR_RNDN);
                    mpfr_init_set_d(boosterMPFR_C, boosterC, MPFR_RNDN);
                    mpfr_init_set_d(boosterMPFR_E, boosterE, MPFR_RNDN);
                    mpfr_init(boostersMPFR);

                    // Multiply Boosters as MPFR variables
                    mpfr_mul(boostersMPFR, boosterMPFR_A, boosterMPFR_B, MPFR_RNDN);
                    mpfr_mul(boostersMPFR, boostersMPFR, boosterMPFR_C, MPFR_RNDN);
                    mpfr_mul(boostersMPFR, boostersMPFR, boosterMPFR_E, MPFR_RNDN);

                    mpfr_clears(boosterMPFR_A, boosterMPFR_B, boosterMPFR_C, boosterMPFR_E, NULL);

                    // Compute 15j Symbol
                    double val15j = J15Symbol(h, &A,
                                              two_l1, two_l2, two_l3, two_l4, two_j5,
                                              two_j6, two_l7, two_l8, two_j9, two_j10,
                                              two_k1, two_k2, two_k3, two_i4, two_k5);

                    // Convert 15j symbol and coherent states to MPFR variables
                    mpfr_t val15jMPFR;
                    mpfr_init_set_d(val15jMPFR, val15j, MPFR_RNDN);

                    // Compute SL2C value as an MPFR variable
                    mpfr_t amplMPFR;
                    mpfr_init(amplMPFR);
                    mpfr_mul(amplMPFR, val15jMPFR, boostersMPFR, MPFR_RNDN);

                    // Back to double and Kahan compensated summation for I, K
                    compsum_mpfr(&err, &ampl, amplMPFR);

                    mpfr_clears(boostersMPFR, val15jMPFR, amplMPFR, NULL);

                } // k2
                } // k3
                } // k5
                } // k1

                free(A);

            } // l8
            } // l7
            } // l4
            } // l3
            } // l2
            } // l1

            clear_wigxjpf_thread();

            } // omp parallel on l

            // normalize amplitude with dimensions of all 5 boundary intertwiners
            ampl = ampl * d(two_i1) * d(two_i2) * d(two_i3) * d(two_i4) * d(two_i5);

            // check if amplitude has consistent value
            if (fabs(ampl) > 1.0) {
                error("amplitude has non-consistent value");
            }


            ////////////////////////////////////////////////////////
            // put key into amplitudes cache to be written to disk
            ////////////////////////////////////////////////////////

            khint_t s = kh_put(HashTableEPRL, h3, keyEPRL, &ret);

            if (ret == -1) {
                error("error inserting key into EPRL hash table");
            }

            if (ret == 1) {

                kh_val(h3, s) = ampl;

                // check if insertion worked
                s = kh_get(HashTableEPRL, h3, keyEPRL);
                if (kh_val(h3, s) == 0 && ampl != 0) {

                    // retry
                    kh_val(h3, s) = ampl;

                    // check again and now if not fail
                    if (kh_val(h3, s) == 0 && ampl != 0) {
                        error("error inserting value into EPRL hashtable");
                    }

                }

            }

            //////////////////////////////////////////////////////////////
            // IF the coupling is IRREDUCIBLE
            // AND IF all the boundary spins are equal
            // store also the companion amplitude obtained by reflecting
            // the pentachoron on the axis crossing i4
            //////////////////////////////////////////////////////////////

            if (get_coupling() == SL2CFOAM_COUPLING_IRREDUCIBLE &&
                two_j1 == two_j2 && two_j2 == two_j3 && two_j3 == two_j4 && two_j4 == two_j5 &&
                two_j5 == two_j6 && two_j6 == two_j7 && two_j7 == two_j8 && two_j8 == two_j9 && two_j9 == two_j10) {

                // relected key i1 <-> i2 , i3 <-> i5
                HashTableEPRL_key_t keyEPRLref = {two_j1, two_j2, two_j3, two_j4, two_j7, two_j8,
                                                  two_i2, two_i1, two_i5, two_i3};

                s = kh_put(HashTableEPRL, h3, keyEPRLref, &ret);

                if (ret == -1) {
                    error("error inserting key into EPRL hash table");
                }

                if (ret == 1) {

                    kh_val(h3, s) = ampl;

                    // check if insertion worked
                    s = kh_get(HashTableEPRL, h3, keyEPRLref);
                    if (kh_val(h3, s) == 0 && ampl != 0) {

                        // retry
                        kh_val(h3, s) = ampl;

                        // check again and now if not fail
                        if (kh_val(h3, s) == 0 && ampl != 0) {
                            error("error inserting value into EPRL hashtable");
                        }

                    }

                }

            }

        } // i1
        } // i2
        } // i3
        } // i5

        // write hash table for this i4 to disk
        kh_write(HashTableEPRL, h3, path_ampl);

        // free memory
        kh_destroy(HashTableJ6, h);
        kh_destroy(HashTableBooster, h2);
        kh_destroy(HashTableEPRL, h3);

    } // i4

    sl2cfoam_free_thread();

    } // omp parallel 

}


void sl2cfoam_hash_four_ampl_BF(dspin two_js[10],
                                dspin two_i1_min, dspin two_i1_max, 
                                dspin two_i2_min, dspin two_i2_max, 
                                dspin two_i3_min, dspin two_i3_max, 
                                dspin two_i4_min, dspin two_i4_max, 
                                dspin two_i5_min, dspin two_i5_max) {                      

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

    // setup data folder
    check_data_4simplex(0);
    struct data_folders* fd = get_data_folders(0);

    //////////////////////////////////////////////////////////////////////////////
    // loop over i4 (safe to be parallelized)
    //////////////////////////////////////////////////////////////////////////////

    #pragma omp parallel
    {

    dspin two_i4;

    // per-thread initialization
    sl2cfoam_init_thread();
        
    #pragma omp for
    for (two_i4 = two_i4_min; two_i4 <= two_i4_max; two_i4 += 2) {

        char filename[1024];
        char path_ampl[1024];
        char path_6j[1024];

        // build path for amplitude values
        sprintf(filename, "%i.%i.%i.%i_%i.bf", two_j5, two_j6, two_j9, two_j10, two_i4);
        strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
        strcat(path_ampl, filename);

        // build path for 6js
        sprintf(filename, "%i.%i.%i.%i_%i_%i.6j", two_j5, two_j6, two_j9, two_j10, two_i4, 0);
        strcpy(path_6j, fd->foursimp_hashtable6j_6j);
        strcat(path_6j, filename);

        // load hash tables for 6Js
        khash_t(HashTableJ6) *h = kh_load(HashTableJ6, path_6j);

        // initialize or load hash table of amplitudes on disk
        kh_HashTableEPRL_t *h3 = NULL;
        if (file_exist(path_ampl)) {

            h3 = kh_load(HashTableEPRL, path_ampl);

        } else {

            // table not found, initialize it
            h3 = kh_init(HashTableEPRL);
            
        }

        // check 6js table: if all 15js zero for this i4
        // write empty EPRL table and continue;
        if (kh_size(h) == 0) {
            
            // write EMPTY hash table for this i4 to disk
            kh_write(HashTableEPRL, h3, path_ampl);

            // free memory
            kh_destroy(HashTableJ6, h);
            kh_destroy(HashTableEPRL, h3);

            continue;
            
        }


        dspin two_i1, two_i2, two_i3, two_i5;

        int ret; // return codes from khash

        for (two_i5 = two_i5_min; two_i5 <= two_i5_max; two_i5 += 2) {
        for (two_i3 = two_i3_min; two_i3 <= two_i3_max; two_i3 += 2) {
        for (two_i2 = two_i2_min; two_i2 <= two_i2_max; two_i2 += 2) {
        for (two_i1 = two_i1_min; two_i1 <= two_i1_max; two_i1 += 2) {


            // key for disk cache of current amplitude
            HashTableEPRL_key_t keyEPRL = {two_j1, two_j2, two_j3, two_j4, two_j7, two_j8,
                                           two_i1, two_i2, two_i3, two_i5};

            // check if value has been already computed and stored
            if (kh_get(HashTableEPRL, h3, keyEPRL) != kh_end(h3)) {

                // value found!
                // go to next intertwiner set
                continue;

            }

            //////////////////////////////////////////////////////////////////////
            // value was not found or table was initialized
            // compute it then
            //////////////////////////////////////////////////////////////////////

            if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

                // with reducible coupling amplitude there are
                // additional triangular inequalities on i1

                dspin two_i1_min_red = (dspin) max(abs(two_i2-two_i3), abs(two_i5-two_i4));
                dspin two_i1_max_red = (dspin) min(two_i2+two_i3, two_i5+two_i4);

                if (!(max(two_i1_min, two_i1_min_red) <= two_i1 && two_i1 <= min(two_i1_max, two_i1_max_red))) {

                    // save 0 and pass to next

                    khint_t s = kh_put(HashTableEPRL, h3, keyEPRL, &ret);

                    if (ret == -1) {
                        error("error inserting key into EPRL hash table");
                    }

                    if (ret == 1) {
                        kh_val(h3, s) = 0.0;
                    }

                    continue;

                }

            } 

            double ampl;

            dspin *A = malloc(6 * sizeof(dspin));

            // Compute 15j Symbol
            ampl = J15Symbol(h, &A,
                             two_j1, two_j2, two_j3, two_j4, two_j5,
                             two_j6, two_j7, two_j8, two_j9, two_j10,
                             two_i1, two_i2, two_i3, two_i4, two_i5);

            free(A);

            // normalize amplitude with dimensions of all 5 boundary intertwiners
            ampl = ampl * d(two_i1) * d(two_i2) * d(two_i3) * d(two_i4) * d(two_i5);


            ////////////////////////////////////////////////////////
            // put key into amplitudes cache to be written to disk
            ////////////////////////////////////////////////////////

            khint_t s = kh_put(HashTableEPRL, h3, keyEPRL, &ret);

            if (ret == -1) {
                error("error inserting key into EPRL hash table");
            }

            if (ret == 1) {

                kh_val(h3, s) = ampl;

                // check if insertion worked
                s = kh_get(HashTableEPRL, h3, keyEPRL);
                if (kh_val(h3, s) == 0 && ampl != 0) {

                    // retry
                    kh_val(h3, s) = ampl;

                    // check again and now if not fail
                    if (kh_val(h3, s) == 0 && ampl != 0) {
                        error("error inserting value into EPRL hashtable");
                    }

                }

            }

            //////////////////////////////////////////////////////////////
            // IF the coupling is IRREDUCIBLE
            // AND IF all the boundary spins are equal
            // store also the companion amplitude obtained by reflecting
            // the pentachoron on the axis crossing i4
            //////////////////////////////////////////////////////////////

            if (get_coupling() == SL2CFOAM_COUPLING_IRREDUCIBLE &&
                two_j1 == two_j2 && two_j2 == two_j3 && two_j3 == two_j4 && two_j4 == two_j5 &&
                two_j5 == two_j6 && two_j6 == two_j7 && two_j7 == two_j8 && two_j8 == two_j9 && two_j9 == two_j10) {

                // relected key i1 <-> i2 , i3 <-> i5
                HashTableEPRL_key_t keyEPRLref = {two_j1, two_j2, two_j3, two_j4, two_j7, two_j8,
                                                  two_i2, two_i1, two_i5, two_i3};

                s = kh_put(HashTableEPRL, h3, keyEPRLref, &ret);

                if (ret == -1) {
                    error("error inserting key into EPRL hash table");
                }

                if (ret == 1) {

                    kh_val(h3, s) = ampl;

                    // check if insertion worked
                    s = kh_get(HashTableEPRL, h3, keyEPRLref);
                    if (kh_val(h3, s) == 0 && ampl != 0) {

                        // retry
                        kh_val(h3, s) = ampl;

                        // check again and now if not fail
                        if (kh_val(h3, s) == 0 && ampl != 0) {
                            error("error inserting value into EPRL hashtable");
                        }

                    }

                }

            }

        } // i1
        } // i2
        } // i3
        } // i5

        // write hash table for this i4 to disk
        kh_write(HashTableEPRL, h3, path_ampl);

        // free memory
        kh_destroy(HashTableJ6, h);
        kh_destroy(HashTableEPRL, h3);

    } // i4

    sl2cfoam_free_thread();

    } // omp parallel 

}