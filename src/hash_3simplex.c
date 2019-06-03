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
#include "utilities.h"
#include "config.h"
#include "hashing.h"
#include "jsymbols.h"
#include "b3function.h"
#include "recursion.h"
#include "error.h"
#include "sl2cfoam.h"

// To hash 6j symbols, boosters and build the 3 simplex amplitude.
// We use this to build the single vertex amplitude and
// the multi vertices amplitudes. it has to be parallelized
// on j1,j2,j3 spins.
static void constructor_three_ampl(kh_HashTableEPRL3_t *h2,
                                   dspin two_j1, dspin two_j2, dspin two_j3, 
                                   dspin two_j4, dspin two_j5, dspin two_j6,
                                   dspin two_Dl, float Immirzi) {
    
    //key for the amplitude
    HashTableEPRL3_key_t keyEPRL = {two_j4, two_j5, two_j6};

    //check if you have already computed the amplitude value
    if (kh_get(HashTableEPRL3, h2, keyEPRL) != kh_end(h2)) {
        // got it, exit from the constructor
        return;
    }

    //data folder check
    check_data_3simplex(Immirzi);

    // build paths
    struct data_folders* fd = get_data_folders(Immirzi);

    char filename[1024];
    char path_boost_1[1024];
    char path_boost_1a[1024];
    char path_boost_1b[1024];
    char path_boost_2[1024];
    char path_boost_2a[1024];
    char path_boost_2b[1024];
    char path_boost_3[1024];
    char path_boost_3a[1024];
    char path_boost_3b[1024];
    char path_6j[1024];
    char path_6j_dli[1024];

    // build path for boosts
    sprintf(filename, "%i.%i.%i_%i.boost",two_j6, two_j5, two_j1, MAX_2D_3SIMPLEX);
    strcpy(path_boost_1, fd->threesimp_imm_boost);
    strcat(path_boost_1, filename);
    sprintf(filename, "%i.%i.%i_%i.boost",two_j1, two_j6, two_j5, MAX_2D_3SIMPLEX);
    strcpy(path_boost_1a, fd->threesimp_imm_boost);
    strcat(path_boost_1a, filename);
    sprintf(filename, "%i.%i.%i_%i.boost",two_j5, two_j1, two_j6, MAX_2D_3SIMPLEX);
    strcpy(path_boost_1b, fd->threesimp_imm_boost);
    strcat(path_boost_1b, filename);

    sprintf(filename, "%i.%i.%i_%i.boost",two_j4, two_j2, two_j6, MAX_2D_3SIMPLEX);
    strcpy(path_boost_2, fd->threesimp_imm_boost);
    strcat(path_boost_2, filename);
    sprintf(filename, "%i.%i.%i_%i.boost",two_j6, two_j4, two_j2, MAX_2D_3SIMPLEX);
    strcpy(path_boost_2a, fd->threesimp_imm_boost);
    strcat(path_boost_2a, filename);
    sprintf(filename, "%i.%i.%i_%i.boost",two_j2, two_j6, two_j4, MAX_2D_3SIMPLEX);
    strcpy(path_boost_2b, fd->threesimp_imm_boost);
    strcat(path_boost_2b, filename);

    sprintf(filename, "%i.%i.%i_%i.boost",two_j3, two_j5, two_j4, MAX_2D_3SIMPLEX);
    strcpy(path_boost_3, fd->threesimp_imm_boost);
    strcat(path_boost_3, filename);
    sprintf(filename, "%i.%i.%i_%i.boost",two_j4, two_j3, two_j5, MAX_2D_3SIMPLEX);
    strcpy(path_boost_3a, fd->threesimp_imm_boost);
    strcat(path_boost_3a, filename);
    sprintf(filename, "%i.%i.%i_%i.boost",two_j5, two_j4, two_j3, MAX_2D_3SIMPLEX);
    strcpy(path_boost_3b, fd->threesimp_imm_boost);
    strcat(path_boost_3b, filename);

    // build path for 6js
    sprintf(filename, "%i.%i.%i_%i.6j", two_j1, two_j2, two_j3, two_Dl);
    strcpy(path_6j, fd->threesimp_hashtable6j);
    strcat(path_6j, filename);

    // Set Hash Tables for B3 
    khash_t(HashTableBooster3) *h_1 = NULL;
    khash_t(HashTableBooster3) *h_2 = NULL;
    khash_t(HashTableBooster3) *h_3 = NULL;

    // now we implement a strategy to use
    // symmetry relations in booster functions
    // to avoid recomputing a seed.

    // indent the three possible cyclic permutations
    // and load table if they are already present


    int permutation_1 = 0;

    if (file_exist(path_boost_1) == 1) { 
        permutation_1 = 1;
        h_1 = kh_load(HashTableBooster3, path_boost_1);
    } else if (file_exist(path_boost_1a) == 1) {
        permutation_1 = 2;
        h_1 = kh_load(HashTableBooster3, path_boost_1a);
    } else if (file_exist(path_boost_1b) == 1) {
        permutation_1 = 3;
        h_1 = kh_load(HashTableBooster3, path_boost_1b);
    }

    int permutation_2 = 0;
    if (file_exist(path_boost_2) == 1) {
        permutation_2 = 1;
        h_2 = kh_load(HashTableBooster3, path_boost_2);
    } else if (file_exist(path_boost_2a) == 1) {
        permutation_2 = 2;
        h_2 = kh_load(HashTableBooster3, path_boost_2a);
    } else if (file_exist(path_boost_2b) == 1) {
        permutation_2 = 3;
        h_2 = kh_load(HashTableBooster3, path_boost_2b);
    }
    
    int permutation_3 = 0;
    if (file_exist(path_boost_3) == 1) {
        permutation_3 = 1;
        h_3 = kh_load(HashTableBooster3, path_boost_3);
    } else if (file_exist(path_boost_3a) == 1) {
        permutation_3 = 2;
        h_3 = kh_load(HashTableBooster3, path_boost_3a);
    } else if (file_exist(path_boost_3b) == 1) {
        permutation_3 = 3;
        h_3 = kh_load(HashTableBooster3, path_boost_3b);
    }

    // if the tables are not present compute values
    // and then load the table
    if (permutation_1 == 0) {     
        b3_hash (two_j6, two_j5, two_j1,Immirzi);
        permutation_1 = 1;
        h_1 = kh_load(HashTableBooster3, path_boost_1);
    }
    if (permutation_2 == 0) {     
        b3_hash (two_j4, two_j2, two_j6,Immirzi);
        permutation_2 = 1;
        h_2 = kh_load(HashTableBooster3, path_boost_2);
    } 
    if (permutation_3 == 0) {     
        b3_hash (two_j3, two_j5, two_j4,Immirzi);
        permutation_3 = 1;
        h_3 = kh_load(HashTableBooster3, path_boost_3);
    }  

    // tables for 6j symbols
    khash_t(HashTableJ6) *h1 = NULL;
    for (int two_i = two_Dl; two_i >= 0; two_i -= 2) {

        sprintf(filename, "%i.%i.%i_%i.6j", two_j1, two_j2, two_j3, two_i);
        strcpy(path_6j_dli, fd->threesimp_hashtable6j);
        strcat(path_6j_dli, filename);

        if (file_exist(path_6j_dli) != 0) {
            h1 = kh_load(HashTableJ6, path_6j_dli);
            break;
        }

    }

    if (h1 == NULL) {
        h1 = kh_init(HashTableJ6);
    }

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    double value = 0.0;
    double err = 0.0;
    double cutoff = 1e-30;

    ///////////////////////////////////
    // Compute J6  and store them
    // in ArrayJ6 as well as in the HashTable
    ///////////////////////////////////

    dspin *A = malloc(6 * sizeof(dspin));

    dspin arrayj6_dim = MAX_2D_3SIMPLEX + 2;
    double ***ArrayJ6 = (double***) malloc(arrayj6_dim * sizeof(double**));
        
    //initialize J6 array
    int i,j;

    for (i = 0; i < CHI_ARRAY; i++) {
        ArrayJ6[i] = (double**) malloc(arrayj6_dim * sizeof(double*));
        for (j = 0; j < CHI_ARRAY; j++) {
            ArrayJ6[i][j] = (double*) malloc(arrayj6_dim * sizeof(double));
        }
    }

    for (dspin two_l4 = two_j4; two_l4 <= two_j4 + two_Dl; two_l4 += 2) {

    for (dspin two_l6 = max(two_j6, abs(two_l4-two_j2));
            two_l6 <= min(two_j6+two_Dl,two_l4+two_j2);
            two_l6 += 2) {

    for (dspin two_l5 = max(two_j5, max(abs(two_l4-two_j3),abs(two_l6-two_j1)));
               two_l5 <= min(two_j5+two_Dl, min(two_l4+two_j3,two_l6+two_j1));
               two_l5 += 2) {

        J6Symbol_Hash(h1, &A,
                      two_j1,  two_j2,  two_j3,
                      two_l4,  two_l5,  two_l6);

        HashTableJ6_key_t keyJ6 = {A[0],A[1],A[2],A[3],A[4],A[5]};
        khint_t s = kh_get(HashTableJ6, h1, keyJ6);

        ArrayJ6[two_l4-two_j4][two_l5-two_j5][two_l6-two_j6] = kh_value(h1,s);

    } // l5
    } // l6
    } // l4

    free(A);

    //Store 6j table
    kh_write(HashTableJ6, h1, path_6j);
    kh_destroy(HashTableJ6, h1);

    //////////////////////////////////////
    // Now we have all elements, start
    // to assembly the value and summing
    //////////////////////////////////////

    // TODO: to parallelize we need to create a J6_ARRAY for
    // 		 each thread.

    #pragma omp parallel reduction (+:value, err)
    {

        #pragma omp for
        for (dspin two_l4 = two_j4; two_l4 <= two_j4 + two_Dl; two_l4 += 2) {

            mpfr_t boosterMPFR_1, boosterMPFR_2, boosterMPFR_3, J6MPFR;

        for (dspin two_l6 = max (two_j6,abs(two_l4-two_j2));
                   two_l6 <= min(two_j6 + two_Dl,two_l4+two_j2);
                   two_l6 +=2) {

            // Retrieve B3_3
            double Boost_2 = b3_get_value( h_2, two_l4, two_j2, two_l6, permutation_2);
            if ( fabs(Boost_2) < cutoff ) continue;
            
            mpfr_init_set_d(boosterMPFR_2, Boost_2, MPFR_RNDN);

            for (dspin two_l5 = max(two_j5,max(abs(two_l4-two_j3),abs(two_l6-two_j1)));
                       two_l5 <= min(two_j5+two_Dl,min(two_l4+two_j3,two_l6+two_j1));
                       two_l5 += 2) {

                // Retrieve B3_1,B3_3
                double Boost_1 = b3_get_value(h_1, two_l6, two_l5, two_j1, permutation_1);
                double Boost_3 = b3_get_value(h_3, two_j3, two_l5, two_l4, permutation_3);
                                        
                if (fabs(Boost_1) < cutoff || fabs(Boost_3) < cutoff) continue;

                //printf("%i %i %i %17g %17g %17g \n",two_l4,two_l5,two_l6, Boost_1,Boost_2,Boost_3);

                double J6 = ArrayJ6[two_l4-two_j4][two_l5-two_j5][two_l6-two_j6];          

                // multiply values as mpfr variables
                mpfr_init_set_d(boosterMPFR_1, Boost_1, MPFR_RNDN);
                mpfr_init_set_d(boosterMPFR_3, Boost_3, MPFR_RNDN);
                mpfr_init_set_d(J6MPFR, J6, MPFR_RNDN);

                mpfr_mul(boosterMPFR_1, boosterMPFR_1, boosterMPFR_2, MPFR_RNDN);
                mpfr_mul(boosterMPFR_1, boosterMPFR_1, boosterMPFR_3, MPFR_RNDN);
                mpfr_mul(boosterMPFR_1, boosterMPFR_1, J6MPFR, MPFR_RNDN);

                // sum via compensated summation
                compsum_mpfr(&err, &value, boosterMPFR_1);

                // clear mpfr variables
                mpfr_clears(boosterMPFR_1, boosterMPFR_3, J6MPFR, NULL);

                }
            mpfr_clear(boosterMPFR_2);
            }
        }
    }
    
    // clear J6array
    for (i = 0; i < arrayj6_dim; i++) {
        for (j = 0; j < arrayj6_dim; j++) {
            free(ArrayJ6[i][j]);
        }
    }

    for (i = 0; i < arrayj6_dim; i++) {
        free(ArrayJ6[i]);
    }

    free(ArrayJ6);

    // check if amplitude has consistent value
    if (fabs(value) > 1.0) {
        error("amplitude has non-consistent value (greater than 1.0)");
    }

    //////////////////// Put key and value in the EPRL3 HashTable ////////////////////

    int ret;
    
    // put key and value in the EPRL HashTable
    khint_t s = kh_put(HashTableEPRL3, h2, keyEPRL, &ret);

    if (ret == -1) {
        error("error inserting key into EPRL hash table");
    }

    if (ret == 1) {

        kh_val(h2, s) = value;

        // check if insertion worked
        s = kh_get(HashTableEPRL3, h2, keyEPRL);
        if (kh_val(h2, s) == 0 && value != 0) {

            // retry
            kh_val(h2, s) = value;

            // check again and now if not fail
            if (kh_val(h2, s) == 0 && value != 0) {
                error("error inserting value into EPRL hashtable");
            }
        }

    }

    // free memory
    kh_destroy(HashTableBooster3, h_1);
    kh_destroy(HashTableBooster3, h_2);
    kh_destroy(HashTableBooster3, h_3);

}

void sl2cfoam_hash_three_ampl(dspin two_j1_min, dspin two_j1_max,
                              dspin two_j2_min, dspin two_j2_max,
                              dspin two_j3_min, dspin two_j3_max,
                              dspin two_j4_min, dspin two_j4_max,
                              dspin two_j5_min, dspin two_j5_max,
                              dspin two_j6_min, dspin two_j6_max,                                                                                                        
                              dspin two_Dl, float Immirzi) {

    // data folder check
    check_data_3simplex(Immirzi);

    // build paths
    struct data_folders* fd = get_data_folders(Immirzi);

    char filename[1024];
    char path_ampl[1024];
    
    for (dspin two_j1 = two_j1_min; two_j1 <= two_j1_max; two_j1 += 2) {
    for (dspin two_j2 = two_j2_min; two_j2 <= two_j2_max; two_j2 += 2) {
    for (dspin two_j3 = two_j3_min; two_j3 <= two_j3_max; two_j3 += 2) {

        // checking non triangular values                                    
        if (two_j3 < abs(two_j1 - two_j2) || two_j3 > two_j1 + two_j2) {    
            continue;
        }

        // build path for amplitude values
        sprintf(filename, "%i.%i.%i_%i.eprl", two_j1, two_j2, two_j3, two_Dl);
        strcpy(path_ampl, fd->threesimp_imm_ampl);
        strcat(path_ampl, filename);

        //initialize eprl table
        khash_t(HashTableEPRL3) *h2 = NULL;

        // Check already computed values or tables
        if (file_exist(path_ampl) != 0) {
            //load the table  
            h2 = kh_load(HashTableEPRL3, path_ampl);
        } else {        
            // table not found, initialize it
            h2 = kh_init(HashTableEPRL3);
        }

        for (dspin two_j4 = two_j4_min; two_j4 <= two_j4_max; two_j4 += 2) {
        for (dspin two_j5 = two_j5_min; two_j5 <= two_j5_max; two_j5 += 2) {
        for (dspin two_j6 = two_j6_min; two_j6 <= two_j6_max; two_j6 += 2) {

        //checking non triangular values                                    
        if (two_j4 < abs(two_j3 - two_j5) || two_j4 > two_j3 + two_j5 ||
            two_j5 < abs(two_j1 - two_j6) || two_j5 > two_j1 + two_j6 ||
            two_j6 < abs(two_j2 - two_j4) || two_j6 > two_j2 + two_j4) {
                continue;
        }

            constructor_three_ampl(h2,
                                    two_j1, two_j2, two_j3, 
                                    two_j4, two_j5, two_j6,
                                    two_Dl, Immirzi);

        } // two_j6
        } // two_j5
        } // two_j4

        // store amplitudes
        kh_write(HashTableEPRL3, h2, path_ampl);

        // and clean
        kh_destroy(HashTableEPRL3, h2);
    
    } // two_j3
    } // two_j2
    } // two_j1
    
}

void sl2cfoam_hashall_three_ampl(dspin two_minj, dspin two_maxj,
                                 dspin two_Dl, float Immirzi) {

    // data folder check
    check_data_3simplex(Immirzi);

    // build paths
    struct data_folders* fd = get_data_folders(Immirzi);

    char filename[1024];
    char path_ampl[1024];

    init_wigxjpf_global();
 
    for (dspin two_j1 = two_minj; two_j1 <= two_maxj; two_j1 +=2) {
    for (dspin two_j2 = two_minj; two_j2 <= two_maxj; two_j2 +=2) {

        dspin two_j3_min = abs(two_j1-two_j2); 
        dspin two_j3_max = two_j1+two_j2;

    for (dspin two_j3 = two_j3_min; two_j3 <= two_j3_max; two_j3 +=2) {

        // build path for amplitude values
        sprintf(filename, "%i.%i.%i_%i.eprl", two_j1, two_j2, two_j3, two_Dl);
        strcpy(path_ampl, fd->threesimp_imm_ampl);
        strcat(path_ampl, filename);

        //initialize eprl table
        khash_t(HashTableEPRL3) *h2 = NULL;

        // Check already computed values or tables
        if (file_exist(path_ampl) != 0) {
            //load the table  
            h2 = kh_load(HashTableEPRL3, path_ampl);
        } else {        
            // table not found, initialize it
            h2 = kh_init(HashTableEPRL3);
        }

        dspin two_j4, two_j5, two_j6;
        dspin two_j5_min, two_j5_max, two_j6_min, two_j6_max;

        for (two_j4 = two_minj; two_j4 <= two_maxj; two_j4 += 2) {

            two_j5_min = abs(two_j4-two_j3); 
            two_j5_max = two_j4+two_j3;

        for (two_j5 = two_j5_min; two_j5 <= two_j5_max; two_j5 += 2) {

            two_j6_min = max(abs(two_j1-two_j5),abs(two_j2-two_j4));
            two_j6_max = min(two_j1+two_j5,two_j2+two_j4);

        for (two_j6 = two_j6_min; two_j6 <= two_j6_max; two_j6 += 2){

            constructor_three_ampl(h2,
                                   two_j1, two_j2, two_j3, 
                                   two_j4, two_j5, two_j6,
                                   two_Dl, Immirzi);
            //printf("%i.%i.%i.%i.%i.%i \n", two_j1, two_j2, two_j3, two_j4, two_j5, two_j6 );

        } // two_j6
        } // two_j5
        } // two_j4

        // store amplitudes
        kh_write(HashTableEPRL3, h2, path_ampl);

        // and clean
        kh_destroy(HashTableEPRL3, h2);
    
    } // two_j3
    } // two_j2
    } // two_j1

    clear_wigxjpf_global();
}