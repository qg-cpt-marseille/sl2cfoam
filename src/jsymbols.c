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

#include <sys/stat.h>
#include <omp.h>

#include "jsymbols.h"
#include "config.h"
#include "utilities.h"
#include "wigxjpf.h"
#include "error.h"
#include "sl2cfoam.h"


void init_wigxjpf_global() {
    wig_table_init(2*MAX_J, 6);
    wig_temp_init(2*MAX_J);
}

void clear_wigxjpf_global() {
    wig_temp_free();
    wig_table_free();
}

void init_wigxjpf_thread() {

    int tid = omp_get_thread_num();
	if (tid != 0) {
		// not in master thread
		wig_thread_temp_init(2*MAX_J);
	}

}

void clear_wigxjpf_thread() {

    int tid = omp_get_thread_num();
	if (tid != 0) {
		// not in master thread
		wig_temp_free();
	}

}


#define COLLECT_NEGATIVE(two_j1,two_j2,two_j3) do {  \
collect_sign |= (two_j1) | (two_j2) | (two_j3);      \
} while (0)

#define COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1,two_j2,two_j3) do {    \
collect_sign |= (((two_j2) + (two_j3)) - (two_j1));                 \
collect_sign |= ((two_j1) - ((two_j2) - (two_j3)));                 \
collect_sign |= ((two_j1) - ((two_j3) - (two_j2)));                 \
} while (0)

void WIGNER6J_REGGE_CANONICALISE(dspin **ret, uint32_t two_j1, uint32_t two_j2, uint32_t two_j3,
                                              uint32_t two_j4, uint32_t two_j5, uint32_t two_j6) {

    uint32_t b1 = (two_j1 + two_j2 + two_j3);
    uint32_t b2 = (two_j1 + two_j5 + two_j6);
    uint32_t b3 = (two_j4 + two_j2 + two_j6);
    uint32_t b4 = (two_j4 + two_j5 + two_j3);

    /* Check trivial-0 */
    uint32_t collect_sign = 0;
    uint32_t collect_odd = 0;

    collect_odd = b1 | b2 | b3 | b4;

    COLLECT_NEGATIVE(two_j1, two_j2, two_j3);
    COLLECT_NEGATIVE(two_j4, two_j5, two_j6);
    COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j2, two_j3);
    COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j5, two_j6);
    COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j4, two_j2, two_j6);
    COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j4, two_j5, two_j3);

    if ((collect_sign & (1 << (sizeof (int) * 8 - 1))) |
            (collect_odd & 1)){
        for (int i=0; i<=5; i++){
                (*ret)[i]=0;
        }
    }

    /* Check trivial-0 end */

    # define SHIFT_DOWN_1(a)  ((a) >> 1)

    b1 = SHIFT_DOWN_1(b1);
    b2 = SHIFT_DOWN_1(b2);
    b3 = SHIFT_DOWN_1(b3);
    b4 = SHIFT_DOWN_1(b4);

    uint32_t a1 = SHIFT_DOWN_1(two_j1 + two_j2 + two_j4 + two_j5);
    uint32_t a2 = SHIFT_DOWN_1(two_j1 + two_j3 + two_j4 + two_j6);
    uint32_t a3 = SHIFT_DOWN_1(two_j2 + two_j3 + two_j5 + two_j6);

    #define SWAP_TO_FIRST_LARGER(filenametype,a,b) do {  \
    filenametype __filename_a = a;				                    \
    filenametype __filename_b = b;				                    \
    a = (__filename_a > __filename_b) ? __filename_a : __filename_b;	  \
    b = (__filename_a < __filename_b) ? __filename_a : __filename_b;	  \
    } while (0)

    SWAP_TO_FIRST_LARGER(uint32_t, a1, a2);
    SWAP_TO_FIRST_LARGER(uint32_t, a2, a3);
    SWAP_TO_FIRST_LARGER(uint32_t, a1, a2);

    SWAP_TO_FIRST_LARGER(uint32_t, b1, b2);
    SWAP_TO_FIRST_LARGER(uint32_t, b2, b3);
    SWAP_TO_FIRST_LARGER(uint32_t, b3, b4);
    SWAP_TO_FIRST_LARGER(uint32_t, b1, b2);
    SWAP_TO_FIRST_LARGER(uint32_t, b2, b3);
    SWAP_TO_FIRST_LARGER(uint32_t, b1, b2);

    // We now have a1 >= a2 >= a3, and b1 >= b2 >= b3 >= b4

    uint32_t S32 = a3 - b1;
    uint32_t B32 = a3 - b2;
    uint32_t T32 = a3 - b3;
    uint32_t X32 = a3 - b4;
    uint32_t L32 = a2 - b4;
    uint32_t E32 = a1 - b4;

    (*ret)[0] = (dspin)E32;
    (*ret)[1] = (dspin)L32;
    (*ret)[2] = (dspin)X32;
    (*ret)[3] = (dspin)T32;
    (*ret)[4] = (dspin)B32;
    (*ret)[5] = (dspin)S32;

}

#undef SWAP_TO_FIRST_LARGER


void J6Symbol_Hash(kh_HashTableJ6_t *h, dspin **A,
                   dspin two_k1, dspin two_k3, dspin two_k2,
                   dspin two_l1, dspin two_l2, dspin two_l3) {

    //////////////////// Convert spins to key ////////////////////

    WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k1, (uint32_t)two_k3, (uint32_t)two_k2, 
                                   (uint32_t)two_l1, (uint32_t)two_l2, (uint32_t)two_l3);

    HashTableJ6_key_t keyA = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};

    // Check we haven't already the function
    if (kh_get(HashTableJ6, h, keyA) == kh_end(h)){

        //////////////////// Save value ////////////////////
        int retA;
        khint_t kA = kh_put(HashTableJ6, h, keyA, &retA);

        // check for no error
        if (retA == -1) {
            error("error inserting key into J6 hash table");
        }

        if ( retA == 1 ){

            double val6jA = wig6jj(two_k1, two_k3, two_k2, two_l1, two_l2, two_l3);
            kh_val(h,kA) = val6jA;
            kA = kh_get(HashTableJ6, h, keyA);
            
            if (kh_val(h,kA) == 0 && val6jA != 0) {

                // retry
                kh_val(h,kA) = val6jA;

                // check again and now if not fail
                if (kh_val(h,kA) == 0 && val6jA != 0) {
                    error("error inserting value into J6 hash table");
                }

            }
        }
    }
}

void J9Symbol_Hash_Sum(kh_HashTableJ6_t *h, dspin **A,
                       dspin two_k2, dspin two_k3, dspin two_k1,
                       dspin two_j10, dspin two_j9, dspin two_i4,
                       dspin two_l7, dspin two_l8, dspin two_k5) {

    dspin two_limx1 = max(max(abs(two_k2-two_k5), abs(two_j10-two_l8)), abs(two_k3-two_i4));
    dspin two_limx2 = min(min(two_k2+two_k5,two_j10+two_l8), two_k3+two_i4);

    if (two_limx2 < two_limx1) {
        return;
    }

    ////////////////////////////
    // 9j via 6j summation
    ////////////////////////////

    for (dspin two_x = two_limx1; two_x <= two_limx2; two_x += 2) {

        //////////////////// Hash all 6j combinations for 9j ////////////////////
        J6Symbol_Hash(h, A,
                      two_k2, two_j10, two_l7,
                      two_l8, two_k5, two_x);

        J6Symbol_Hash(h, A,
                      two_k3, two_j9, two_l8,
                      two_j10, two_x, two_i4 );

        J6Symbol_Hash(h, A,
                      two_k1, two_i4, two_k5,
                      two_x, two_k2, two_k3);

    }
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void J15Symbol_Hash(dspin two_j1, dspin two_j2, dspin two_j3,
                    dspin two_j4, dspin two_j5, dspin two_j6,
                    dspin two_j7, dspin two_j8, dspin two_j9,
                    dspin two_j10, dspin two_i4, dspin two_Dl) {

    dspin two_limi4_1, two_limi4_2;

    two_limi4_1 = max(abs(two_j5-two_j6),abs(two_j10-two_j9));
    two_limi4_2 = min(two_j5+two_j6,two_j10+two_j9);

    if (two_i4 < two_limi4_1 && two_i4 > two_limi4_2) {
        error("error with intertwiners on boundary");
    }

    // get subfolders structure
    struct data_folders* fd = get_data_folders(1.20); // fictitious Immirzi, not needed here

    char filename[1024];
    char path_6j[1024];
    char path_6j_dli[1024];
    char path_key15j[1024];
    char path_key15j_dli[1024];

    // build path for 6j symbols
    sprintf(filename, "%i.%i.%i.%i_%i_%i.6j", two_j5, two_j6, two_j9, two_j10, two_i4, two_Dl);
    strcpy(path_6j, fd->foursimp_hashtable6j_6j);
    strcat(path_6j, filename);

    // build path for 15j keys
    sprintf(filename, "%i.%i.%i.%i_%i_%i.k15j", two_j5, two_j6, two_j9, two_j10, two_i4, two_Dl);
    strcpy(path_key15j, fd->foursimp_aux_hashtable15j_15j);
    strcat(path_key15j, filename);

    //////////////////// Hash Table initialization ///////////////////

    // check for previously computed tables.
    khash_t(HashTableJ6) *h = NULL;
    int data_yes = 0;
    int data_check = 0;

    for (int two_i = two_Dl;  two_i >= 0; two_i -= 2) {

        sprintf(filename, "%i.%i.%i.%i_%i_%i.6j", two_j5, two_j6, two_j9, two_j10, two_i4, two_i);
        strcpy(path_6j_dli, fd->foursimp_hashtable6j_6j);
        strcat(path_6j_dli, filename);

        if (file_exist(path_6j_dli) != 0) {

            // Data checks are used in case a Key hash table is present
            // but the actual data have been removed for some reason
            data_yes = 1;
            data_check = two_i;
            h = kh_load(HashTableJ6, path_6j_dli);
            break;

        }

    }

    if (h == NULL) {
        h = kh_init(HashTableJ6);
    }

    // initialize a table with keys to know
    // what 15j pieces have been already computed

    khash_t(HashTableKeyJ15) *h2 = NULL;

    for (int two_i = two_Dl; two_i >= 0; two_i -= 2) {

        sprintf(filename, "%i.%i.%i.%i_%i_%i.k15j", two_j5, two_j6, two_j9, two_j10, two_i4, two_i);
        strcpy(path_key15j_dli, fd->foursimp_aux_hashtable15j_15j);
        strcat(path_key15j_dli, filename);

        if (file_exist(path_key15j_dli) != 0 ) {
            h2 = kh_load(HashTableKeyJ15, path_key15j_dli);
            break;
        }

    }

    if (h2 == NULL) {
        h2 = kh_init(HashTableKeyJ15);
    }

    //////////////////// Check if I have the {6j}s. If not compute ////////////////////
    HashTableKeyJ15_key_t key = {two_j1, two_j2, two_j3, two_j4, two_j7, two_j8};

    if (data_yes == 1 && data_check == two_Dl &&
        kh_get(HashTableKeyJ15, h2, key) != kh_end(h2)) {
        
        kh_destroy(HashTableKeyJ15, h2);
        kh_destroy(HashTableJ6, h);
        return;

    }

    int ret;
    khint_t k = kh_put(HashTableKeyJ15, h2, key, &ret);

    // check for no error
    if (ret == -1) {
        error("error inserting key into aux J15 hash table");
    }


    //////////////////// Initialize pointer for J6 symbols ////////////////////

    dspin* A = malloc(6 * sizeof(dspin));

    //////////////////// Start cycling all possible L combinations ////////////////////

    dspin two_l1, two_l2, two_l3, two_l4, two_l7, two_l8,
          two_k2, two_k3, two_k5, two_k1;

    for(two_l1 = two_j1; two_l1 <= two_j1+two_Dl; two_l1 += 2) {
    for(two_l2 = two_j2; two_l2 <= two_j2+two_Dl; two_l2 += 2) {
    for(two_l3 = two_j3; two_l3 <= two_j3+two_Dl; two_l3 += 2) {
    for(two_l4 = two_j4; two_l4 <= two_j4+two_Dl; two_l4 += 2) {
    for(two_l7 = two_j7; two_l7 <= two_j7+two_Dl; two_l7 += 2) {
    for(two_l8 = two_j8; two_l8 <= two_j8+two_Dl; two_l8 += 2) {

        dspin two_limk1_1, two_limk1_2, two_limk2_1, 
              two_limk2_2, two_limk3_1, two_limk3_2,
              two_limk5_1, two_limk5_2;


        if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

            ///////////////////////////////////////
            // using REDUCIBLE coupling
            ///////////////////////////////////////

            two_limk2_1 = max(abs(two_j10-two_l7),abs(two_l1-two_l2));
            two_limk2_2 = min(two_j10+two_l7,two_l1+two_l2);
            two_limk3_1 = max(abs(two_j9-two_l8),abs(two_l3-two_l1));
            two_limk3_2 = min(two_j9+two_l8,two_l3+two_l1);
            two_limk5_1 = max(abs(two_l4-two_j6),abs(two_l7-two_l8));
            two_limk5_2 = min(two_l4+two_j6,two_l7+two_l8);

            //////////////////// Start cycling all possible K combinations ////////////////////

            for(two_k2 = two_limk2_1; two_k2 <= two_limk2_2; two_k2+=2){
            for(two_k3 = two_limk3_1; two_k3 <= two_limk3_2; two_k3+=2){
            for(two_k5 = two_limk5_1; two_k5 <= two_limk5_2; two_k5+=2){
            for(two_k1 = max(max(abs(two_l3-two_l2),abs(two_j5-two_l4)),max(abs(two_k3-two_k2),abs(two_k5-two_i4)));
                two_k1 <= min(min(two_l3+two_l2,two_j5+two_l4),min(two_k2+two_k3,two_k5+two_i4)); two_k1+=2){

                ///////////////////////////////////////
                // 6j Symbols - Key conversion and save
                ///////////////////////////////////////

                J6Symbol_Hash(h, &A,
                              two_k1, two_k3, two_k2,
                              two_l1, two_l2, two_l3);
                J6Symbol_Hash(h, &A,
                              two_k1, two_i4, two_k5,
                              two_j6, two_l4, two_j5);

                ///////////////////////////////////////////////////////
                // 9j Symbols via 6j Summation - Key conversion and save
                ///////////////////////////////////////////////////////

                J9Symbol_Hash_Sum(h, &A,
                                  two_k2, two_k3, two_k1,
                                  two_j10, two_j9, two_i4,
                                  two_l7, two_l8, two_k5);

            }
            }
            }
            }

        } else {

            ///////////////////////////////////////
            // using IRREDUCIBLE coupling
            ///////////////////////////////////////

            two_limk1_1 = max(abs(two_l2-two_l3), abs(two_j5-two_l4));
            two_limk1_2 = min(two_l2+two_l3, two_j5+two_l4);
            two_limk2_1 = max(abs(two_l1-two_j10), abs(two_l7-two_l2));
            two_limk2_2 = min(two_l1+two_j10, two_l7+two_l2);
            two_limk3_1 = max(abs(two_j9-two_l8), abs(two_l3-two_l1));
            two_limk3_2 = min(two_j9+two_l8, two_l3+two_l1);
            two_limk5_1 = max(abs(two_l4-two_l7), abs(two_l8-two_j6));
            two_limk5_2 = min(two_l4+two_l7, two_l8+two_j6);

            dspin two_limx1, two_limx2;

            for(two_k1 = two_limk1_1 ;two_k1 <= two_limk1_2; two_k1+=2){
            for(two_k2 = two_limk2_1; two_k2 <= two_limk2_2; two_k2+=2){
            for(two_k3 = two_limk3_1; two_k3 <= two_limk3_2; two_k3+=2){
            for(two_k5 = two_limk5_1; two_k5 <= two_limk5_2; two_k5+=2){

                two_limx1 = max(abs(two_k1-two_l7), abs(two_j5-two_k5));
                two_limx1 = max(two_limx1, abs(two_i4-two_l8));
                two_limx1 = max(two_limx1, abs(two_j10-two_k3));
                two_limx1 = max(two_limx1, abs(two_k2-two_l3));

                two_limx2 = min(two_k1+two_l7, two_j5+two_k5);
                two_limx2 = min(two_limx2, two_i4+two_l8);
                two_limx2 = min(two_limx2, two_j10+two_k3);
                two_limx2 = min(two_limx2, two_k2+two_l3);

                if (two_limx2 < two_limx1) {
                    continue;
                }

                dspin two_x;
                for (two_x = two_limx1; two_x <= two_limx2; two_x += 2) {

                    J6Symbol_Hash(h, &A,
                                  two_k1, two_l7, two_x, 
                                  two_k5, two_j5, two_l4 );
                    J6Symbol_Hash(h, &A,
                                  two_j5, two_k5, two_x, 
                                  two_l8, two_i4, two_j6);                       
                    J6Symbol_Hash(h, &A,
                                  two_i4, two_l8, two_x, 
                                  two_k3, two_j10, two_j9);             
                    J6Symbol_Hash(h, &A,
                                  two_j10, two_k3, two_x, 
                                  two_l3, two_k2, two_l1);
                    J6Symbol_Hash(h, &A,
                                  two_k2, two_l3, two_x, 
                                  two_k1, two_l7, two_l2);

                }

            }
            }
            }
            }

        }

    } // l8
    } // l7
    } // l4
    } // l3
    } // l2
    } // l1

    free(A);

    // write hashtables to disk
    kh_write(HashTableKeyJ15, h2, path_key15j);
    kh_write(HashTableJ6, h, path_6j);

    // free memory
    kh_destroy(HashTableJ6, h);
    kh_destroy(HashTableKeyJ15, h2);

}


double J15Symbol(const kh_HashTableJ6_t *h, dspin **A,
                 dspin two_l1, dspin two_l2, dspin two_l3, dspin two_l4, dspin two_j5,
                 dspin two_j6, dspin two_l7, dspin two_l8, dspin two_j9, dspin two_j10,
                 dspin two_k1, dspin two_k2, dspin two_k3, dspin two_i4, dspin two_k5) {

    double val15j = 0.0;
    double cutoff = 1e-40;

    ////////////////////////////////
    // Recover 6j and Compute 15j
    ////////////////////////////////

    double val6jA, val6jB, val6jC, val6jD, val6jE;
    khiter_t kA, kB, kC, kD, kE;

    if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

        ///////////////////////////////////////
        // using REDUCIBLE coupling
        ///////////////////////////////////////

        WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k1, (uint32_t)two_k3, (uint32_t)two_k2,
                                       (uint32_t)two_l1, (uint32_t)two_l2, (uint32_t)two_l3);
        HashTableJ6_key_t keyA = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
        kA = kh_get(HashTableJ6, h, keyA);
        val6jA = kh_value(h, kA);

        if (fabs(val6jA) < cutoff) return 0;

        WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k1, (uint32_t)two_i4, (uint32_t)two_k5,
                                       (uint32_t)two_j6, (uint32_t)two_l4, (uint32_t)two_j5);
        HashTableJ6_key_t keyB = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
        kB = kh_get(HashTableJ6, h, keyB);
        val6jB = kh_value(h, kB);

        if (fabs(val6jB) < cutoff) return 0;

        dspin two_limx1 = max(max(abs(two_k2-two_k5),abs(two_j10-two_l8)),abs(two_k3-two_i4));
        dspin two_limx2 = min(min(two_k2+two_k5,two_j10+two_l8),two_k3+two_i4);

        if (two_limx2 < two_limx1) {
            return 0.0;
        }

        ////////////////////////////
        // 9j via 6j summation
        ////////////////////////////

        double val9j = 0.0;

        for (dspin two_x = two_limx1; two_x <= two_limx2; two_x += 2){

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k2, (uint32_t)two_j10, (uint32_t)two_l7,
                                           (uint32_t)two_l8, (uint32_t)two_k5,  (uint32_t)two_x);
            HashTableJ6_key_t keyC = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kC = kh_get(HashTableJ6, h, keyC);
            val6jC = kh_value(h, kC);

            if (fabs(val6jC) < cutoff) continue;

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k3, (uint32_t)two_j9, (uint32_t)two_l8,
                                           (uint32_t)two_j10,(uint32_t)two_x,  (uint32_t)two_i4);
            HashTableJ6_key_t keyD= {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kD = kh_get(HashTableJ6, h, keyD);
            val6jD = kh_value(h, kD);

            if (fabs(val6jD) < cutoff) continue;

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k1, (uint32_t)two_i4, (uint32_t)two_k5, 
                                           (uint32_t)two_x,  (uint32_t)two_k2, (uint32_t)two_k3);
            HashTableJ6_key_t keyE = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kE = kh_get(HashTableJ6, h, keyE);
            val6jE = kh_value(h, kE);

            if (fabs(val6jE) < cutoff) continue;

            val9j += real_negpow(2 * two_x) * d(two_x) * val6jC * val6jD * val6jE;

        }

        val15j = real_negpow(2*(two_l4+two_l3+two_k1)-(two_k2+two_k3+two_i4+two_k5))
                 * val6jA * val6jB * val9j;

    } else {

        ///////////////////////////////////////
        // using IRREDUCIBLE coupling
        ///////////////////////////////////////

        dspin two_x, two_limx1, two_limx2;

        two_limx1 = max(abs(two_k1-two_l7), abs(two_j5-two_k5));
        two_limx1 = max(two_limx1, abs(two_i4-two_l8));
        two_limx1 = max(two_limx1, abs(two_j10-two_k3));
        two_limx1 = max(two_limx1, abs(two_k2-two_l3));

        two_limx2 = min(two_k1+two_l7, two_j5+two_k5);
        two_limx2 = min(two_limx2, two_i4+two_l8);
        two_limx2 = min(two_limx2, two_j10+two_k3);
        two_limx2 = min(two_limx2, two_k2+two_l3);

        if (two_limx2 < two_limx1) {
            return 0.0;
        }

        for (two_x = two_limx1; two_x <= two_limx2; two_x += 2) {

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k1, (uint32_t)two_l7, (uint32_t)two_x,
                                           (uint32_t)two_k5, (uint32_t)two_j5, (uint32_t)two_l4);

            HashTableJ6_key_t keyA = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kA = kh_get(HashTableJ6, h, keyA);
            val6jA = kh_value(h, kA);

            if (fabs(val6jA) < cutoff) continue;

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_j5, (uint32_t)two_k5, (uint32_t)two_x, 
                                           (uint32_t)two_l8, (uint32_t)two_i4, (uint32_t)two_j6);  

            HashTableJ6_key_t keyB = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kB = kh_get(HashTableJ6, h, keyB);
            val6jB = kh_value(h, kB);

            if (fabs(val6jB) < cutoff) continue;

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_i4, (uint32_t)two_l8,  (uint32_t)two_x, 
                                           (uint32_t)two_k3, (uint32_t)two_j10, (uint32_t)two_j9);

            HashTableJ6_key_t keyC = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kC = kh_get(HashTableJ6, h, keyC);
            val6jC = kh_value(h, kC);

            if (fabs(val6jC) < cutoff) continue;

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_j10, (uint32_t)two_k3, (uint32_t)two_x, 
                                           (uint32_t)two_l3,  (uint32_t)two_k2, (uint32_t)two_l1);

            HashTableJ6_key_t keyD = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kD = kh_get(HashTableJ6, h, keyD);
            val6jD = kh_value(h, kD);

            if (fabs(val6jD) < cutoff) continue;

            WIGNER6J_REGGE_CANONICALISE(A, (uint32_t)two_k2, (uint32_t)two_l3, (uint32_t)two_x, 
                                           (uint32_t)two_k1, (uint32_t)two_l7, (uint32_t) two_l2);

            HashTableJ6_key_t keyE = {(*A)[0],(*A)[1],(*A)[2],(*A)[3],(*A)[4],(*A)[5]};
            kE = kh_get(HashTableJ6, h, keyE);
            val6jE = kh_value(h, kE);

            if (fabs(val6jE) < cutoff) continue;

            val15j += d(two_x) * val6jA * val6jB * val6jC * val6jD * val6jE;

        }

        int sign = real_negpow(two_l1 + two_l2 + two_l3 + two_l4 + two_j5  + 
                               two_j6 + two_l7 + two_l8 + two_j9 + two_j10 +
                               two_k1 + two_k2 + two_k3 + two_i4 + two_k5);

        val15j = sign * val15j;   

    }

    return val15j;

}