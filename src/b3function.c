/* Copyright 2018 Giorgio Sarno, Pietro Donà and Francesco Gozzini */

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

#include "recouplingSL2C.h"
#include "recursion.h"
#include "config.h"
#include "utilities.h"

#include "b3function.h"

double B3(dspin two_j1, dspin two_j2, dspin two_j3,
          dspin two_l1, dspin two_l2, dspin two_l3,
          float Immirzi) {

	if (two_l1 < two_j1 || two_l2 < two_j2 || two_l3 < two_j3) {
    return 0.0;
  }

  double phase1 = real_negpow((int)(two_j1 - two_j2 + two_j3));
  double phase2 = real_negpow((int)(two_l1 - two_l2 + two_l3));

  double complex X1 = sqrt(d(two_j3 ))*conj(chi_simp(two_j1, two_j2, two_j3,
                                                     two_j1, two_j2, two_j3,
                                                     Immirzi));

  double complex X2 = sqrt(d(two_l3 ))*chi_simp(two_l1, two_l2, two_l3,
                                                two_j1, two_j2, two_j3,
                                                Immirzi);

	mpc_t z,w;
	mpc_init2(z, MPBITS);
	mpc_init2(w, MPBITS);
	mpc_set_dc(z, X1, MPC_RNDNN);

	mpc_set_dc(w, X2, MPC_RNDNN);
	mpc_mul(z, z, w,MPC_RNDNN);

    double complex result = phase1*phase2*mpc_get_dc(z, MPC_RNDNN);

	mpc_clear(z);
	mpc_clear(w);

    return creal(result);

}

//function to build the b3 function 
//from the chi already computed
//and to store them in an hash table
void b3_recursion(  kh_HashTableBooster3_t *h, mpfr_t*** chi_array, 
                    mpfr_t *b3_value,
                    dspin two_j1, dspin two_j2, dspin two_j3,
                    dspin two_l1, dspin two_l2, dspin two_l3,
                    float Immirzi) {
 
    int ret;
    HashTableBooster3_key_t key = {two_l1, two_l2, two_l3};

    if (kh_get(HashTableBooster3, h, key) != kh_end(h)) return;

    khint_t k = kh_put(HashTableBooster3, h, key, &ret);

    double cutoff = 1e-30;

    // check for error in hashing the ket
    if (ret == -1) {
        error("error inserting key into Booster3 hash table");
    }

    // if the key is correct put the value
    if (ret == 1) {

      double phase1 = real_negpow((int)(two_j1 - two_j2 + two_j3));
      double phase2 = real_negpow((int)(two_l1 - two_l2 + two_l3));
      double dim = sqrt(d(two_j3)*d(two_l3));
      
      mpfr_set(*b3_value,chi_array[0][0][0],MPFR_RNDN);
      mpfr_mul(*b3_value,*b3_value,chi_array[two_l1-two_j1][two_l2-two_j2][two_l3-two_j3],MPFR_RNDN);
      mpfr_mul_d(*b3_value,*b3_value,dim,MPFR_RNDN);

      double result = phase1 * phase2 *  mpfr_get_d(*b3_value, MPC_RNDNN);
      if (fabs(result) < cutoff){
        kh_val(h,k) = 0.; 
        //printf("%i %i %i %17g \n",two_l1,two_l2,two_l3,kh_val(h,k));
      } else {    
        kh_val(h,k) = result;
        //printf("%i %i %i %17g \n",two_l1,two_l2,two_l3,kh_val(h,k));
      }
    } 
  

}

void b3_hash ( dspin two_j1, dspin two_j2, dspin two_j3,
               float Immirzi) {
    
  //////////////////// data folder check ////////////////////

  check_data_3simplex(Immirzi);

  struct data_folders* fd = get_data_folders(Immirzi);

  char filename[1024];
  char path_boost[1024];
  char path_boost_i[1024];

  dspin two_Dl =  MAX_2D_3SIMPLEX;

  // build path for b3 values
  sprintf(filename, "%i.%i.%i_%i.boost", two_j1, two_j2, two_j3, two_Dl);
  strcpy(path_boost, fd->threesimp_imm_boost);
  strcat(path_boost, filename);

  //////////////////// Initialize or load hash table ////////////////////

  if (file_exist(path_boost) == 0) {

      //if the table is not found initialize a new one
      khash_t(HashTableBooster3) *h =  NULL;

      h = kh_init(HashTableBooster3);
  
      // to hold chi values at arbitrary precision
      mpfr_t ***chi_array = (mpfr_t***) malloc(CHI_ARRAY*sizeof(mpfr_t**));
        
      //initialize mpfr array
      int i,j,k;

      for (i = 0; i < CHI_ARRAY; i++) {
          chi_array[i] = (mpfr_t**) malloc(CHI_ARRAY*sizeof(mpfr_t*));
          for (j = 0; j < CHI_ARRAY; j++) {
              chi_array[i][j] = (mpfr_t*) malloc(CHI_ARRAY*sizeof(mpfr_t));
          }
      }

      for (i = 0; i < CHI_ARRAY; i++) {
          for (j = 0; j < CHI_ARRAY; j++) {
              for (k = 0; k < CHI_ARRAY; k++) {
                mpfr_init2(chi_array[i][j][k], MPBITS);
                mpfr_set_zero(chi_array[i][j][k], 1);

              }
          }
      }

      // we start the computation of the seed    

      chi_mpc(&chi_array [0][0][0], two_j1, two_j2, two_j3,
                                    two_j1, two_j2, two_j3,
                                    Immirzi*two_j1, Immirzi*two_j2, Immirzi*two_j3);                   

      float two_rho_1 =  Immirzi*two_j1;
      float two_rho_2 =  Immirzi*two_j2;
      float two_rho_3 =  Immirzi*two_j3;

      //Initialize mpfr variable to store che new
      //chi value   
      mpfr_t chi_value; 
      mpfr_init2(chi_value, MPBITS);

      // write just one to file
 
      //We start to compute all possible chi combination
      //They have to be done in a precise order: the next
      //order is built starting from the previous one  
      for (dspin two_tower = 2; two_tower <=  two_Dl; two_tower += 2) {

          for (dspin two_l1 = two_j1; two_l1 <= two_j1+two_tower; two_l1 += 2) {
          for (dspin two_l2 = two_j2; two_l2 <= two_j2+two_tower; two_l2 += 2) {
          for (dspin two_l3 = max(abs(two_l1-two_l2), two_j3); two_l3 <= min(two_l1+two_l2,two_j3+two_tower); two_l3 += 2) {

              //checking for zero and for already computed values
              if (two_l1 - two_j1 == 0 && two_l2 - two_j2 == 0 && two_l3 - two_j3 == 0) continue;   
              if (mpfr_get_d(chi_array[two_l1 - two_j1][two_l2 - two_j2][two_l3 - two_j3],MPFR_RNDN) != 0) continue;

              //start the recursion to get our number  
              recursion ( &chi_value, chi_array, 
                          two_l1, two_l2, two_l3,
                          two_j1, two_j2, two_j3,
                          two_rho_1, two_rho_2, two_rho_3,
                          Immirzi );
              //printf("R %i %i %i %i %i %i %17g\n", two_j1, two_j2, two_j3,two_l1, two_l2, two_l3, mpfr_get_d(chi_value,MPFR_RNDN));
              
              //set the new value in the mpfr array
              mpfr_set(chi_array[two_l1 - two_j1][two_l2 - two_j2][two_l3 - two_j3], chi_value, MPFR_RNDN);          

          } // l1
          } // l2
          } // l3

      } // tower

      mpfr_clear(chi_value);

      //now that we have all possible chi values
      // we compute and store the b3 functions

      mpfr_t b3_value;
      mpfr_init2(b3_value,MPBITS);

      for (dspin two_l1 = two_j1; two_l1 <= two_j1+two_Dl; two_l1 += 2) {
          for (dspin two_l2 = two_j2; two_l2 <= two_j2+two_Dl; two_l2 += 2) {
              for (dspin two_l3 = max(abs(two_l1-two_l2), two_j3); two_l3 <= min(two_l1+two_l2,two_j3+two_Dl); two_l3 += 2){
            
                b3_recursion(   h, chi_array, &b3_value,
                                two_j1, two_j2, two_j3,
                                two_l1, two_l2, two_l3,
                                Immirzi) ;           

                /*printf("%i %i %i %i %i %i %17g \n",
                                two_j1, two_j2, two_j3,
                                two_l1, two_l2, two_l3,mpfr_get_d(b3_value,MPFR_RNDN));*/                     
              } //l1
          } //l2
      } //l3

    mpfr_clear(b3_value);

    //clear array
    for (i = 0; i < CHI_ARRAY; i++) {
        for (j = 0; j < CHI_ARRAY; j++) {
            for (k = 0; k < CHI_ARRAY; k++) {
              mpfr_clear(chi_array[i][j][k]);
            }
        }
    }

    for (i = 0; i < CHI_ARRAY; i++) {
        for (j = 0; j < CHI_ARRAY; j++) {
            free(chi_array[i][j]);
        }
    }

    for (i = 0; i < CHI_ARRAY; i++) {
        free(chi_array[i]);
    }

    free(chi_array);
    
    //write table and clear
    kh_write(HashTableBooster3,h,path_boost);
    kh_destroy(HashTableBooster3, h);

  } 
}

//function to get the value
//of a b3 in an hash table
//It considers permutations.
double b3_get_value(  kh_HashTableBooster3_t *h_1, 
                      dspin two_l1, dspin two_l2, dspin two_l3,
                      int permutation) {

    HashTableBooster3_key_t boost_key_1;

    //Three possible cyclical permutations
    // We have to get the l n the right order                
    if (permutation == 1) { boost_key_1 = (HashTableBooster3_key_t) {two_l1,two_l2,two_l3}; }
    if (permutation == 2) { boost_key_1 = (HashTableBooster3_key_t) {two_l3,two_l1,two_l2}; }
    if (permutation == 3) { boost_key_1 = (HashTableBooster3_key_t) {two_l2,two_l3,two_l1}; }

    khint_t s_1 = kh_get(HashTableBooster3, h_1, boost_key_1);  
    double Boost_1;
    if (s_1 != kh_end(h_1)) {
        Boost_1 = kh_val(h_1, s_1);
    } else{
        error("cannot find 3 booster key");
    }    

    return Boost_1;
    
}