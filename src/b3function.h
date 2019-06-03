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

#ifndef __SL2CFOAM_B3FUNCTION_H__
#define __SL2CFOAM_B3FUNCTION_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "hashing.h"

////////////////////////////////////////////////////////////////
// Boosters3 function to be used to compute
// 3-simplices.
////////////////////////////////////////////////////////////////

// Computes the B3 booster function using
// chi and k formulas
double B3(dspin two_j1, dspin two_j2, dspin two_j3,
          dspin two_l1, dspin two_l2, dspin two_l3,
          float Immirzi);

// Hash the B3 booster function using recursion relations    
void b3_hash ( dspin two_j1, dspin two_j2, dspin two_j3,
               float Immirzi) ;

//function to get the value
//of a b3 in an hash table
//It considers permutations.
double b3_get_value(  kh_HashTableBooster3_t *h_1, 
                      dspin two_l1, dspin two_l2, dspin two_l3,
                      int permutation);        

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_B3FUNCTION_H__*/
