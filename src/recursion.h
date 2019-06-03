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

#ifndef __SL2CFOAM_RECURSION_H__
#define __SL2CFOAM_RECURSION_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"

////////////////////////////////////////////////////////////////
// Hashing of seeds up to a certain spin
// The seeds can then by use by recursion relations.
////////////////////////////////////////////////////////////////

// Define the maximum number of shell
#define MAX_2D_3SIMPLEX 100

// Define the array dimension
#define CHI_ARRAY (MAX_2D_3SIMPLEX+2)

// We use this function to create chi functions via
// recursion relations with precomputed seeds.
void recursion(mpfr_t* chi_value, mpfr_t*** chi_array,  
               dspin two_l1, dspin two_l2, dspin two_l3,
               dspin two_k1, dspin two_k2, dspin two_k3,
               float two_rho1, float two_rho2, float two_rho3,
               float Immirzi );

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_RECURSION_H__*/
