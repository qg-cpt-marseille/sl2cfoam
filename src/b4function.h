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

#ifndef __SL2CFOAM_B4FUNCTION_H__
#define __SL2CFOAM_B4FUNCTION_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"

////////////////////////////////////////////////////////////////
// Number of points for integration. Default is N = 3000.
// Set 1000 for low precision or 5000 for high precision.
////////////////////////////////////////////////////////////////

  #define B4_INTEGRATION_POINTS 3000 // default precision
//#define B4_INTEGRATION_POINTS 1000 // low precison
//#define B4_INTEGRATION_POINTS 10000 // high precision

////////////////////////////////////////////////////////////////
/// Functions to compute booster functions with 4 legs.
////////////////////////////////////////////////////////////////

// Hashing function for B4s 
void B4_Hash(dspin two_j1,  dspin two_j2, dspin two_j3,
             dspin two_j4,  dspin two_j5, dspin two_j6,
             dspin two_j7,  dspin two_j8, dspin two_j9,
             dspin two_j10, dspin two_Dl, float Immirzi);

// Computes the functions using Collet's code.
long double **B4Function(int two_k1, int two_k2, int two_k3, int two_k4,
                         float two_rho1, float two_rho2, float two_rho3, float two_rho4,
                         int two_j1, int two_j2, int two_j3, int two_j4,
                         int two_l1, int two_l2, int two_l3, int two_l4);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_B4FUNCTION_H__*/
