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


#ifndef __SL2CFOAM_COHERENTSTATES_H__
#define __SL2CFOAM_COHERENTSTATES_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"

//////////////////// Coherent Tetrahedron at intertwiner I ////////////////////

//This function requieres four addition int as input: they have to be +1 or -1 depending on
//the strand's orientation. An "ingoing" strand going from the boundary node
//to the intertwiner requires a -1 while an "outgoing" strand, from the intertwiner
// to the boundary node requires a +1.

double complex CoherentStateI(dspin two_i, dspin two_j1, dspin two_j2,dspin two_j3,dspin two_j4,
                              double f1, double t1, double f2, double t2, double f3, double t3, double f4, double t4,
                              int sign1, int sign2, int sign3, int sign4 );

// Computes SU(2) Wigner D-matrix elements.
double complex WignerD(dspin two_j, int two_m, double f1, double t1, int sign);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_COHERENTSTATES_H__*/
