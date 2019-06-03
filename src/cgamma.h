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

/**********************************************************************/

#ifndef __SL2CFOAM_CGAMMA_H__
#define __SL2CFOAM_CGAMMA_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"

////////////////////////////////////////////////////////////////////////
// Computes the logarithm of complex Gamma function
// in arbitrary precision (ab) using Lanczos methos.
////////////////////////////////////////////////////////////////////////

// Computes the logarithm of complex Gamma function of op and stores
// result in rop.
int complex_lngamma(mpc_t rop, mpc_t op);

// Always call this before using complex gamma.
void init_complex_gamma();

// Always call this when finished using complex gamma.
void clear_complex_gamma();

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_CGAMMA_H__*/
