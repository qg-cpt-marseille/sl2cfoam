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

#ifndef __SL2CFOAM_COMMON_H__
#define __SL2CFOAM_COMMON_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

///////////////////////////////////////////////////////////
// Common includes, macros and defines for the library.
///////////////////////////////////////////////////////////

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <complex.h>
#include <assert.h>
#include <math.h>
#include <mpfr.h>
#include <mpc.h>


// Simple macros.
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define d(tj) ((int)tj+1)

// Bool types.
typedef enum { false, true } bool;

// 32 bytes integer
typedef __uint32_t uint32_t;

// Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

///////////////////////////////////////////////////////////
// Spin native machine types.
//
// By convention, all spin arguments to functions are
// (positive) integers of the form 2 * spin.
// Nomenclature is as follows: a dspin variable that starts
// with 'two_' means that the variable holds 2 * spin value.
// Example: j == 0.5, two_j == 1.
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// IMPORTANT: replicate the same definitions in 
//            external include header sl2cfoam.h.
///////////////////////////////////////////////////////////             

// Machine type for double-spins.
typedef unsigned short dspin;

// Machine type for half-integer spins.
typedef float spin;

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_COMMON_H__*/
