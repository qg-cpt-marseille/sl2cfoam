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


#ifndef __SL2CFOAM_UTILITIES_H__
#define __SL2CFOAM_UTILITIES_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"
#include "error.h"

/////////////////////////////////////////////////////////////////////////
// Filesystem utilities.
/////////////////////////////////////////////////////////////////////////

// Checks if files exists.
int file_exist(char* filename);

// Data Folder Check for 4-Simplex
void check_data_4simplex(float Immirzi);

// Data Folder Check for 4-Simplex
void check_data_3simplex(float Immirzi);


/////////////////////////////////////////////////////////////////////////
// Math utilities.
/////////////////////////////////////////////////////////////////////////

// Checks that the dspin value corresponds to an integer spin.
#define ensure_integer_spin(two_j) \
    { if ((two_j) % 2 != 0) { error("integer check failed"); } }

// Checks that the dspins satisfy the triangle inequality.
#define ensure_triangle(two_j1, two_j2, two_j3) \
	{ if (!(abs(two_j1 - two_j2) <= two_j3 && two_j3 <= (two_j1 + two_j2))) { error("triangle check failed"); } }

// Checks that the 3 dspins sum to an integer.
#define ensure_integer_3sum(two_j1, two_j2, two_j3) \
    { if ((two_j1+two_j2+two_j3) % 2 != 0) { error("integer sum check failed"); } }

// Multiplication of complex values as MPFR variables.
double complex mpfr_complex_mul(double complex coherentA, double complex coherentB);

// Division of complex values as MPFR variables.
double complex mpfr_complex_div(double complex coherentA, double complex coherentB);

// Compensated summation functions.
void compsum_mpfr(double* err, double* sum, mpfr_t toadd);

void compsum_mpfr_complex(double complex* err, double complex* sum,
                          mpfr_t toadd_Re, mpfr_t toadd_Im);


// Computes (-1)^j for semi-integer j.
static inline double complex complex_negpow(int two_j) {

    int k = two_j % 2; // i factor
    int j = two_j / 2;

    if (k == 1) {
        if (j % 2 == 0) {
                return I;
        }
        return -I;
    }

    if (k == -1) {
        if (j % 2 == 0) {
                return -I;
        }
        return I;
    }

    if (j % 2 == 0) {
            return 1;
    }
    return -1;

}

// Computes (-1)^j for integer j.
// For semi-integer j it returns the real part.
static inline int real_negpow(int tj) {

    ensure_integer_spin(tj);

    int j = tj / 2;

    if (j % 2 == 0) {
            return 1;
    }
    return -1;

}

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_UTILITIES_H__*/
