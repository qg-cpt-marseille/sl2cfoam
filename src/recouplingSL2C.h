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

#ifndef __SL2CFOAM_RECOUPLINGSL2C_H__
#define __SL2CFOAM_RECOUPLINGSL2C_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"

////////////////////////////////////////////////////////////////
// SL(2,C) Clebsch Gordan coefficients and symbols
// in arbitrary precision.
//
// IMPORTANT: initialize complex gamma when program starts.
////////////////////////////////////////////////////////////////

// SL2C C-G coefficient Chi.
double complex chi(dspin two_j1, dspin two_j2, dspin two_j3,
                   dspin two_k1, dspin two_k2, dspin two_k3,
                   double two_rho1, double two_rho2, double two_rho3);
// for chi at arbitrary precision
void chi_mpc(mpfr_t *chi_value, dspin two_j1, dspin two_j2, dspin two_j3,
                              dspin two_k1, dspin two_k2, dspin two_k3,
                              double two_rho1, double two_rho2, double two_rho3); 

// for chi at arbitrary precision with zeros
void chi_mpc_simp0( mpfr_t *chi_value, dspin two_j, dspin two_k,
                    float Immirzi);                                 

// Chi for simplified series.
double complex chi_simp(dspin two_j1, dspin two_j2, dspin two_j3,
                        dspin two_k1, dspin two_k2, dspin two_k3,
                        float Immirzi);

// Computes 3j symbol for SL(2,C).
double  wigner_3j_sl2c(dspin two_j1, dspin two_j2, dspin two_j3,
                       dspin two_k1, dspin two_k2, dspin two_k3,
                       float two_rho1, float two_rho2, float two_rho3);

// Computes 6j symbol for SL(2,C).
double  wigner_6j_sl2c(dspin two_k1, dspin two_k2, dspin two_k3,
                      dspin two_k4, dspin two_k5, dspin two_k6,
                      float two_rho1, float two_rho2, float two_rho3,
                      float two_rho4, float two_rho5, float two_rho6,
                      dspin two_Dj);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_RECOUPLINGSL2C_H__*/
