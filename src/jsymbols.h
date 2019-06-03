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

#ifndef __SL2CFOAM_JSYMBOLS_H__
#define __SL2CFOAM_JSYMBOLS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"
#include "hashing.h"

///////////////////////////////////////////////////////////
// Initialization utilities for wigxjpf 3j, 6j library.
///////////////////////////////////////////////////////////

// Maximum allowed spin for recoupling symbols.
#define MAX_J 1000

// Global init, call when program starts.
void init_wigxjpf_global();
void clear_wigxjpf_global();

// Call every time you enter a thread and need wigxjpf.
void init_wigxjpf_thread();
void clear_wigxjpf_thread();

///////////////////////////////////////////////////////////
// Various function to compute and store 6j, 9j, 15j symbols.
///////////////////////////////////////////////////////////

// Canonicalization function for 6j symbols.
void WIGNER6J_REGGE_CANONICALISE(dspin **ret, uint32_t two_j1, uint32_t two_j2, uint32_t two_j3,
                                             uint32_t two_j4, uint32_t two_j5, uint32_t two_j6);

// 6J symbol Hash function.
void J6Symbol_Hash(kh_HashTableJ6_t *h, dspin **A,
                   dspin two_k1, dspin two_k3, dspin two_k2,
                   dspin two_l1, dspin two_l2, dspin two_l3);

// 9J symbol Hash function via 6j summation.
void J9Symbol_Hash_Sum(kh_HashTableJ6_t *h, dspin **A,
                       dspin two_k2, dspin two_k3, dspin two_k1,
                       dspin two_j10, dspin two_j9, dspin two_i4,
                       dspin two_l7, dspin two_l8, dspin two_k5);

// 15J symbol Hash function.
void J15Symbol_Hash(dspin two_j1, dspin two_j2, dspin two_j3,
                    dspin two_j4, dspin two_j5, dspin two_j6,
                    dspin two_j7, dspin two_j8, dspin two_j9,
                    dspin two_j10, dspin two_i4, dspin two_Dl);

// 15J (reducible) symbol function.
double J15Symbol(const kh_HashTableJ6_t *h, dspin **A,
                 dspin two_l1, dspin two_l2, dspin two_l3, dspin two_l4, dspin two_j5,
                 dspin two_j6, dspin two_l7, dspin two_l8, dspin two_j9, dspin two_j10,
                 dspin two_k1, dspin two_k2, dspin two_k3, dspin two_i4, dspin two_k5);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_JSYMBOLS_H__*/
