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

#ifndef __SL2CFOAM_HASHING_H__
#define __SL2CFOAM_HASHING_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "khash.h"

///////////////////////////////////////////////////////////
// Definition of Hash tables.
///////////////////////////////////////////////////////////

// 6 spins keys with double values for 6j symbols
ARRAY_TABLE_INIT(HashTableJ6, 6, double);

// 10 spins keys with double values for boosters
ARRAY_TABLE_INIT(HashTableBooster, 10, double);

// 6 spins keys for auxiliary 15j pieces keys
ARRAY_TABLE_INIT(HashTableKeyJ15, 6, int);

// 6 spins keys for auxiliary booster keys
ARRAY_TABLE_INIT(HashTableKeyBooster, 6, int);

// 10 spins keys with double values for eprl 4-simplices
ARRAY_TABLE_INIT(HashTableEPRL, 10, double);

// 3 spins keys with double values for eprl 3-simplices
ARRAY_TABLE_INIT(HashTableEPRL3, 3, double);

// 3 spins keys with double values for 3-booster keys
ARRAY_TABLE_INIT(HashTableBooster3, 3, double);


/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_HASHING_H__*/
