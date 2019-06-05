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

#ifndef __SL2CFOAM_LIB_H__
#define __SL2CFOAM_LIB_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

// Machine type for double-spins.
typedef unsigned short dspin;

// Machine type for half-integer spins.
typedef float spin;

////////////////////////////////////////////////////////////////////////
// Configuration options
////////////////////////////////////////////////////////////////////////

// Coupling types.
#define SL2CFOAM_COUPLING_REDUCIBLE   0
#define SL2CFOAM_COUPLING_IRREDUCIBLE 1

// Verbosity levels.
#define SL2CFOAM_VERBOSE_OFF   0
#define SL2CFOAM_VERBOSE_LOW   1
#define SL2CFOAM_VERBOSE_HIGH  2

// Contains general parameters for setup of the library.
struct sl2cfoam_config {
    char* data_folder;
    int   coupling;
    int   verbosity;
    int   store_ampls;
};

////////////////////////////////////////////////////////////////////////
// Global setup functions.
////////////////////////////////////////////////////////////////////////

// Call this function to initialize the library
// at the start of the program.
void sl2cfoam_init_conf(struct sl2cfoam_config* conf);

// Initalize with default configurations.
void sl2cfoam_init();

// Call this function when finished using the library
// at the end of the program.
void sl2cfoam_free();

// If writing multi-threaded programs, call this
// function when a thread starts.
void sl2cfoam_init_thread();

// If writing multi-threaded programs, call this
// function when a thread ends.
void sl2cfoam_free_thread();

///////////////////////////////////////////////////////////////////////////
// FOUR-simplex functions.
///////////////////////////////////////////////////////////////////////////

// Self-contained computation of the vertex amplitude of a four simplex.
// WARNING: not safe to be parallelized.
double sl2cfoam_four_simplex(dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5,
                             dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10,
                             dspin two_i1, dspin two_i2, dspin two_i3, dspin two_i4, dspin two_i5,
                             dspin two_Dl, float Immirzi);

// Precomputation of 6js and boosters function for
// a given intertwiner i4.
void sl2cfoam_hash_symbols(dspin two_js[10],
                           dspin two_i4_min, dspin two_i4_max,
                           dspin two_Dl, float Immirzi);

// Computes all the amplitudes in the given intertwiner range
// and stores the results on disk. Parallelized over i4.
void sl2cfoam_hash_four_ampl(dspin two_js[10],
                             dspin two_i1_min, dspin two_i1_max, 
                             dspin two_i2_min, dspin two_i2_max, 
                             dspin two_i3_min, dspin two_i3_max, 
                             dspin two_i4_min, dspin two_i4_max, 
                             dspin two_i5_min, dspin two_i5_max, 
                             dspin two_Dl, float Immirzi);

// Loads a precomputed amplitudes table from disk
// given i4, two Delta-l and Immirzi.
void sl2cfoam_four_ampl_load(dspin two_js[10], dspin i4,
                             dspin two_Dl, float Immirzi);

// Unloads the previously loaded amplitudes table.
void sl2cfoam_four_ampl_unload();

// Takes the value of an amplitude from the loaded table.
double sl2cfoam_four_ampl_get(dspin two_js[10],
                              dspin two_i1, dspin two_i2, dspin two_i3,
                              dspin two_i4, dspin two_i5,
                              dspin two_Dl, float Immirzi);


///////////////////////////////////////////////////////////////////////////
// FOUR-simplex functions, BF theory.
///////////////////////////////////////////////////////////////////////////

double sl2cfoam_four_simplex_BF(dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5,
                                dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10,
                                dspin two_i1, dspin two_i2, dspin two_i3, dspin two_i4, dspin two_i5);

void sl2cfoam_hash_symbols_BF(dspin two_js[10],
                              dspin two_i4_min, dspin two_i4_max);

void sl2cfoam_hash_four_ampl_BF(dspin two_js[10],
                                dspin two_i1_min, dspin two_i1_max, 
                                dspin two_i2_min, dspin two_i2_max, 
                                dspin two_i3_min, dspin two_i3_max, 
                                dspin two_i4_min, dspin two_i4_max, 
                                dspin two_i5_min, dspin two_i5_max);

void sl2cfoam_four_ampl_load_BF(dspin two_js[10], dspin i4);

double sl2cfoam_four_ampl_get_BF(dspin two_js[10],
                                dspin two_i1, dspin two_i2, dspin two_i3,
                                dspin two_i4, dspin two_i5);


///////////////////////////////////////////////////////////////////////////
// THREE-simplex functions.
///////////////////////////////////////////////////////////////////////////

double sl2cfoam_three_simplex( dspin two_j1, dspin two_j2, dspin two_j3, 
                               dspin two_j4, dspin two_j5, dspin two_j6,
                               dspin two_Dl, float Immirzi);

void sl2cfoam_hash_three_ampl(dspin two_j1_min, dspin two_j1_max,
                              dspin two_j2_min, dspin two_j2_max,
                              dspin two_j3_min, dspin two_j3_max,
                              dspin two_j4_min, dspin two_j4_max,
                              dspin two_j5_min, dspin two_j5_max,
                              dspin two_j6_min, dspin two_j6_max,                                                                                                        
                              dspin two_Dl, float Immirzi);
                              
void sl2cfoam_hashall_three_ampl(dspin two_minj, dspin two_maxj,
                                 dspin two_Dl, float Immirzi); 

void sl2cfoam_three_ampl_load(dspin two_j1, dspin two_j2, dspin two_j3,
                              dspin two_Dl, float Immirzi);        

void sl2cfoam_three_ampl_unload();  

double sl2cfoam_three_ampl_get(dspin two_j1, dspin two_j2, dspin two_j3,
                               dspin two_j4, dspin two_j5, dspin two_j6,
                               dspin two_Dl, float Immirzi);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_LIB_H__*/
