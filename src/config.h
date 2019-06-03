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


#ifndef __SL2CFOAM_CONFIG_H__
#define __SL2CFOAM_CONFIG_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

////////////////////////////////////////////////////////////////////////
// Global compile-time parameters.
////////////////////////////////////////////////////////////////////////

// Compile or not debug code.
#define DEBUG_ON 0

// Set high bits precision for MPFR/MPC.
// Low default precision set in lib initialization.
#define MPBITS 1256

////////////////////////////////////////////////////////////////////////
// Setters and getters for runtime configuration.
////////////////////////////////////////////////////////////////////////

// Contains the folder structure for a given Immirzi.
struct data_folders {

    char* foursimp;
    char* foursimp_bf;
    char* foursimp_bf_hashtableampl;
    char* foursimp_bf_hashtableampl_ampl;
    char* foursimp_imm;
    char* foursimp_imm_boost;
    char* foursimp_imm_hashtableampl;
    char* foursimp_imm_hashtableampl_ampl;
    char* foursimp_hashtable6j;
    char* foursimp_hashtable6j_6j;
    char* foursimp_aux;
    char* foursimp_aux_imm;
    char* foursimp_aux_imm_boost;
    char* foursimp_aux_hashtable15j;
    char* foursimp_aux_hashtable15j_15j;

    char* threesimp;
    char* threesimp_hashtable6j;
    char* threesimp_imm;
    char* threesimp_imm_boost;
    char* threesimp_imm_ampl;

};

// Get/set main data folder.
void set_data_root(char* root_folder);
char* get_data_root();

// Returns the folder structure for a given Immirzi.
// Cached version.
struct data_folders* get_data_folders(float Immirzi);

// Constructs all the paths for subfolders.
void build_data_folders(struct data_folders* df, float Immirzi);

// Call to free memory of subfolders paths.
void free_data_folders(struct data_folders* df);

// Get/set coupling type.
void set_coupling(int coupling);
int get_coupling();

// Get/set verbosity.
void set_verbosity(int verbosity);
int get_verbosity();


/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_CONFIG_H__*/
