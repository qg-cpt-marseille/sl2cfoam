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

#include <string.h>

#include "common.h"
#include "config.h"
#include "sl2cfoam.h"

// Runtime internal configuration variables.
static char data_root[1024];
static int  coupling;
static int  verbosity;

// Cached data folder structure
static struct data_folders df;
static float df_Immirzi;


void set_data_root(char* df) {
    strcpy(data_root, df);
}

char* get_data_root() {
    return data_root;
}

struct data_folders* get_data_folders(float Immirzi) {

    // search for cached
    // Immirzi stored with 2-digits precision
    // exclude BF case (Immirzi = 0)
    if (fabs(Immirzi) > 0 && fabs(Immirzi - df_Immirzi) < 1e-3) {
        return &df;
    }

    // build and cache for this Immirzi
    df_Immirzi = Immirzi;
    build_data_folders(&df, Immirzi);

    return &df;

}

void build_data_folders(struct data_folders* df, float Immirzi) {

    df->foursimp = (char*)malloc(1024*sizeof(char));
    df->foursimp_bf = (char*)malloc(1024*sizeof(char));
    df->foursimp_bf_hashtableampl = (char*)malloc(1024*sizeof(char));
    df->foursimp_bf_hashtableampl_ampl = (char*)malloc(1024*sizeof(char));
    df->foursimp_imm = (char*)malloc(1024*sizeof(char));
    df->foursimp_imm_boost = (char*)malloc(1024*sizeof(char));
    df->foursimp_imm_hashtableampl = (char*)malloc(1024*sizeof(char));
    df->foursimp_imm_hashtableampl_ampl = (char*)malloc(1024*sizeof(char));
    df->foursimp_hashtable6j = (char*)malloc(1024*sizeof(char));
    df->foursimp_hashtable6j_6j = (char*)malloc(1024*sizeof(char));

    df->foursimp_aux = (char*)malloc(1024*sizeof(char));
    df->foursimp_aux_imm = (char*)malloc(1024*sizeof(char));
    df->foursimp_aux_imm_boost = (char*)malloc(1024*sizeof(char));
    df->foursimp_aux_hashtable15j = (char*)malloc(1024*sizeof(char));
    df->foursimp_aux_hashtable15j_15j = (char*)malloc(1024*sizeof(char));

    df->threesimp = (char*)malloc(1024*sizeof(char));
    df->threesimp_imm = (char*)malloc(1024*sizeof(char));
    df->threesimp_imm_boost = (char*)malloc(1024*sizeof(char));
    df->threesimp_imm_ampl = (char*)malloc(1024*sizeof(char));
    df->threesimp_hashtable6j = (char*)malloc(1024*sizeof(char));
    
    strcpy(df->foursimp, data_root);
    strcpy(df->foursimp_bf, data_root);
    strcpy(df->foursimp_bf_hashtableampl, data_root);
    strcpy(df->foursimp_bf_hashtableampl_ampl, data_root);
    strcpy(df->foursimp_imm, data_root);
    strcpy(df->foursimp_imm_boost, data_root);
    strcpy(df->foursimp_imm_hashtableampl, data_root);
    strcpy(df->foursimp_imm_hashtableampl_ampl, data_root);
    strcpy(df->foursimp_hashtable6j, data_root);
    strcpy(df->foursimp_hashtable6j_6j, data_root);
    strcpy(df->foursimp_aux, data_root);
    strcpy(df->foursimp_aux_imm, data_root);
    strcpy(df->foursimp_aux_imm_boost, data_root);
    strcpy(df->foursimp_aux_hashtable15j, data_root);
    strcpy(df->foursimp_aux_hashtable15j_15j, data_root);

    strcpy(df->threesimp, data_root);
    strcpy(df->threesimp_imm, data_root);
    strcpy(df->threesimp_imm_boost, data_root);
    strcpy(df->threesimp_imm_ampl, data_root);
    strcpy(df->threesimp_hashtable6j, data_root);

    char tmp[1024];

    // 4 simplex folders

    strcat(df->foursimp, "/4simplex/");

    strcat(df->foursimp_bf, "/4simplex/bf/");

    strcat(df->foursimp_bf_hashtableampl, "/4simplex/bf/hashtablesBF/");

    sprintf(tmp, "/4simplex/immirzi_%.2f/", Immirzi);
    strcat(df->foursimp_imm, tmp);

    sprintf(tmp, "/4simplex/immirzi_%.2f/hashtablesBooster/", Immirzi);
    strcat(df->foursimp_imm_boost, tmp);

    sprintf(tmp, "/4simplex/immirzi_%.2f/hashtablesEPRL/", Immirzi);
    strcat(df->foursimp_imm_hashtableampl, tmp);

    strcat(df->foursimp_aux, "/4simplex/aux/");

    sprintf(tmp, "/4simplex/aux/immirzi_%.2f/", Immirzi);
    strcat(df->foursimp_aux_imm, tmp);

    sprintf(tmp, "/4simplex/aux/immirzi_%.2f/hashtablesKeyBooster/", Immirzi);
    strcat(df->foursimp_aux_imm_boost, tmp);

    strcat(df->foursimp_hashtable6j, "/4simplex/hashtablesJ6/");

    strcat(df->foursimp_aux_hashtable15j, "/4simplex/aux/hashtablesKeyJ15/");

    if (coupling == SL2CFOAM_COUPLING_REDUCIBLE) {

        strcat(df->foursimp_hashtable6j_6j, "/4simplex/hashtablesJ6/reducibleJ6/");

        strcat(df->foursimp_bf_hashtableampl_ampl, "/4simplex/bf/hashtablesBF/reducibleBF/");

        sprintf(tmp, "/4simplex/immirzi_%.2f/hashtablesEPRL/reducibleEPRL/", Immirzi);
        strcat(df->foursimp_imm_hashtableampl_ampl, tmp);

        strcat(df->foursimp_aux_hashtable15j_15j, "/4simplex/aux/hashtablesKeyJ15/reducibleKeyJ15/");

    } else {

        strcat(df->foursimp_hashtable6j_6j, "/4simplex/hashtablesJ6/irreducibleJ6/");

        strcat(df->foursimp_bf_hashtableampl_ampl, "/4simplex/bf/hashtablesBF/irreducibleBF/");

        sprintf(tmp, "/4simplex/immirzi_%.2f/hashtablesEPRL/irreducibleEPRL/", Immirzi);
        strcat(df->foursimp_imm_hashtableampl_ampl, tmp);

        strcat(df->foursimp_aux_hashtable15j_15j, "/4simplex/aux/hashtablesKeyJ15/irreducibleKeyJ15/");

    }

    // 3 simplex folders

    strcat(df->threesimp, "/3simplex/");

    strcat(df->threesimp_hashtable6j, "/3simplex/hashtablesJ6/");

    sprintf(tmp, "/3simplex/immirzi_%.2f/", Immirzi);
    strcat(df->threesimp_imm, tmp);

    sprintf(tmp, "/3simplex/immirzi_%.2f/hashtablesBooster/", Immirzi);
    strcat(df->threesimp_imm_boost, tmp);

    sprintf(tmp, "/3simplex/immirzi_%.2f/hashtablesEPRL/", Immirzi);
    strcat(df->threesimp_imm_ampl, tmp);

}

void free_data_folders(struct data_folders* df) {

    free(df->foursimp);
    free(df->foursimp_hashtable6j);
    free(df->foursimp_hashtable6j_6j);
    free(df->foursimp_imm);
    free(df->foursimp_imm_hashtableampl);
    free(df->foursimp_imm_hashtableampl_ampl);
    free(df->foursimp_bf);
    free(df->foursimp_bf_hashtableampl);
    free(df->foursimp_bf_hashtableampl_ampl);
    free(df->foursimp_imm_boost);
    free(df->foursimp_aux);
    free(df->foursimp_aux_hashtable15j);
    free(df->foursimp_aux_hashtable15j_15j);
    free(df->foursimp_aux_imm);
    free(df->foursimp_aux_imm_boost);
    free(df->threesimp);
    free(df->threesimp_hashtable6j);
    free(df->threesimp_imm);
    free(df->threesimp_imm_boost);
    free(df->threesimp_imm_ampl);

}


// Get/set coupling type.
void set_coupling(int c) {
    coupling = c;
}

int get_coupling() {
    return coupling;
}

// Get/set verbosity.
void set_verbosity(int v) {
    verbosity = v;
}

int get_verbosity() {
    return verbosity;
}
