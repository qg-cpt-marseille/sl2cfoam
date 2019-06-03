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

#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

#include "utilities.h"
#include "config.h"
#include "error.h"

int file_exist(char* filename) {
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}

double complex mpfr_complex_mul(double complex coherentA, double complex coherentB){

    mpfr_t  statesMPFR_A_Re, statesMPFR_B_Re, statesMPFR_A_Im, statesMPFR_B_Im,
            statesMPFR_Re1, statesMPFR_Re2, statesMPFR_Im1, statesMPFR_Im2,
            statesMPFR1, statesMPFR2;

    mpfr_init_set_d(statesMPFR_A_Re, creal(coherentA), MPFR_RNDN);
    mpfr_init_set_d(statesMPFR_B_Re, creal(coherentB), MPFR_RNDN);
    mpfr_init_set_d(statesMPFR_A_Im, cimag(coherentA), MPFR_RNDN);
    mpfr_init_set_d(statesMPFR_B_Im, cimag(coherentB), MPFR_RNDN);

    mpfr_inits(statesMPFR_Re1, statesMPFR_Re2, statesMPFR_Im1, statesMPFR_Im2, NULL);

    mpfr_mul(statesMPFR_Re1, statesMPFR_A_Re, statesMPFR_B_Re, MPFR_RNDN);
    mpfr_mul(statesMPFR_Re2, statesMPFR_A_Im, statesMPFR_B_Im, MPFR_RNDN);
    mpfr_mul(statesMPFR_Im1, statesMPFR_A_Im, statesMPFR_B_Re, MPFR_RNDN);
    mpfr_mul(statesMPFR_Im2, statesMPFR_A_Re, statesMPFR_B_Im, MPFR_RNDN);

    mpfr_inits(statesMPFR1, statesMPFR2, NULL);
    mpfr_sub(statesMPFR1, statesMPFR_Re1, statesMPFR_Re2, MPFR_RNDN);
    mpfr_add(statesMPFR2, statesMPFR_Im1, statesMPFR_Im2, MPFR_RNDN);

    double complex result = mpfr_get_d (statesMPFR1, MPFR_RNDN) + I*mpfr_get_d (statesMPFR2, MPFR_RNDN);

    mpfr_clears(statesMPFR_A_Re, statesMPFR_B_Re, statesMPFR_A_Im, statesMPFR_B_Im,
                statesMPFR_Re1,statesMPFR_Re2, statesMPFR_Im1, statesMPFR_Im2,
                statesMPFR1, statesMPFR2, NULL);

    return result;
}

double complex mpfr_complex_div(double complex coherentA, double complex coherentB) {

    mpfr_t  statesMPFR_A_Re, statesMPFR_B_Re, statesMPFR_A_Im, statesMPFR_B_Im,
            statesMPFR_Den1,statesMPFR_Den2,
            statesMPFR_Re1,statesMPFR_Re2, statesMPFR_Im1,statesMPFR_Im2,
            statesMPFR1, statesMPFR2;

    mpfr_rnd_t MPFR_RNDN;

    mpfr_init_set_d(statesMPFR_A_Re, creal(coherentA), MPFR_RNDN);
    mpfr_init_set_d(statesMPFR_A_Im, cimag(coherentA), MPFR_RNDN);

    mpfr_init_set_d(statesMPFR_B_Re, creal(coherentB), MPFR_RNDN);
    mpfr_init_set_d(statesMPFR_B_Im, cimag(coherentB), MPFR_RNDN);

    mpfr_inits(statesMPFR_Re1, statesMPFR_Re2, statesMPFR_Im1, statesMPFR_Im2,
                         statesMPFR_Den1,statesMPFR_Den2, NULL);

    mpfr_mul(statesMPFR_Re1, statesMPFR_A_Re, statesMPFR_B_Re, MPFR_RNDN);
    mpfr_mul(statesMPFR_Re2, statesMPFR_A_Im, statesMPFR_B_Im, MPFR_RNDN);

    mpfr_mul(statesMPFR_Im1, statesMPFR_A_Im, statesMPFR_B_Re, MPFR_RNDN);
    mpfr_mul(statesMPFR_Im2, statesMPFR_A_Re, statesMPFR_B_Im, MPFR_RNDN);

    mpfr_pow_ui(statesMPFR_Den1, statesMPFR_B_Re, 2, MPFR_RNDN);
    mpfr_pow_ui(statesMPFR_Den2, statesMPFR_B_Im, 2, MPFR_RNDN);
    mpfr_add(statesMPFR_Den1, statesMPFR_Den1, statesMPFR_Den2, MPFR_RNDN );

    mpfr_div(statesMPFR_Re1, statesMPFR_Re1, statesMPFR_Den1, MPFR_RNDN);
    mpfr_div(statesMPFR_Re2, statesMPFR_Re2, statesMPFR_Den1, MPFR_RNDN);
    mpfr_div(statesMPFR_Im1, statesMPFR_Im1, statesMPFR_Den1, MPFR_RNDN);
    mpfr_div(statesMPFR_Im2, statesMPFR_Im2, statesMPFR_Den1, MPFR_RNDN);

    mpfr_inits(statesMPFR1,statesMPFR2, NULL);
    mpfr_add(statesMPFR1, statesMPFR_Re1, statesMPFR_Re2, MPFR_RNDN);
    mpfr_sub(statesMPFR2, statesMPFR_Im1, statesMPFR_Im2, MPFR_RNDN);

    double complex result = mpfr_get_d (statesMPFR1, MPFR_RNDN) + I*mpfr_get_d (statesMPFR2, MPFR_RNDN);

    mpfr_clears(statesMPFR_A_Re, statesMPFR_B_Re, statesMPFR_A_Im, statesMPFR_B_Im,
                statesMPFR_Re1, statesMPFR_Re2, statesMPFR_Im1, statesMPFR_Im2,
                statesMPFR1, statesMPFR2, NULL );

    return result;
}

void compsum_mpfr(double* err, double* sum, mpfr_t toadd) {

    double c, y, t;
    c = mpfr_get_d (toadd, MPFR_RNDN) ;
    y = *sum + c;
    t = (y - *sum) - c;
    *err = t;
    *sum = y;

}

void compsum_mpfr_complex(double complex* err, double complex* sum,
                          mpfr_t toadd_Re, mpfr_t toadd_Im) {

    double complex c , y, t;

    c = mpfr_get_d(toadd_Re, MPFR_RNDN) + mpfr_get_d(toadd_Im, MPFR_RNDN)*I;
    y = *sum + c;
    t = (y - *sum) - c;
    *err = t;
    *sum =  y;

}

void check_data_4simplex(float Immirzi){

    struct stat st;

    char* root_folder = get_data_root();

    if (mkdir(root_folder, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create root directory"); }
    }

    struct data_folders* fd = get_data_folders(Immirzi);

    // make all the directories, or do nothing
    // if they already exist
    // exit if any other error

    if (mkdir(fd->foursimp, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->foursimp_bf, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->foursimp_bf_hashtableampl, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->foursimp_bf_hashtableampl_ampl, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (Immirzi != 0) {

        if (mkdir(fd->foursimp_imm, 0755) == -1) {
            if (errno != EEXIST) { error("cannot create directory"); }
        }

        if (mkdir(fd->foursimp_imm_boost, 0755) == -1) {
            if (errno != EEXIST) { error("cannot create directory"); }
        }

        if (mkdir(fd->foursimp_imm_hashtableampl, 0755) == -1) {
            if (errno != EEXIST) { error("cannot create directory"); }
        }

        if (mkdir(fd->foursimp_imm_hashtableampl_ampl, 0755) == -1) {
            if (errno != EEXIST) { error("cannot create directory"); }
        }

    }


    if (mkdir(fd->foursimp_hashtable6j, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->foursimp_hashtable6j_6j, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->foursimp_aux, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (Immirzi != 0) {

        if (mkdir(fd->foursimp_aux_imm, 0755) == -1) {
            if (errno != EEXIST) { error("cannot create directory"); }
        }

        if (mkdir(fd->foursimp_aux_imm_boost, 0755) == -1) {
            if (errno != EEXIST) { error("cannot create directory"); }
        }

    }

    if (mkdir(fd->foursimp_aux_hashtable15j, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->foursimp_aux_hashtable15j_15j, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

}

void check_data_3simplex(float Immirzi){

    struct stat st;

    char* root_folder = get_data_root();

    if (mkdir(root_folder, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create root directory"); }
    }

    struct data_folders* fd = get_data_folders(Immirzi);

    // make all the directories, or do nothing
    // if they already exist
    // exit if any other error

    if (mkdir(fd->threesimp, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->threesimp_hashtable6j, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->threesimp_imm, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->threesimp_imm_ampl, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(fd->threesimp_imm_boost, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

}



