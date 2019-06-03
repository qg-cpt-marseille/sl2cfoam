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

#include "coherentstates.h"
#include "utilities.h"
#include "wigxjpf.h"

double complex CoherentStateI(dspin two_i, dspin two_j1, dspin two_j2,dspin two_j3,dspin two_j4,
                              double f1, double t1, double f2, double t2, double f3, double t3, double f4, double t4,
                              int sign1, int sign2, int sign3, int sign4 ){

    double complex wignerD1;
    double complex wignerD2;
    double complex wignerD3;
    double complex wignerD4;

    double complex coherentStatei;
    double J4symbol;

    coherentStatei = 0 + 0*I ;

    for (int two_m1 = -(int)two_j1; two_m1 <= (int)two_j1; two_m1 += 2) {
      for (int two_m2 = max(-(int)two_i-(int)two_m1,-(int)two_j2); two_m2 <= min((int)two_i-(int)two_m1, (int)two_j2); two_m2 += 2) {
        for (int two_m3 = max(-(int)two_j4-(int)two_m1-(int)two_m2,-(int)two_j3);
                 two_m3 <= min ((int)two_j4-(int)two_m1-(int)two_m2, (int)two_j3); two_m3 += 2) {

          wignerD1 = WignerD(two_j1, two_m1, f1, t1, sign1);
          wignerD2 = WignerD(two_j2, two_m2, f2, t2, sign2);
          wignerD3 = WignerD(two_j3, two_m3, f3, t3, sign3);
          wignerD4 = WignerD(two_j4, -two_m1-two_m2-two_m3, f4, t4, sign4);

          J4symbol = real_negpow(two_i -(-two_m1 -two_m2))
                   * wig3jj(two_j1, two_j2, two_i, two_m1, two_m2, -two_m1-two_m2)
                   * wig3jj(two_j3, two_j4, two_i, two_m3, -two_m1-two_m2-two_m3, two_m1+two_m2);

          coherentStatei += J4symbol * wignerD1 * wignerD2 * wignerD3 * wignerD4;

        }
      }
    }

  return(coherentStatei);

}



double complex WignerD(dspin two_j, int two_m, double f1, double t1, int sign) {

    mpz_t rop1;
    mpz_t rop2;
    mpz_t rop3;

    double j = (float)two_j/2;
    double m = (float)two_m/2;

    // Initialize integer variables for gmp
    mpz_init(rop1);
    mpz_init(rop2);
    mpz_init(rop3);

    unsigned long int op1,op2,op3;
    op1 = 2*j; 
    op2 = j+m;
    op3 = j-m;

    // Factorials
    mpz_fac_ui(rop1, op1);
    mpz_fac_ui(rop2, op2);
    mpz_fac_ui(rop3, op3);

    mpz_div(rop1,rop1,rop2);
    mpz_div(rop1,rop1,rop3);

    // inizialize float mpfr variables to account of the square root of factorials
    mpfr_t sqf;
    mpfr_init(sqf);

    // compute sqrt
    mpfr_set_z(sqf, rop1, MPFR_RNDN);
    mpfr_sqrt(sqf, sqf, MPFR_RNDN);

    // Clear mpz variables
    mpz_clear(rop1);
    mpz_clear(rop2);
    mpz_clear(rop3);

    mpfr_t th1;
    mpfr_t c1;
    mpfr_t s1;

    // Theta Angle for the D small. Divided by two.
    t1 = t1/2;

    // set up a mpfr variable th1 with the value of the half angle
    mpfr_init_set_d(th1, t1, MPFR_RNDN);

    mpfr_init(s1);
    mpfr_init(c1);

    // compute cosine and sine of t1
    mpfr_cos(c1, th1, MPFR_RNDN);
    mpfr_sin(s1, th1, MPFR_RNDN);

    if (sign == -1) {
        mpfr_mul_d(s1, s1, -1,  MPFR_RNDN);
    }

    mpfr_t result1;
    mpfr_t result2;

    mpfr_init(result1);
    mpfr_init(result2);

    if (sign == -1) {
        mpfr_pow_ui(result1, c1, j-m, MPFR_RNDN);
    } else {   
        // Exp of cos
        mpfr_pow_ui(result1, c1, j+m, MPFR_RNDN);
    }

    if (sign == -1) {
        mpfr_pow_ui(result2, s1, j+m, MPFR_RNDN);
    } else {
        mpfr_pow_ui(result2, s1, j-m, MPFR_RNDN);
    }

    mpfr_clear(th1);
    mpfr_clear(c1);
    mpfr_clear(s1);

    mpfr_t result3;
    mpfr_init(result3);

    // Create the D small function
    mpfr_mul(result3, result1, result2, MPFR_RNDN);
    mpfr_mul(result3, result3, sqf, MPFR_RNDN);

    double result = mpfr_get_d(result3,  MPFR_RNDN);

    mpfr_clear(sqf);
    mpfr_clear(result1);
    mpfr_clear(result2);
    mpfr_clear(result3);

    double complex wignerD;

    if (sign == -1) {
        wignerD = cexp(-I*m*(f1))*result*cexp(I*j*(-f1));
    } else {
        wignerD = cexp(-I*m*(f1))*result*cexp(-I*j*(-f1));
    }

    return wignerD;

}
