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

#include "recouplingSL2C.h"
#include "cgamma.h"
#include "config.h"
#include "utilities.h"
#include "wigxjpf.h"
#include "error.h"

// Internal function to compute intermediate coefficient k.
static void kappa(mpc_t rop, int two_j1, int two_j2, int two_j3,
                  int two_k1, int two_k2, int two_k3,
									double two_rho1, double two_rho2, double two_rho3) {

  // check for integers
  ensure_integer_spin(two_j1+two_j2-two_j3);
  ensure_integer_spin(two_j2+two_j3-two_j1);
  ensure_integer_spin(two_j1+two_j3-two_j2);
  ensure_integer_spin(two_j1+two_j2+two_j3);

	double j1 = ((double)two_j1)/2.0;
 	double j2 = ((double)two_j2)/2.0;
 	double j3 = ((double)two_j3)/2.0;

 	double k1 = ((double)two_k1)/2.0;
 	double k2 = ((double)two_k2)/2.0;
 	double k3 = ((double)two_k3)/2.0;
 	double K = k1+k2+k3;

 	double rho_1 = two_rho1/2.0;
 	double rho_2 = two_rho2/2.0;
 	double rho_3 = two_rho3/2.0;

  mpc_t result;
  mpc_init2(result, MPBITS);

  // real ap vars
  mpfr_t x, y;
  mpfr_init2(x, MPBITS);
  mpfr_init2(y, MPBITS);

  // complex ap vars
  mpc_t z, w;
  mpc_init2(z, MPBITS);
  mpc_init2(w, MPBITS);

  // first factors in the product, before the sums
	double complex phase = complex_negpow(-two_k1-two_k2+two_j1-two_j2+two_j3);
  mpc_set_dc(result, phase, MPC_RNDNN);

  mpfr_set_ui(x, (unsigned long int)(d(two_j3)), MPFR_RNDN);
  mpfr_sqrt(y, x, MPFR_RNDN);
  mpc_div_fr(result, result, y, MPC_RNDNN);

  // phase part
  mpc_set_d_d(z, 1.0 + j3, rho_3, MPC_RNDNN);
  complex_lngamma(w, z);

  mpc_set_d_d(z, 1.0 + j1, -rho_1, MPC_RNDNN);
  complex_lngamma(z, z);
  mpc_add(w, w, z, MPC_RNDNN);

  mpc_set_d_d(z, 1.0 + j2, -rho_2, MPC_RNDNN);
  complex_lngamma(z, z);
  mpc_add(w, w, z, MPC_RNDNN);

  mpc_exp(w, w, MPC_RNDNN);

  mpc_abs(x, w, MPC_RNDNN);
  mpc_div_fr(w, w, x, MPC_RNDNN);
	//printf ("G %i %i %i %17g %17g \n",two_j1,two_j2,two_j3,mpc_get_dc(w,MPC_RNDNN));
  mpc_mul(result, result, w, MPC_RNDNN);

  // check the factorial are non-negative
  assert((two_j1 - two_k1) >= 0);
  assert((two_j2 + two_k2) >= 0);
  assert((two_j1 + two_k1) >= 0);
  assert((two_j2 - two_k2) >= 0);

  // sqrt of factorial before the sums
  mpfr_fac_ui(x, (dspin)(two_j1 - two_k1)/2, MPFR_RNDN);
  mpfr_fac_ui(y, (dspin)(two_j2 + two_k2)/2, MPFR_RNDN);
  mpfr_mul(x, x, y, MPFR_RNDN);
  mpfr_fac_ui(y, (dspin)(two_j1 + two_k1)/2, MPFR_RNDN);
  mpfr_div(x, x, y, MPFR_RNDN);
  mpfr_fac_ui(y, (dspin)(two_j2 - two_k2)/2, MPFR_RNDN);
  mpfr_div(x, x, y, MPFR_RNDN);
  mpfr_sqrt(x, x, MPFR_RNDN);
  mpc_mul_fr(result, result, x, MPC_RNDNN);

  // now the sums remain ...
  mpc_t nsum, ssum, s1prod, s2prod;
  mpc_init2(nsum, MPBITS);
  mpc_init2(ssum, MPBITS);
  mpc_init2(s1prod, MPBITS);
  mpc_init2(s2prod, MPBITS);
  mpc_set_ui(nsum, 0, MPC_RNDNN);

  mpfr_t spref;
  mpfr_init2(spref, MPBITS);

  double wig;
  int two_n, two_s1, two_s2;
  float s1, s2;

  for (two_n = -two_j1; two_n <= min(two_j1, two_k3 + two_j2); two_n += 2) {

    // check for allowed summands
    if (two_j1-two_n < 0) {
        continue;
    }
    if (two_j2+two_k3-two_n < 0) {
        continue;
    }
    if (two_j1+two_n < 0) {
        continue;
    }
    if (two_j2-two_k3+two_n < 0) {
        continue;
    }

    // sqrt of factorial into sum over n

    mpfr_fac_ui(x, (dspin)(two_j1 + two_n)/2, MPFR_RNDN);
    mpfr_fac_ui(y, (dspin)(two_j2 - two_k3 + two_n)/2, MPFR_RNDN);
    mpfr_mul(x, x, y, MPFR_RNDN);

    mpfr_fac_ui(y, (unsigned long int)(two_j1 - two_n)/2, MPFR_RNDN);
    mpfr_div(x, y, x, MPFR_RNDN);
    mpfr_fac_ui(y, (dspin)(two_j2 + two_k3 - two_n)/2, MPFR_RNDN);
    mpfr_mul(x, x, y, MPFR_RNDN);

    mpfr_sqrt(spref, x, MPFR_RNDN);

    wig = wig3jj(	two_j1, two_j2, two_j3,
        					two_n,  two_k3-two_n, -two_k3);

    mpfr_mul_d(spref, spref, wig, MPFR_RNDN);
    mpc_set_ui(ssum, 0, MPC_RNDNN);

    for (two_s1 = max(two_k1, two_n); two_s1 <= two_j1; two_s1 += 2) {

      s1 = ((double)two_s1) / 2.0;

      // check for allowed summands
      if (two_j1+two_s1 < 0) {
          continue;
      }
      if (two_j1-two_s1 < 0) {
          continue;
      }
      if (two_s1-two_k1 < 0) {
          continue;
      }
      if (two_s1-two_n < 0) {
          continue;
      }

      mpfr_fac_ui(x, (dspin)(two_j1-two_s1)/2, MPFR_RNDN);
      mpfr_fac_ui(y, (dspin)(two_s1-two_k1)/2, MPFR_RNDN);
      mpfr_mul(x, x, y, MPFR_RNDN);
      mpfr_fac_ui(y, (dspin)(two_s1-two_n)/2, MPFR_RNDN);
      mpfr_mul(x, x, y, MPFR_RNDN);

      mpfr_fac_ui(y, (dspin)(two_j1+two_s1)/2, MPFR_RNDN);
      mpfr_div(x, y, x, MPFR_RNDN);

      mpc_set_fr(s1prod, x, MPC_RNDNN);

      mpc_set_d_d(z, (1.0-K+two_s1)/2.0, -(rho_1-rho_2-rho_3)/2.0, MPC_RNDNN);
      complex_lngamma(w, z);

      mpc_set_d_d(z, 1.0+s1, -rho_1, MPC_RNDNN);
      complex_lngamma(z, z);
      mpc_sub(w, w, z, MPC_RNDNN);

      mpc_exp(w, w, MPC_RNDNN);
      mpc_mul(s1prod, s1prod, w, MPC_RNDNN);

      for (two_s2 = max(-two_k2, two_n-two_k3); two_s2 <= two_j2; two_s2 += 2) {

        s2 = ((double)two_s2) / 2.0;

        // check for allowed summands
        if (two_j2+two_s2 < 0) {
            continue;
        }
        if (two_j2-two_s2 < 0) {
            continue;
        }
        if (two_k2+two_s2 < 0) {
            continue;
        }
        if (two_k3-two_n+two_s2 < 0) {
            continue;
        }

        mpc_set(s2prod, s1prod, MPC_RNDNN);
				double complex phase_s = complex_negpow(two_s1+two_s2-two_k1+two_k2);
        mpc_set_dc(z, phase_s , MPC_RNDNN);
        mpc_mul(s2prod, s2prod, z, MPC_RNDNN);

        mpfr_fac_ui(x, (dspin)(two_j2-two_s2)/2, MPFR_RNDN);
        mpfr_fac_ui(y, (dspin)(two_s2+two_k2)/2, MPFR_RNDN);
        mpfr_mul(x, x, y, MPFR_RNDN);
        mpfr_fac_ui(y, (dspin)(two_k3-two_n+two_s2)/2, MPFR_RNDN);
        mpfr_mul(x, x, y, MPFR_RNDN);

        mpfr_fac_ui(y, (dspin)(two_j2+two_s2)/2, MPFR_RNDN);
        mpfr_div(x, y, x, MPFR_RNDN);

        mpc_mul_fr(s2prod, s2prod, x, MPC_RNDNN);

        mpc_set_d_d(z, (1.0+K+two_s2)/2.0, +(rho_1-rho_2+rho_3)/2.0, MPC_RNDNN);
        complex_lngamma(w, z);

        mpc_set_d_d(z, (1.0-k1+k2+k3-two_n+two_s1+two_s2)/2.0, -(rho_1+rho_2-rho_3)/2.0, MPC_RNDNN);
        complex_lngamma(z, z);
        mpc_add(w, w, z, MPC_RNDNN);

        mpc_set_d_d(z, 1.0+s2, -rho_2, MPC_RNDNN);
        complex_lngamma(z, z);
        mpc_sub(w, w, z, MPC_RNDNN);

        mpc_set_d_d(z, 1.0+s1+s2, rho_3, MPC_RNDNN);
        complex_lngamma(z, z);
        mpc_sub(w, w, z, MPC_RNDNN);

        mpc_exp(w, w, MPC_RNDNN);
        mpc_mul(s2prod, s2prod, w, MPC_RNDNN);

        mpc_add(ssum, ssum, s2prod, MPC_RNDNN);
      }
    }

    mpc_set_d_d(z, (1.0-k1+k2+k3-two_n)/2.0, -(rho_1+rho_2+rho_3)/2.0, MPC_RNDNN);
    complex_lngamma(w, z);
    mpc_exp(w, w, MPC_RNDNN);
    mpc_div(ssum, ssum, w, MPC_RNDNN);

    mpc_mul_fr(ssum, ssum, spref, MPC_RNDNN);
    mpc_add(nsum, nsum, ssum, MPC_RNDNN);

  }

  mpc_mul(result, result, nsum, MPC_RNDNN);
	mpc_set(rop, result, MPC_RNDNN);

  // clear variables
  mpfr_clear(x);
  mpfr_clear(y);
  mpfr_clear(spref);
  mpc_clear(result);
  mpc_clear(z);
  mpc_clear(w);
  mpc_clear(nsum);
  mpc_clear(ssum);
  mpc_clear(s1prod);
  mpc_clear(s2prod);

}

double complex chi(dspin two_j1, dspin two_j2, dspin two_j3,
                   dspin two_k1, dspin two_k2, dspin two_k3,
                   double two_rho1, double two_rho2, double two_rho3) {

  ensure_integer_spin(two_j1 + two_j2 + two_j3 + two_k1 + two_k2 + two_k3);

  spin j1, j2, j3, k1, k2, k3, rho1, rho2, rho3;

  k1 = ((spin)(two_k1))/2.0; k2 = ((spin)(two_k2))/2.0; k3 = ((spin)(two_k3))/2.0;
  j1 = ((spin)(two_j1))/2.0; j2 = ((spin)(two_j2))/2.0; j3 = ((spin)(two_j3))/2.0;
  rho1 = two_rho1/2.0; rho2 = two_rho2/2.0; rho3 = two_rho3/2.0;

  double complex phase;
  phase = complex_negpow((int)(j1 + j2 + j3 + k1 + k2 + k3)) / (4 * sqrt(2*M_PI));

	mpfr_t x;
	mpfr_init2(x, MPBITS);

	mpc_t z,w,y;
	mpc_init2(z, MPBITS);
	mpc_init2(w, MPBITS);
	mpc_init2(y, MPBITS);

	mpc_set_d_d(z, (1-k1-k2-k3)/2, (rho1+rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(w, z);

	mpc_set_d_d(z, (1-k1+k2+k3)/2, (rho1-rho2-rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

	mpc_set_d_d(z, (1+k1-k2+k3)/2, -(rho1-rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

	mpc_set_d_d(z, (1-k1-k2+k3)/2, (rho1+rho2-rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

	mpc_exp(w, w, MPC_RNDNN);

	mpc_abs(x, w, MPC_RNDNN);
	mpc_div_fr(w, w, x, MPC_RNDNN);

	mpc_set_d_d(z, (1+k1+k2+k3)/2, -(rho1+rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(y, z);

	mpc_set_d_d(z, (1-k1-k2-k3)/2, -(rho1+rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(y, y, z, MPC_RNDNN);

	mpc_exp(y, y, MPC_RNDNN);
	mpc_mul(z, y, w, MPC_RNDNN );

	mpc_t kappa_res;
	mpc_init2(kappa_res, MPBITS);

	kappa (kappa_res, two_j1, two_j2, two_j3,
				 two_k1, two_k2, two_k3,
				 two_rho1, two_rho2, two_rho3);

	mpc_mul (z,z,kappa_res,MPC_RNDNN );

	double complex result = phase * sqrt(d(two_j1)*d(two_j2)*d(two_j3))
													      * mpc_get_dc(z, MPC_RNDNN);

	mpc_clear(z);
	mpc_clear(y);
	mpc_clear(w);
	mpc_clear(kappa_res);
	mpfr_clear(x);

  return result;

}

double complex chi_simp(dspin two_l1, dspin two_l2, dspin two_l3,
                        dspin two_j1, dspin two_j2, dspin two_j3,
                        float Immirzi) {

  assert(two_l1 >= two_j1);
  assert(two_l2 >= two_j2);
  assert(two_l3 >= two_j3);                          

  double complex result =  chi(two_l1, two_l2, two_l3,
                               two_j1, two_j2, two_j3,
                               Immirzi*two_j1, Immirzi*two_j2, Immirzi*two_j3);

  return result;

}

void chi_mpc( mpfr_t *chi_value, dspin two_j1, dspin two_j2, dspin two_j3,
                                 dspin two_k1, dspin two_k2, dspin two_k3,
                                 double two_rho1, double two_rho2, double two_rho3) {
                                 
  ensure_integer_spin(two_j1 + two_j2 + two_j3 + two_k1 + two_k2 + two_k3);

  spin j1, j2, j3, k1, k2, k3, rho1, rho2, rho3;

  k1 = ((spin)(two_k1))/2.0; k2 = ((spin)(two_k2))/2.0; k3 = ((spin)(two_k3))/2.0;
  j1 = ((spin)(two_j1))/2.0; j2 = ((spin)(two_j2))/2.0; j3 = ((spin)(two_j3))/2.0;
  rho1 = two_rho1/2.0; rho2 = two_rho2/2.0; rho3 = two_rho3/2.0;

  double complex prefactor =  sqrt(d(two_j1)*d(two_j2)*d(two_j3))*complex_negpow((int)(j1 + j2 + j3 + k1 + k2 + k3)) / (4 * sqrt(2*M_PI));

  mpc_t prefactor_mpc;
  mpc_init2(prefactor_mpc, MPBITS);
  mpc_set_dc(prefactor_mpc, prefactor, MPC_RNDNN);

	mpfr_t x;
	mpfr_init2(x, MPBITS);

	mpc_t z,w,y;
	mpc_init2(z, MPBITS);
	mpc_init2(w, MPBITS);
	mpc_init2(y, MPBITS);

	mpc_set_d_d(z, (1-k1-k2-k3)/2, (rho1+rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(w, z);

	mpc_set_d_d(z, (1-k1+k2+k3)/2, (rho1-rho2-rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

	mpc_set_d_d(z, (1+k1-k2+k3)/2, -(rho1-rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

	mpc_set_d_d(z, (1-k1-k2+k3)/2, (rho1+rho2-rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

	mpc_exp(w, w, MPC_RNDNN);

	mpc_abs(x, w, MPC_RNDNN);
	mpc_div_fr(w, w, x, MPC_RNDNN);

	mpc_set_d_d(z, (1+k1+k2+k3)/2, -(rho1+rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(y, z);

	mpc_set_d_d(z, (1-k1-k2-k3)/2, -(rho1+rho2+rho3)/2, MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(y, y, z, MPC_RNDNN);

	mpc_exp(y, y, MPC_RNDNN);
	mpc_mul(z, y, w, MPC_RNDNN );

	mpc_t kappa_res;
	mpc_init2(kappa_res, MPBITS);

	kappa (kappa_res, two_j1, two_j2, two_j3,
				 two_k1, two_k2, two_k3,
				 two_rho1, two_rho2, two_rho3);

	mpc_mul (z,z,kappa_res,MPC_RNDNN );

  mpc_mul(z,z,prefactor_mpc,MPC_RNDNN );

  mpc_real(*chi_value, z, MPC_RNDNN ) ; 

  mpc_clear(prefactor_mpc);
	mpc_clear(z);
	mpc_clear(y);
	mpc_clear(w);
	mpc_clear(kappa_res);
	mpfr_clear(x);


}

// a function to compute a simple chi function involving
// zeros. we need this function in order to solve 
// a problem within recursion relations.
void chi_mpc_simp0( mpfr_t *chi_value, dspin two_j, dspin two_k,
                    float Immirzi) {
                                     
  assert((two_j - two_k) >= 0);

  spin j, k;

  k = ((spin)(two_k))/2.0;
  j = ((spin)(two_j))/2.0;

  double complex prefactordc =  complex_negpow(two_k)/(4*sqrt(2*M_PI));
  
  mpc_t prefactor;
  mpc_init2(prefactor, MPBITS);
  mpc_set_dc(prefactor, prefactordc, MPC_RNDNN);

  mpfr_t prefactor_1, prefactor_2;
  mpfr_init2(prefactor_1, MPBITS);
  mpfr_init2(prefactor_2, MPBITS);

  mpfr_fac_ui(prefactor_2, (dspin)(two_j + two_k)/2, MPFR_RNDN);
  mpfr_set(prefactor_1, prefactor_2, MPFR_RNDN);

  mpfr_fac_ui(prefactor_2, (dspin)(two_j - two_k)/2, MPFR_RNDN);
  mpfr_div(prefactor_1, prefactor_1, prefactor_2, MPFR_RNDN);
  mpc_mul_fr(prefactor, prefactor, prefactor_1, MPFR_RNDN);

  mpfr_t a;
  mpfr_init2(a,MPBITS);
  
	mpc_t z,w,q;
	mpc_init2(z, MPBITS);
	mpc_init2(w, MPBITS);
	mpc_init2(q, MPBITS);

	mpc_set_d_d(z, 1+j, -(Immirzi*k), MPC_RNDNN);
	complex_lngamma(w, z);

	mpc_set_d_d(z, 1+j, (Immirzi*k), MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

  mpc_set(q,w,MPC_RNDNN);

	mpc_set_d_d(z, (1.0-2*k)/2.0, -(Immirzi*k), MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

	mpc_set_d_d(z, (1.0+2*k)/2.0, -(Immirzi*k), MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);
	mpc_add(w, w, z, MPC_RNDNN);

  mpc_add(q,q,z,MPC_RNDNN);

  mpc_set_d_d(z, (1.0-2.0*k)/2.0, (Immirzi*k), MPC_RNDNN);
	complex_lngamma(z, z);
	mpc_add(w, w, z, MPC_RNDNN);

  mpc_add(q,q,z,MPC_RNDNN);

	mpc_exp(w, w, MPC_RNDNN);
  mpc_exp(q,q, MPC_RNDNN);
  mpc_abs(a,q,MPFR_RNDN);

  mpc_div_fr(w,w,a,MPC_RNDNN);

	mpc_t result;

  mpc_init2(result, MPBITS);
  mpc_set(result, prefactor, MPC_RNDNN);
  mpc_mul(result, result, w, MPC_RNDNN);

  // now we have to perform the summation

  int two_s2;
  double s2;

  mpc_t s2prod, s2sum;

  mpc_init2(s2prod, MPBITS);
  mpc_init2(s2sum, MPBITS);
  mpc_set_d_d(s2sum,0,0,MPC_RNDNN);

	mpfr_t x,y;
	mpfr_init2(x, MPBITS);
  mpfr_init2(y, MPBITS);

  for (two_s2 = - two_k; two_s2 <= two_j; two_s2 += 2) {

      s2 = ((double)two_s2) / 2.0;

      mpfr_set_d(prefactor_1, real_negpow(two_k + two_s2), MPFR_RNDN);

      mpfr_fac_ui(x, (dspin)(two_j+two_s2)/2, MPFR_RNDN);
      mpfr_fac_ui(y, (dspin)(two_j-two_s2)/2, MPFR_RNDN);
      mpfr_div(x, x, y, MPFR_RNDN);
      mpfr_fac_ui(y, (dspin)(two_k+two_s2)/2, MPFR_RNDN);
      mpfr_div(x, x, y, MPFR_RNDN);
      mpfr_div(x, x, y, MPFR_RNDN);
      mpfr_mul(x,x,prefactor_1,MPFR_RNDN);

      mpc_set_fr(s2prod, x, MPC_RNDNN);

      mpc_set_d_d(z, (1.0+2.0*k+2.0*s2)/2.0, 0.0, MPC_RNDNN);
      complex_lngamma(w, z);

      mpc_set_d_d(z, (1.0+2.0*k+2.0*s2)/2.0, 0.0, MPC_RNDNN);
      complex_lngamma(z, z);
      mpc_add(w, w, z, MPC_RNDNN);

      mpc_set_d_d(z, (1.0-2.0*k)/2.0, Immirzi*k, MPC_RNDNN);
      complex_lngamma(z, z);
      mpc_add(w, w, z, MPC_RNDNN);

      mpc_set_d_d(z, 1.0+s2, -Immirzi*k, MPC_RNDNN);
      complex_lngamma(z, z);
      mpc_sub(w, w, z, MPC_RNDNN);

      mpc_set_d_d(z, 1.0+s2, Immirzi*k, MPC_RNDNN);
      complex_lngamma(z, z);
      mpc_sub(w, w, z, MPC_RNDNN);

      mpc_set_d_d(z, (1.0+2.0*k)/2.0, -Immirzi*k, MPC_RNDNN);
      complex_lngamma(z, z);
      mpc_sub(w, w, z, MPC_RNDNN);

      mpc_exp(w, w, MPC_RNDNN);

      mpc_mul(s2prod, s2prod, w, MPC_RNDNN);
      mpc_add( s2sum, s2sum, s2prod, MPC_RNDNN);
  
  }

  mpc_mul(result,result,s2sum,MPC_RNDNN);

  mpc_real(*chi_value, result, MPC_RNDNN ); 

  mpfr_clear(prefactor_1);
  mpfr_clear(prefactor_2);
  mpc_clear(prefactor);
  mpc_clear(s2prod);
  mpc_clear(s2sum);
  mpc_clear(result);
	mpc_clear(z);
	mpc_clear(q);
	mpc_clear(w);
  mpfr_clear(a);
	mpfr_clear(x);
  mpfr_clear(y);


}

double wigner_3j_sl2c(dspin two_j1, dspin two_j2, dspin two_j3,
                      dspin two_k1, dspin two_k2, dspin two_k3,
                      float two_rho1, float two_rho2, float two_rho3) {

   assert(two_j1 >= two_k1);
   assert(two_j2 >= two_k2);
   assert(two_j3 >= two_k3);

  double cutoff = 1e-20;
  int phase = real_negpow(two_j1 -two_j2 +two_j3);

  double value;
  value = phase * sqrt(d(two_j3)) * chi(two_j1, two_j2, two_j3,
                                        two_k1, two_k2, two_k3,
                                        two_rho1, two_rho2, two_rho3);

  if (fabs(value) < cutoff) {
    return 0.0;
  }

  return value;

}

double wigner_6j_sl2c(dspin two_k1, dspin two_k2, dspin two_k3,
                      dspin two_k4, dspin two_k5, dspin two_k6,
                      float two_rho1, float two_rho2, float two_rho3,
                      float two_rho4, float two_rho5, float two_rho6,
                      dspin two_Dj) {

  assert(two_Dj >= two_k1);
  assert(two_Dj >= two_k2);
  assert(two_Dj >= two_k3);
  assert(two_Dj >= two_k4);
  assert(two_Dj >= two_k5);
  assert(two_Dj >= two_k6);

  double value = 0.0 + 0.0*I;
  double cutoff = 1e-25;

  #pragma omp parallel reduction (+:value)
  {
    #pragma omp for collapse(2)
    for (dspin two_j1 = two_k1; two_j1 <= two_Dj; two_j1 +=2) {
    for (dspin two_j2 = two_k2; two_j2 <= two_Dj; two_j2 +=2) {
    for (dspin two_j3 = max(two_k3, abs(two_j1 - two_j2));
                      two_j3 <= min(two_Dj , two_j1 + two_j2); two_j3 +=2) {
    for (dspin two_j4 = two_k4; two_j4 <= two_Dj; two_j4 +=2) {
    for (dspin two_j5 = max (two_k5, abs(two_j4-two_j3));
                      two_j5 <= min (two_Dj , two_j4+two_j3); two_j5 +=2) {
    for (dspin two_j6 = max(abs(two_j4-two_j2), max(two_k6, abs(two_j1 - two_j5)));
                      two_j6 <= min(two_j4+two_j2, min ( two_Dj ,two_j1 + two_j5)); two_j6 +=2) {

      double complex wig6j, wig1, wig2, wig3, wig4;

      wig6j = wig6jj(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6);

      wig1 = wigner_3j_sl2c(two_j1, two_j2, two_j3,
                             two_k1, two_k2, two_k3,
                             two_rho1, two_rho2, two_rho3);
      if (fabs(wig1) < cutoff) continue;

      wig2 = wigner_3j_sl2c(two_j1, two_j5, two_j6,
                             two_k1, two_k5, two_k6,
                             two_rho1, two_rho5, two_rho6);
      if (fabs(wig2) < cutoff) continue;

      wig3 = wigner_3j_sl2c(two_j4, two_j2, two_j6,
                             two_k4, two_k2, two_k6,
                             two_rho4, two_rho2, two_rho6);
      if (fabs(wig3) < cutoff) continue;

      wig4 = wigner_3j_sl2c(two_j4, two_j5, two_j3,
                             two_k4, two_k5, two_k3,
                             two_rho4, two_rho5, two_rho3);
      if (fabs(wig4) < cutoff) continue;

      value += wig6j * wig1 * wig2 * wig3 * wig4;

    }
    }
    }
    }
    }
    }
  }

  if (fabs(value) > 1.0) {
    error("Non-sense value for 6j symbol.");
  }

  if (fabs(value) < cutoff ) {
    return 0.0;
  }

  return value;

}
