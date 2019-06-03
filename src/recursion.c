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
#include "config.h"
#include "utilities.h"
#include "hashing.h"
#include "jsymbols.h"
#include "recursion.h"

// function to be used in the
// recursion formulas
void rprefactor (mpc_t *r_factor, mpfr_t l, mpfr_t k, mpfr_t rho ) {

    // r prefactor: 
    // -2*I / l *sqrt((l*l-k*k)*(l*l+rho*rho));

    mpfr_t aux1,aux2,aux3;

    mpfr_init2(aux1,MPBITS);                    
    mpfr_init2(aux2,MPBITS);                    
    mpfr_init2(aux3,MPBITS); 

    mpfr_set(aux1,l,MPFR_RNDN);
    mpfr_mul(aux1,aux1,l,MPFR_RNDN);

    mpfr_set(aux2,k,MPFR_RNDN);
    mpfr_mul(aux2,aux2,k,MPFR_RNDN);
    

    mpfr_set(aux3,rho,MPFR_RNDN);
    mpfr_mul(aux3,aux3,rho,MPFR_RNDN);

    mpfr_sub(aux2,aux1,aux2,MPFR_RNDN);

    mpfr_add(aux3,aux1,aux3,MPFR_RNDN);

    mpfr_mul(aux2,aux2,aux3,MPFR_RNDN);
    mpfr_sqrt(aux2,aux2,MPFR_RNDN);

    mpc_set_fr(*r_factor, aux2, MPC_RNDNN);
    mpc_mul_i(*r_factor,*r_factor,-1,MPC_RNDNN);
    mpc_mul_ui(*r_factor,*r_factor,2,MPC_RNDNN);

    if ( mpfr_get_d(l,MPFR_RNDN) > 0.1){
        mpc_div_fr(*r_factor,*r_factor,l,MPC_RNDNN);
    }

    mpfr_clears(aux1,aux2,aux3,NULL);


}

// this function gets the previously computed chi value     
// that is stored in chi_array and we put it in chi_value       
void chi_value_get (mpfr_t *chi_value, mpfr_t*** chi_array,
                dspin two_l1, dspin two_l2, dspin two_l3,
                dspin two_k1, dspin two_k2, dspin two_k3) { 

        // check for zeros            
        if (two_l1 < two_k1 || two_l2 < two_k2 || two_l3 < two_k3 || 
            two_l3 < abs(two_l1-two_l2) || two_l3 > two_l1+two_l2 || 
            two_l2 < abs(two_l1-two_l3) || two_l2 > two_l1+two_l3 || 
            two_l1 < abs(two_l3-two_l2) || two_l1 > two_l3+two_l2){
            mpfr_set_ui(*chi_value,0,MPFR_RNDN);

        } else {   
            mpfr_set(*chi_value, chi_array[two_l1-two_k1][two_l2-two_k2][two_l3-two_k3],MPFR_RNDN);
        }

}

//first recursion relation    
void recursion_1 (mpfr_t *chi_value, mpfr_t*** chi_array,
                  mpfr_t l1, mpfr_t l2, mpfr_t l3,
                  mpfr_t k1, mpfr_t k2, mpfr_t k3,
                  mpfr_t rho1, mpfr_t rho2, mpfr_t rho3, 
                  dspin two_l1, dspin two_l2, dspin two_l3,
                  dspin two_k1, dspin two_k2, dspin two_k3,
                  float two_rho1, float two_rho2, float two_rho3,
                  float Immirzi ) {              

    mpfr_t aux1,aux2,aux3,aux4,aux5,aux6;

    mpfr_init2(aux1, MPBITS);                    
    mpfr_init2(aux2, MPBITS);                    
    mpfr_init2(aux3, MPBITS);                    
    mpfr_init2(aux4, MPBITS);
    mpfr_init2(aux5, MPBITS);
    mpfr_init2(aux6, MPBITS);

    mpc_t r_factor;

    mpc_init2(r_factor,MPBITS); 

    mpc_t gf, fi ,se, th, fo;
    mpc_init2 (gf, MPBITS);
    mpc_init2 (fi, MPBITS);
    mpc_init2 (se, MPBITS);
    mpc_init2 (th, MPBITS);
    mpc_init2 (fo, MPBITS);

    // start with the computation of the global factor
    // global_factor_1 = sqrt((2*l1-1)*(2*l1+1))/rprefactor(l1,k1,rho1)*
    //                   1/sqrt((l1+l2+l3+1)*(l1+l2+l3)*(l1+l2-l3)*(l1+l3-l2));

    //Everything has to be done at arbitrary precision
    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);
    mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

    mpfr_set(aux2,l1,MPFR_RNDN);
    mpfr_mul_ui(aux2,aux2,2,MPFR_RNDN);
    mpfr_add_ui(aux2,aux2,1,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);

    mpc_set_fr(gf,aux1,MPC_RNDNN);
    
    // compute the first r_prefactor
    // the value is stored in r_factor
    rprefactor (&r_factor, l1, k1, rho1);  

    mpc_div(gf,gf,r_factor,MPC_RNDNN);

    // start with the denominator
    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_add(aux1,aux1,l2,MPFR_RNDN);
    mpfr_add(aux1,aux1,l3,MPFR_RNDN);
    mpfr_add_ui(aux1,aux1,1,MPFR_RNDN);

    mpfr_set(aux2,l1,MPFR_RNDN);
    mpfr_add(aux2,aux2,l2,MPFR_RNDN);
    mpfr_add(aux2,aux2,l3,MPFR_RNDN);

    mpfr_set(aux3,l1,MPFR_RNDN);
    mpfr_add(aux3,aux3,l2,MPFR_RNDN);
    mpfr_sub(aux3,aux3,l3,MPFR_RNDN);

    mpfr_set(aux4,l1,MPFR_RNDN);
    mpfr_add(aux4,aux4,l3,MPFR_RNDN);
    mpfr_sub(aux4,aux4,l2,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
    mpfr_mul(aux1,aux1,aux4,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);

    mpc_div_fr(gf,gf,aux1,MPC_RNDNN);

    // end of the global factor. 

    // starting the first factor
    // first_1 = 2*sqrt((l1+l2+l3)*(l3+l2-l1+1))* (k3*rho3/l3-k2*rho2/l2-k1*rho1/l1*(l3-l2)/(l1-1))*I;

    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_add(aux1,aux1,l2,MPFR_RNDN);
    mpfr_add(aux1,aux1,l3,MPFR_RNDN);

    mpfr_set(aux2,l3,MPFR_RNDN);
    mpfr_add(aux2,aux2,l2,MPFR_RNDN);
    mpfr_sub(aux2,aux2,l1,MPFR_RNDN);
    mpfr_add_ui(aux2,aux2,1,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);
    mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);

    mpc_set_fr(fi,aux1,MPC_RNDNN);

    if (two_l3 == 0 && two_k3 == 0) {
        mpfr_set(aux1,rho3,MPFR_RNDN);
    } else {
        mpfr_set(aux1,k3,MPFR_RNDN);
        mpfr_mul(aux1,aux1,rho3,MPFR_RNDN);
        mpfr_div(aux1,aux1,l3,MPFR_RNDN);
    }

    if (two_l2 == 0 && two_k2 == 0) {
        mpfr_set(aux2,rho2,MPFR_RNDN);
    } else {
        mpfr_set(aux2,k2,MPFR_RNDN);
        mpfr_mul(aux2,aux2,rho2,MPFR_RNDN);
        mpfr_div(aux2,aux2,l2,MPFR_RNDN);
    }

    if (two_k1 == 0 && two_l1 == 0) {        
            mpfr_set(aux3,rho1,MPFR_RNDN);
            mpfr_mul_si(aux3,aux3,-1,MPFR_RNDN);
            mpfr_set_ui(aux5,1,MPFR_RNDN);
    } else if (two_k1 == 0 && two_l1 == 2) {
            mpfr_set(aux3,rho1,MPFR_RNDN);
            mpfr_set_ui(aux5,1,MPFR_RNDN);
    } else {
        mpfr_set(aux3,k1,MPFR_RNDN);
        mpfr_mul(aux3,aux3,rho1,MPFR_RNDN);
        mpfr_div(aux3,aux3,l1,MPFR_RNDN);
        mpfr_set(aux5,l1,MPFR_RNDN);
        mpfr_sub_ui(aux5,aux5,1,MPFR_RNDN);
    }
    
    mpfr_set(aux4,l3,MPFR_RNDN);
    mpfr_sub(aux4,aux4,l2,MPFR_RNDN);
    mpfr_div(aux4,aux4,aux5,MPFR_RNDN);

    mpfr_mul(aux3,aux3,aux4,MPFR_RNDN);

    mpfr_sub(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sub(aux1,aux1,aux3,MPFR_RNDN);

    mpc_mul_fr(fi,fi,aux1,MPC_RNDNN);
    mpc_mul_i(fi,fi,1,MPC_RNDNN);

    // get the chi value for the first factor
    chi_value_get(chi_value, chi_array, two_l1-2, two_l2, two_l3, two_k1, two_k2, two_k3);
   
    mpc_mul_fr (fi,fi,*chi_value,MPC_RNDNN);


    // end of the first factor. 

    // starting the second factor                     
    // second_1 =  rprefactor(l3,k3,rho3)*sqrt((l3+l1-l2-1)*(l1+l2-l3));

    rprefactor (&r_factor,l3,k3,rho3);

    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_add(aux1,aux1,l3,MPFR_RNDN);
    mpfr_sub(aux1,aux1,l2,MPFR_RNDN);
    mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

    mpfr_set(aux2,l1,MPFR_RNDN);
    mpfr_add(aux2,aux2,l2,MPFR_RNDN);
    mpfr_sub(aux2,aux2,l3,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);

    mpc_set_fr(se,aux1,MPC_RNDNN);
    mpc_mul(se,se,r_factor,MPC_RNDNN);

    chi_value_get(chi_value, chi_array, two_l1-2, two_l2, two_l3-2, two_k1, two_k2, two_k3);

    mpc_mul_fr (se,se,*chi_value, MPC_RNDNN);             
    

    // end of the second factor. 

    // starting the third factor
    // third_1 = rprefactor(l2,k2,rho2)*sqrt((2*l2+1)*(l1+l3-l2)*(l1+l2-l3-1)/(2*l2-1));

    if (two_l2 > 1) {
        
        rprefactor (&r_factor,l2,k2,rho2);

        mpfr_set(aux1,l2,MPFR_RNDN);
        mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);
        mpfr_add_ui(aux1,aux1,1,MPFR_RNDN);

        mpfr_set(aux2,l1,MPFR_RNDN);
        mpfr_add(aux2,aux2,l3,MPFR_RNDN);
        mpfr_sub(aux2,aux2,l2,MPFR_RNDN);

        mpfr_set(aux3,l1,MPFR_RNDN);
        mpfr_add(aux3,aux3,l2,MPFR_RNDN);
        mpfr_sub(aux3,aux3,l3,MPFR_RNDN);
        mpfr_sub_ui(aux3,aux3,1,MPFR_RNDN);

        mpfr_set(aux4,l2,MPFR_RNDN);
        mpfr_mul_ui(aux4,aux4,2,MPFR_RNDN);
        mpfr_sub_ui(aux4,aux4,1,MPFR_RNDN);

        mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
        mpfr_div(aux1,aux1,aux4,MPFR_RNDN);
        
        mpfr_sqrt(aux1,aux1,MPFR_RNDN);

        mpc_set_fr(th,aux1,MPC_RNDNN);
        mpc_mul(th,th,r_factor,MPC_RNDNN);

        chi_value_get(chi_value, chi_array, two_l1-2, two_l2-2, two_l3, two_k1, two_k2, two_k3);

        mpc_mul_fr(th,th,*chi_value, MPC_RNDNN);                         

    } else {
        mpc_set_dc(th, 0, MPC_RNDNN );
    }              

    // end of the third factor. 

    // starting the fourth factor
    // fourth_1 = -rprefactor(l1-1,k1,rho1)*sqrt((l3+l1-l2-1)*(l3+l2-l1+2)*(l3+l2-l1+1)*(l1-1+l2-l3)/((2*l1-1)*(2*l1-3)));       

    if (two_l1 > 3) {

        mpfr_set(aux1,l1,MPFR_RNDN);
        mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);
        rprefactor (&r_factor,aux1,k1,rho1);

        mpfr_set(aux1,l3,MPFR_RNDN);
        mpfr_add(aux1,aux1,l1,MPFR_RNDN);
        mpfr_sub(aux1,aux1,l2,MPFR_RNDN);
        mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

        mpfr_set(aux2,l3,MPFR_RNDN);
        mpfr_add(aux2,aux2,l2,MPFR_RNDN);
        mpfr_sub(aux2,aux2,l1,MPFR_RNDN);
        mpfr_add_ui(aux2,aux2,2,MPFR_RNDN);

        mpfr_set(aux3,l3,MPFR_RNDN);
        mpfr_add(aux3,aux3,l2,MPFR_RNDN);
        mpfr_sub(aux3,aux3,l1,MPFR_RNDN);
        mpfr_add_ui(aux3,aux3,1,MPFR_RNDN);

        mpfr_set(aux4,l1,MPFR_RNDN);
        mpfr_add(aux4,aux4,l2,MPFR_RNDN);
        mpfr_sub(aux4,aux4,l3,MPFR_RNDN);
        mpfr_sub_ui(aux4,aux4,1,MPFR_RNDN);

        mpfr_set(aux5,l1,MPFR_RNDN);
        mpfr_mul_ui(aux5,aux5,2,MPFR_RNDN);
        mpfr_sub_ui(aux5,aux5,1,MPFR_RNDN);

        mpfr_set(aux6,l1,MPFR_RNDN);
        mpfr_mul_ui(aux6,aux6,2,MPFR_RNDN);
        mpfr_sub_ui(aux6,aux6,3,MPFR_RNDN);      

        mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux4,MPFR_RNDN);
        mpfr_div(aux1,aux1,aux5,MPFR_RNDN);
        mpfr_div(aux1,aux1,aux6,MPFR_RNDN);
      
        mpfr_sqrt(aux1,aux1,MPFR_RNDN);
        mpfr_mul_si(aux1,aux1,-1,MPFR_RNDN);
        
        mpc_set_fr(fo,aux1,MPC_RNDNN);
        mpc_mul(fo,fo,r_factor,MPC_RNDNN);

        chi_value_get(chi_value, chi_array, two_l1-4, two_l2, two_l3, two_k1, two_k2, two_k3);

        mpc_mul_fr(fo,fo,*chi_value, MPC_RNDNN);             
                    
    } else {
        mpc_set_dc(fo, 0, MPC_RNDNN );
    }          

    mpc_add (fi, fi, se, MPC_RNDNN);      
    mpc_add (fi, fi, th, MPC_RNDNN);
    mpc_add (fi, fi, fo, MPC_RNDNN);
    mpc_mul (fi, fi, gf, MPC_RNDNN);            

    mpc_real(*chi_value,fi,MPC_RNDNN);             

    mpc_clear (gf);
    mpc_clear (fi);
    mpc_clear (se);
    mpc_clear (th);
    mpc_clear (fo);

    mpc_clear(r_factor);
    mpfr_clears(aux1,aux2,aux3,aux4,aux5,aux6,NULL);

}

// second recursion relation    
void recursion_2 (mpfr_t *chi_value, mpfr_t*** chi_array, 
                  mpfr_t l1, mpfr_t l2, mpfr_t l3,
                  mpfr_t k1, mpfr_t k2, mpfr_t k3,
                  mpfr_t rho1, mpfr_t rho2, mpfr_t rho3,
                  dspin two_l1, dspin two_l2, dspin two_l3,
                  dspin two_k1, dspin two_k2, dspin two_k3,
                  float two_rho1, float two_rho2, float two_rho3,
                  float Immirzi ) {

    mpfr_t aux1,aux2,aux3,aux4,aux5,aux6;

    mpfr_init2(aux1, MPBITS);                    
    mpfr_init2(aux2, MPBITS);                    
    mpfr_init2(aux3, MPBITS);                    
    mpfr_init2(aux4, MPBITS);
    mpfr_init2(aux5, MPBITS);
    mpfr_init2(aux6, MPBITS);

    mpc_t r_factor;

    mpc_init2(r_factor, MPBITS); 

    mpc_t gf, fi ,se, th, fo;
    mpc_init2 (gf, MPBITS);
    mpc_init2 (fi, MPBITS);
    mpc_init2 (se, MPBITS);
    mpc_init2 (th, MPBITS);
    mpc_init2 (fo, MPBITS);

    
    // start with the computation of the global factor
    // global_factor_2 = sqrt((2*l2-1)*(2*l2+1))/rprefactor(l2,k2,rho2)*
    //                   1/sqrt((l1+l2+l3+1)*(l1+l2+l3)*(l1+l2-l3)*(l3+l2-l1));

    mpfr_set(aux1,l2,MPFR_RNDN);
    mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);
    mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

    mpfr_set(aux2,l2,MPFR_RNDN);
    mpfr_mul_ui(aux2,aux2,2,MPFR_RNDN);
    mpfr_add_ui(aux2,aux2,1,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);

    mpc_set_fr(gf,aux1,MPC_RNDNN);

    // compute the first r_prefactor
    // the value is stored in r_factor
    rprefactor (&r_factor, l2, k2, rho2);

    mpc_div(gf,gf,r_factor,MPC_RNDNN);

    // start with the denominator
    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_add(aux1,aux1,l2,MPFR_RNDN);
    mpfr_add(aux1,aux1,l3,MPFR_RNDN);
    mpfr_add_ui(aux1,aux1,1,MPFR_RNDN);

    mpfr_set(aux2,l1,MPFR_RNDN);
    mpfr_add(aux2,aux2,l2,MPFR_RNDN);
    mpfr_add(aux2,aux2,l3,MPFR_RNDN);

    mpfr_set(aux3,l1,MPFR_RNDN);
    mpfr_add(aux3,aux3,l2,MPFR_RNDN);
    mpfr_sub(aux3,aux3,l3,MPFR_RNDN);

    mpfr_set(aux4,l3,MPFR_RNDN);
    mpfr_add(aux4,aux4,l2,MPFR_RNDN);
    mpfr_sub(aux4,aux4,l1,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
    mpfr_mul(aux1,aux1,aux4,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);

    mpc_div_fr(gf,gf,aux1,MPC_RNDNN);

    // end of the global factor

    // starting the first factor
    // formula: 
    // first_2 = -2*sqrt((l1+l2+l3)*(l3+l1-l2+1))* (k3*rho3/l3-k1*rho1/l1-k2*rho2/l2*(l3-l1)/(l2-1))*I; 

    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_add(aux1,aux1,l2,MPFR_RNDN);
    mpfr_add(aux1,aux1,l3,MPFR_RNDN);

    mpfr_set(aux2,l3,MPFR_RNDN);
    mpfr_add(aux2,aux2,l1,MPFR_RNDN);
    mpfr_sub(aux2,aux2,l2,MPFR_RNDN);
    mpfr_add_ui(aux2,aux2,1,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);
    mpfr_mul_si(aux1,aux1,-2,MPFR_RNDN);

    mpc_set_fr(fi,aux1,MPC_RNDNN);

    if (two_l3 == 0 && two_k3 == 0) {
        mpfr_set(aux1,rho3,MPFR_RNDN);
    } else {
        mpfr_set(aux1,k3,MPFR_RNDN);
        mpfr_mul(aux1,aux1,rho3,MPFR_RNDN);
        mpfr_div(aux1,aux1,l3,MPFR_RNDN);
    }

    if (two_l1 == 0 && two_k1 == 0) {
        mpfr_set(aux2,rho1,MPFR_RNDN);
    } else {
        mpfr_set(aux2,k1,MPFR_RNDN);
        mpfr_mul(aux2,aux2,rho1,MPFR_RNDN);
        mpfr_div(aux2,aux2,l1,MPFR_RNDN);
    }

    if (two_k2 == 0 && two_l2 == 0) {        
            mpfr_set(aux3,rho2,MPFR_RNDN);
            mpfr_mul_si(aux3,aux3,-1,MPFR_RNDN);
            mpfr_set_ui(aux5,1,MPFR_RNDN);

    } else if (two_k2 == 0 && two_l2 == 2) {
            mpfr_set(aux3,rho2,MPFR_RNDN);
            mpfr_set_ui(aux5,1,MPFR_RNDN);
    } else {

        mpfr_set(aux3,k2,MPFR_RNDN);
        mpfr_mul(aux3,aux3,rho2,MPFR_RNDN);
        mpfr_div(aux3,aux3,l2,MPFR_RNDN);
        mpfr_set(aux5,l2,MPFR_RNDN);
        mpfr_sub_ui(aux5,aux5,1,MPFR_RNDN);

    }

    mpfr_set(aux4,l3,MPFR_RNDN);
    mpfr_sub(aux4,aux4,l1,MPFR_RNDN);
    mpfr_div(aux4,aux4,aux5,MPFR_RNDN);

    mpfr_mul(aux3,aux3,aux4,MPFR_RNDN);
    mpfr_sub(aux1,aux1,aux2,MPFR_RNDN);

    mpfr_sub(aux1,aux1,aux3,MPFR_RNDN);

    mpc_mul_fr(fi,fi,aux1,MPC_RNDNN);
    mpc_mul_i(fi,fi,1,MPC_RNDNN);

    // get the chi value for the first factor
    chi_value_get(chi_value, chi_array, two_l1, two_l2-2, two_l3, two_k1, two_k2, two_k3);

    mpc_mul_fr (fi,fi,*chi_value,MPC_RNDNN);

    // end of the first factor
        
    // starting the second factor       
    // formula:              
    // second_2 = rprefactor(l3,k3,rho3)*sqrt((l3-l1+l2-1)*(l1+l2-l3));

    rprefactor (&r_factor,l3,k3,rho3);

    mpfr_set(aux1,l3,MPFR_RNDN);
    mpfr_add(aux1,aux1,l2,MPFR_RNDN);
    mpfr_sub(aux1,aux1,l1,MPFR_RNDN);
    mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

    mpfr_set(aux2,l1,MPFR_RNDN);
    mpfr_add(aux2,aux2,l2,MPFR_RNDN);
    mpfr_sub(aux2,aux2,l3,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);

    mpc_set_fr(se,aux1,MPC_RNDNN);
    mpc_mul(se,se,r_factor,MPC_RNDNN);

    chi_value_get(chi_value, chi_array, two_l1, two_l2-2, two_l3-2, two_k1, two_k2, two_k3);

    mpc_mul_fr (se,se,*chi_value, MPC_RNDNN);

    // end of the second factor

    // starting the third factor
    // formula:
    // third_2 = rprefactor(l1,k1,rho1)*sqrt((2*l1+1)*(l3-l1+l2)*(l1+l2-l3-1)/(2*l1-1));
            
    if (two_l1 > 1) {

        rprefactor (&r_factor,l1,k1,rho1);

        mpfr_set(aux1,l1,MPFR_RNDN);
        mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);
        mpfr_add_ui(aux1,aux1,1,MPFR_RNDN);

        mpfr_set(aux2,l3,MPFR_RNDN);
        mpfr_add(aux2,aux2,l2,MPFR_RNDN);
        mpfr_sub(aux2,aux2,l1,MPFR_RNDN);

        mpfr_set(aux3,l1,MPFR_RNDN);
        mpfr_add(aux3,aux3,l2,MPFR_RNDN);
        mpfr_sub(aux3,aux3,l3,MPFR_RNDN);
        mpfr_sub_ui(aux3,aux3,1,MPFR_RNDN);

        mpfr_set(aux4,l1,MPFR_RNDN);
        mpfr_mul_ui(aux4,aux4,2,MPFR_RNDN);
        mpfr_sub_ui(aux4,aux4,1,MPFR_RNDN);

        mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
        mpfr_div(aux1,aux1,aux4,MPFR_RNDN);
        
        mpfr_sqrt(aux1,aux1,MPFR_RNDN);

        mpc_set_fr(th,aux1,MPC_RNDNN);
        mpc_mul(th,th,r_factor,MPC_RNDNN);

        chi_value_get(chi_value, chi_array, two_l1-2, two_l2-2, two_l3,  two_k1, two_k2, two_k3);

        mpc_mul_fr(th,th,*chi_value, MPC_RNDNN);  
        
    } else {
        mpc_set_dc(th, 0, MPC_RNDNN );
    }        

    // end of the third factor

    // starting the fourth factor 
    // formula: 
    // fourth_2 = -rprefactor(l2-1,k2,rho2)*sqrt((l3+l2-l1-1)*(l3+l1-l2+2)*(l3+l1-l2+1)*(l2-1+l1-l3)/((2*l2-1)*(2*l2-3)));
 
    if (two_l2 > 3) {  

        mpfr_set(aux1,l2,MPFR_RNDN);
        mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);
        rprefactor (&r_factor,aux1,k2,rho2);

        mpfr_set(aux1,l3,MPFR_RNDN);
        mpfr_add(aux1,aux1,l2,MPFR_RNDN);
        mpfr_sub(aux1,aux1,l1,MPFR_RNDN);
        mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

        mpfr_set(aux2,l3,MPFR_RNDN);
        mpfr_add(aux2,aux2,l1,MPFR_RNDN);
        mpfr_sub(aux2,aux2,l2,MPFR_RNDN);
        mpfr_add_ui(aux2,aux2,2,MPFR_RNDN);

        mpfr_set(aux3,l3,MPFR_RNDN);
        mpfr_add(aux3,aux3,l1,MPFR_RNDN);
        mpfr_sub(aux3,aux3,l2,MPFR_RNDN);
        mpfr_add_ui(aux3,aux3,1,MPFR_RNDN);

        mpfr_set(aux4,l2,MPFR_RNDN);
        mpfr_add(aux4,aux4,l1,MPFR_RNDN);
        mpfr_sub(aux4,aux4,l3,MPFR_RNDN);
        mpfr_sub_ui(aux4,aux4,1,MPFR_RNDN);

        mpfr_set(aux5,l2,MPFR_RNDN);
        mpfr_mul_ui(aux5,aux5,2,MPFR_RNDN);
        mpfr_sub_ui(aux5,aux5,1,MPFR_RNDN);

        mpfr_set(aux6,l2,MPFR_RNDN);
        mpfr_mul_ui(aux6,aux6,2,MPFR_RNDN);
        mpfr_sub_ui(aux6,aux6,3,MPFR_RNDN);   

        mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);  
        mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);  
        mpfr_mul(aux1,aux1,aux4,MPFR_RNDN);
        mpfr_mul(aux5,aux5,aux6,MPFR_RNDN);
        mpfr_div(aux1,aux1,aux5,MPFR_RNDN);
      
        mpfr_sqrt(aux1,aux1,MPFR_RNDN);
        mpfr_mul_si(aux1,aux1,-1,MPFR_RNDN);
        
        mpc_set_fr(fo,aux1,MPC_RNDNN);
        mpc_mul(fo,fo,r_factor,MPC_RNDNN);
        
        chi_value_get(chi_value, chi_array, two_l1, two_l2-4, two_l3,two_k1, two_k2, two_k3);

        mpc_mul_fr(fo,fo,*chi_value, MPC_RNDNN);    

    } else {
        mpc_set_dc(fo, 0, MPC_RNDNN );
    }

    // end of the fourth factor                                         
    
    mpc_add (fi, fi, se, MPC_RNDNN);      
    mpc_add (fi, fi, th, MPC_RNDNN);
    mpc_add (fi, fi, fo, MPC_RNDNN);
    mpc_mul (fi, fi, gf, MPC_RNDNN);            

    mpc_real(*chi_value,fi,MPC_RNDNN);             

    mpc_clear (gf);
    mpc_clear (fi);
    mpc_clear (se);
    mpc_clear (th);
    mpc_clear (fo);

    mpc_clear(r_factor);
    mpfr_clears(aux1,aux2,aux3,aux4,aux5,aux6,NULL);

}

// third recursion relation    
void recursion_3 (mpfr_t *chi_value, mpfr_t*** chi_array, 
                  mpfr_t l1, mpfr_t l2, mpfr_t l3,
                  mpfr_t k1, mpfr_t k2, mpfr_t k3,
                  mpfr_t rho1, mpfr_t rho2, mpfr_t rho3,dspin two_l1, dspin two_l2, dspin two_l3,
                  dspin two_k1, dspin two_k2, dspin two_k3,
                  float two_rho1, float two_rho2, float two_rho3,
                  float Immirzi ) {                

    mpfr_t aux1,aux2,aux3,aux4,aux5,aux6;

    mpfr_init2(aux1,MPBITS);                    
    mpfr_init2(aux2,MPBITS);                    
    mpfr_init2(aux3,MPBITS);                    
    mpfr_init2(aux4,MPBITS);
    mpfr_init2(aux5,MPBITS);
    mpfr_init2(aux6,MPBITS);

    mpc_t r_factor;

    mpc_init2(r_factor,MPBITS); 

    mpc_t gf, fi ,se, th, fo;
    mpc_init2 (gf, MPBITS);
    mpc_init2 (fi, MPBITS);
    mpc_init2 (se, MPBITS);
    mpc_init2 (th, MPBITS);
    mpc_init2 (fo, MPBITS);
    
    // start with the computation of the global factor  
    // global_factor_3 = (2*l3-1)/(rprefactor(l3,k3,rho3)*sqrt((l1+l2+l3+1)*(l1+l2+l3)*(l3+l1-l2)*(l3-l1+l2)));

    mpfr_set(aux1,l3,MPFR_RNDN);
    mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);
    mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

    mpc_set_fr(gf,aux1,MPC_RNDNN);

    // compute the first r_prefactor
    // the value is stored in r_factor
    rprefactor (&r_factor, l3, k3, rho3);

    mpc_div(gf,gf,r_factor,MPC_RNDNN);

    // start with the denominator
    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_add(aux1,aux1,l2,MPFR_RNDN);
    mpfr_add(aux1,aux1,l3,MPFR_RNDN);
    mpfr_add_ui(aux1,aux1,1,MPFR_RNDN);

    mpfr_set(aux2,l1,MPFR_RNDN);
    mpfr_add(aux2,aux2,l2,MPFR_RNDN);
    mpfr_add(aux2,aux2,l3,MPFR_RNDN);

    mpfr_set(aux3,l3,MPFR_RNDN);
    mpfr_add(aux3,aux3,l1,MPFR_RNDN);
    mpfr_sub(aux3,aux3,l2,MPFR_RNDN);

    mpfr_set(aux4,l3,MPFR_RNDN);
    mpfr_add(aux4,aux4,l2,MPFR_RNDN);
    mpfr_sub(aux4,aux4,l1,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
    mpfr_mul(aux1,aux1,aux4,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);

    mpc_div_fr(gf,gf,aux1,MPC_RNDNN);

    // end of the global factor. 
    // starting the first factor
    // first_3 = 2*sqrt((l1+l2+l3)*(-l3+l1+l2+1))* (k3*rho3/(l3*(l3-1))*(l1-l2)-k1*rho1/l1+k2*rho2/l2)*I;

    mpfr_set(aux1,l1,MPFR_RNDN);
    mpfr_add(aux1,aux1,l2,MPFR_RNDN);
    mpfr_add(aux1,aux1,l3,MPFR_RNDN);

    mpfr_set(aux2,l1,MPFR_RNDN);
    mpfr_add(aux2,aux2,l2,MPFR_RNDN);
    mpfr_sub(aux2,aux2,l3,MPFR_RNDN);
    mpfr_add_ui(aux2,aux2,1,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_sqrt(aux1,aux1,MPFR_RNDN);
    mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);

    mpc_set_fr(fi,aux1,MPC_RNDNN);

    if (two_k3 == 0 && two_l3 == 0) {
        mpfr_set(aux1,rho3,MPFR_RNDN);
        mpfr_mul_si(aux1,aux1,-1,MPFR_RNDN);
    } else if (two_k3 == 0 && two_l3 == 2) {
            mpfr_set(aux1,rho3,MPFR_RNDN);
    } else {
        mpfr_set(aux1,k3,MPFR_RNDN);
        mpfr_mul(aux1,aux1,rho3,MPFR_RNDN);
        mpfr_set(aux4,l3,MPFR_RNDN);
        mpfr_set(aux5,l3,MPFR_RNDN);
        mpfr_sub_ui(aux5,aux5,1,MPFR_RNDN);
        mpfr_mul(aux4,aux4,aux5,MPFR_RNDN);

        mpfr_div(aux1,aux1,aux4,MPFR_RNDN);
    }
   
    mpfr_set(aux5,l1,MPFR_RNDN);
    mpfr_sub(aux5,aux5,l2,MPFR_RNDN);

    mpfr_mul(aux1,aux1,aux5,MPFR_RNDN);

    mpfr_set(aux2,k1,MPFR_RNDN);
    mpfr_mul(aux2,aux2,rho1,MPFR_RNDN);
    mpfr_div(aux2,aux2,l1,MPFR_RNDN);

    mpfr_set(aux3,k2,MPFR_RNDN);
    mpfr_mul(aux3,aux3,rho2,MPFR_RNDN);
    mpfr_div(aux3,aux3,l2,MPFR_RNDN);

    mpfr_sub(aux1,aux1,aux2,MPFR_RNDN);
    mpfr_add(aux1,aux1,aux3,MPFR_RNDN);

    mpc_mul_fr(fi,fi,aux1,MPC_RNDNN);
    mpc_mul_i(fi,fi,1,MPC_RNDNN);

    // get the chi value for the first factor
    chi_value_get(chi_value, chi_array, two_l1, two_l2, two_l3-2, two_k1, two_k2, two_k3);

    mpc_mul_fr (fi,fi,*chi_value,MPC_RNDNN);

    // end of the first factor.            
    // starting the second factor
    // second_3 = -rprefactor(l3-1,k3,rho3)*sqrt((l1-l2+l3-1)*(-l1+l2+l3-1)*(l1+l2-l3+1)*(l1+l2-l3+2))/(2*l3-1);                     

    if (two_l3 > 1) {  

        mpfr_set(aux1,l3,MPFR_RNDN);
        mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);      

        rprefactor (&r_factor,aux1,k3,rho3);

        mpfr_set(aux1,l1,MPFR_RNDN);
        mpfr_add(aux1,aux1,l3,MPFR_RNDN);
        mpfr_sub(aux1,aux1,l2,MPFR_RNDN);
        mpfr_sub_ui(aux1,aux1,1,MPFR_RNDN);

        mpfr_set(aux2,l2,MPFR_RNDN);
        mpfr_add(aux2,aux2,l3,MPFR_RNDN);
        mpfr_sub(aux2,aux2,l1,MPFR_RNDN);
        mpfr_sub_ui(aux2,aux2,1,MPFR_RNDN);

        mpfr_set(aux3,l1,MPFR_RNDN);
        mpfr_add(aux3,aux3,l2,MPFR_RNDN);
        mpfr_sub(aux3,aux3,l3,MPFR_RNDN);
        mpfr_add_ui(aux3,aux3,1,MPFR_RNDN);

        mpfr_set(aux4,l1,MPFR_RNDN);
        mpfr_add(aux4,aux4,l2,MPFR_RNDN);
        mpfr_sub(aux4,aux4,l3,MPFR_RNDN);
        mpfr_add_ui(aux4,aux4,2,MPFR_RNDN);

        mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux4,MPFR_RNDN);
        mpfr_sqrt(aux1,aux1,MPFR_RNDN);

        mpfr_set(aux5,l3,MPFR_RNDN);
        mpfr_mul_ui(aux5,aux5,2,MPFR_RNDN);
        mpfr_sub_ui(aux5,aux5,1,MPFR_RNDN);

        mpfr_div(aux1,aux1,aux5,MPFR_RNDN);

        mpfr_mul_si(aux1,aux1,-1,MPFR_RNDN);

        mpc_set_fr(se,aux1,MPC_RNDNN);
        mpc_mul(se,se,r_factor,MPC_RNDNN);

        chi_value_get(chi_value, chi_array, two_l1, two_l2, two_l3-4, two_k1, two_k2, two_k3);

        mpc_mul_fr (se,se,*chi_value, MPC_RNDNN);             
       
    } else {
        mpc_set_dc(se, 0, MPC_RNDNN );
    }       
        // end of the second factor.

        // starting the third factor
        // third_3 = rprefactor(l1,k1,rho1)*sqrt((2*l1+1)*(l3-l1+l2)*(l1+l3-l2-1)/(2*l1-1));


    if (two_l1 > 1) {

        rprefactor (&r_factor,l1,k1,rho1);

        mpfr_set(aux1,l1,MPFR_RNDN);
        mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);
        mpfr_add_ui(aux1,aux1,1,MPFR_RNDN);

        mpfr_set(aux2,l3,MPFR_RNDN);
        mpfr_add(aux2,aux2,l2,MPFR_RNDN);
        mpfr_sub(aux2,aux2,l1,MPFR_RNDN);

        mpfr_set(aux3,l1,MPFR_RNDN);
        mpfr_add(aux3,aux3,l3,MPFR_RNDN);
        mpfr_sub(aux3,aux3,l2,MPFR_RNDN);
        mpfr_sub_ui(aux3,aux3,1,MPFR_RNDN);

        mpfr_set(aux4,l1,MPFR_RNDN);
        mpfr_mul_ui(aux4,aux4,2,MPFR_RNDN);
        mpfr_sub_ui(aux4,aux4,1,MPFR_RNDN);

        mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
        mpfr_div(aux1,aux1,aux4,MPFR_RNDN);
        
        mpfr_sqrt(aux1,aux1,MPFR_RNDN);

        mpc_set_fr(th,aux1,MPC_RNDNN);
        mpc_mul(th,th,r_factor,MPC_RNDNN);
        
        chi_value_get(chi_value, chi_array, two_l1-2, two_l2, two_l3-2, two_k1, two_k2, two_k3);

        mpc_mul_fr(th,th,*chi_value, MPC_RNDNN);   
                        
    } else {
        mpc_set_dc(th, 0, MPC_RNDNN );
    }   

    // end of the third factor.

    // starting the fourth factor
    // fourth_3 = +rprefactor(l2,k2,rho2)*sqrt((2*l2+1)*(l3-l2+l1)*(l3+l2-l1-1)/(2*l2-1));

    if (two_l2 > 1) { 

        rprefactor (&r_factor,l2,k2,rho2);

        mpfr_set(aux1,l2,MPFR_RNDN);
        mpfr_mul_ui(aux1,aux1,2,MPFR_RNDN);
        mpfr_add_ui(aux1,aux1,1,MPFR_RNDN);

        mpfr_set(aux2,l3,MPFR_RNDN);
        mpfr_add(aux2,aux2,l1,MPFR_RNDN);
        mpfr_sub(aux2,aux2,l2,MPFR_RNDN);

        mpfr_set(aux3,l3,MPFR_RNDN);
        mpfr_add(aux3,aux3,l2,MPFR_RNDN);
        mpfr_sub(aux3,aux3,l1,MPFR_RNDN);
        mpfr_sub_ui(aux3,aux3,1,MPFR_RNDN);

        mpfr_set(aux4,l2,MPFR_RNDN);
        mpfr_mul_ui(aux4,aux4,2,MPFR_RNDN);
        mpfr_sub_ui(aux4,aux4,1,MPFR_RNDN);    

        mpfr_mul(aux1,aux1,aux2,MPFR_RNDN);
        mpfr_mul(aux1,aux1,aux3,MPFR_RNDN);
        mpfr_div(aux1,aux1,aux4,MPFR_RNDN);
      
        mpfr_sqrt(aux1,aux1,MPFR_RNDN);
        
        mpc_set_fr(fo,aux1,MPC_RNDNN);
        mpc_mul(fo,fo,r_factor,MPC_RNDNN);

        chi_value_get(chi_value, chi_array, two_l1, two_l2-2, two_l3-2, two_k1, two_k2, two_k3);

        mpc_mul_fr(fo,fo,*chi_value, MPC_RNDNN);  
            
    } else {
        mpc_set_dc(fo, 0, MPC_RNDNN );
    }               
                                
    mpc_add (fi, fi, se, MPC_RNDNN);      
    mpc_add (fi, fi, th, MPC_RNDNN);
    mpc_add (fi, fi, fo, MPC_RNDNN);
    mpc_mul (fi, fi, gf, MPC_RNDNN);            

    mpc_real(*chi_value,fi,MPC_RNDNN);             

    mpc_clear (gf);
    mpc_clear (fi);
    mpc_clear (se);
    mpc_clear (th);
    mpc_clear (fo);

    mpc_clear(r_factor);

    mpfr_clears(aux1,aux2,aux3,aux4,aux5,aux6,NULL); 

}

void recursion (mpfr_t* chi_value, mpfr_t*** chi_array,  
                dspin two_l1, dspin two_l2, dspin two_l3,
                dspin two_k1, dspin two_k2, dspin two_k3,
                float two_rho1, float two_rho2, float two_rho3,
                float Immirzi ) {      

    // check for zeros            
    if (two_l1 < two_k1 || two_l2 < two_k2 || two_l3 < two_k3 || 
        two_l3 < abs(two_l1-two_l2) || two_l3 > two_l1+two_l2 || 
        two_l2 < abs(two_l1-two_l3) || two_l2 > two_l1+two_l3 || 
        two_l1 < abs(two_l3-two_l2) || two_l1 > two_l3+two_l2){
        
        mpfr_set_ui(*chi_value,0,MPFR_RNDN);
        return;

    }

    // set all mpfr variables and initialize them                
    mpfr_t l1,l2,l3,k1,k2,k3,rho1,rho2,rho3;

    mpfr_init2(l1,MPBITS);                    
    mpfr_init2(l2,MPBITS);                    
    mpfr_init2(l3,MPBITS);                    
    mpfr_init2(k1,MPBITS);                    
    mpfr_init2(k2,MPBITS);                    
    mpfr_init2(k3,MPBITS);                    
    mpfr_init2(rho1,MPBITS);                    
    mpfr_init2(rho2,MPBITS);                    
    mpfr_init2(rho3,MPBITS);  

    mpfr_set_d(l1,(float)two_l1/2,MPFR_RNDN);
    mpfr_set_d(l2,(float)two_l2/2,MPFR_RNDN);   
    mpfr_set_d(l3,(float)two_l3/2,MPFR_RNDN);   
    mpfr_set_d(k1,(float)two_k1/2,MPFR_RNDN);   
    mpfr_set_d(k2,(float)two_k2/2,MPFR_RNDN);   
    mpfr_set_d(k3,(float)two_k3/2,MPFR_RNDN);   
    mpfr_set_d(rho1,(float)two_rho1/2,MPFR_RNDN);   
    mpfr_set_d(rho2,(float)two_rho2/2,MPFR_RNDN);   
    mpfr_set_d(rho3,(float)two_rho3/2,MPFR_RNDN);     

    int aux = 0;    

    // depending on the spins that you have
    // you use one of the three formulae
    // and you store the value in chi_value
    if (-two_l3+two_l2+two_l1 != 0 && two_l3+two_l1-two_l2 !=0 && two_l1 != two_k1 ){

        aux = 1;
        recursion_1 (chi_value, chi_array,
                    l1,l2,l3,
                    k1,k2,k3,
                    rho1,rho2,rho3,
                    two_l1, two_l2, two_l3,
                    two_k1, two_k2, two_k3,
                    two_rho1, two_rho2, two_rho3,
                    Immirzi );
  
        goto clear;  

    }

    if (two_l3+two_l2-two_l1 != 0 && -two_l3+two_l1+two_l2 != 0 && two_l2 != two_k2){

        aux = 1;
        recursion_2 (chi_value, chi_array,
                    l1,l2,l3,
                    k1,k2,k3,
                    rho1,rho2,rho3,
                    two_l1, two_l2, two_l3,
                    two_k1, two_k2, two_k3,
                    two_rho1, two_rho2, two_rho3,
                    Immirzi ); 
        goto clear;              
            
    }                                      

    if (two_l3+two_l2-two_l1 != 0 && two_l3+two_l1-two_l2 !=0 && two_l3 != two_k3){

        aux = 1;
        recursion_3 (chi_value, chi_array,
                    l1,l2,l3,
                    k1,k2,k3,
                    rho1,rho2,rho3,
                    two_l1, two_l2, two_l3,
                    two_k1, two_k2, two_k3,
                    two_rho1, two_rho2, two_rho3,
                    Immirzi );   
        goto clear;      
            
    }   

    // there are cases for which the previous recursion
    // relations are not valid, they give 0 = 0 
    // these cases are: 0 0 0, 0 L L and permutations
    // and 0 J J 0 L L and permutations.
    // we implement an indipendet simplified chi formula
    // for these cases
    if (aux == 0){
        
        if (two_k1 == 0 && two_k2 == 0 && two_k3 == 0) {
            if (two_l1 == 0) {
                chi_mpc_simp0(chi_value, two_l2, two_k2, Immirzi);
            } else if (two_l2 == 0) {
                chi_mpc_simp0(chi_value, two_l3, two_k3, Immirzi);
            } else if (two_l3 == 0) {

                chi_mpc_simp0(chi_value, two_l1, two_k1, Immirzi);

                // to ensure simmetries we multiply for dimensional 
                // factor. the B3 is symmetric while the chi not
                // completely

                mpfr_mul_d(*chi_value,*chi_value, sqrt(d(two_l1)), MPFR_RNDN);
            } else {
                error("recursion problem");
            }
            goto clear;
        }
        
        // case 0 J J 0 L L
        if (two_k1 == 0) {
            chi_mpc_simp0(chi_value, two_l2, two_k2, Immirzi);
        } else if (two_k2 == 0) {
            chi_mpc_simp0(chi_value, two_l3, two_k3, Immirzi);
        } else if (two_k3 == 0) {
            chi_mpc_simp0(chi_value, two_l1, two_k1, Immirzi);
            mpfr_mul_d(*chi_value,*chi_value, sqrt(d(two_l1)), MPFR_RNDN);
        } else {
            error("recursion problem");
        }
        goto clear;
    }

clear:
    mpfr_clears(l1,l2,l3,k1,k2,k3,rho1,rho2,rho3,NULL);

}
