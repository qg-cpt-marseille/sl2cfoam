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

////////////////////////////////////////////////////////////////
// CREDITS: 
// Functions for Booster Functions via Collet's formula.
//
// These functions have been written by François Collet and
// they have been studied in a related paper (to appear).
// Collet's formula is written using Ruhl conventions.
//
// The following is the original code plus modifications
// for interfacing with the library.
// More details at https://arxiv.org/abs/1504.08329
////////////////////////////////////////////////////////////////

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>

#include "b4function.h"
#include "config.h"
#include "hashing.h"
#include "utilities.h"
#include "wigxjpf.h"
#include "error.h"
#include "sl2cfoam.h"

// Various internal utilities.
static float maxF(float e,float f);
static float minF(float e,float f);
static long double fact(float j);
static void f(mpfr_t f_res,int l,float J,float j,int s,float p);
static void alpha(mpc_t alpha_res,int l,float J,float rho,
                                    float j2,float j1,int s2,int s1,float p,int m);
static void alpha_0(mpfr_t alpha_res,int l,float J,float j2,
                                        float j1,int s2,int s1,float p,int m);
static void sigma(mpfr_t sigma_res,int l,float J,float j2,
                                    float j1,int s2,int s1,float p,float ms);
static void sqrtj(mpfr_t sqrtj_res,int l,float J,float j,float p);
static void Y(mpc_t Y_res,int l,float J,float rho,
                            float j2,float j1,float p,int m);
static void Z(mpfr_t Z_res,int l,float J,
                            float j2,float j1,float p,int m);
static void W(mpfr_t W_res,int l,float J,
                            float j2,float j1,float p,float ms);
static long double complex D(mpc_t Ym[],mpc_t Yn[],int l,float J,float rho,
                                                         float j2,float j1,float p,long double X);
static long double complex D_0(mpfr_t Zm[],mpfr_t Zn[],mpfr_t W[],int l,float J,
                                                            float j2,float j1,float p,long double X);
static long double mesure_x(long double X);

// Integration.
static void I4_Interval(long double complex *I4_inf,long double complex *I4_sup,int N,
                                 long double complex *D1,long double complex *D2,long double complex *D3,long double complex *D4);

// Intertwiners.
static long double intertwiner_j(float Ji,long double trois_j_part1,long double trois_j_part2,
                                                    float m1,float m2,float m3,float m4);

// Complex Gamma (machine precision) to rephase.
static long double complex ComplexGamma(double complex z);
static long double complex dPhase (float j, float l, float rho);

static long double complex B4_Interval(long double complex *B4_inf,long double complex *B4_sup,
                                                                 int Ji2_int,float Ji2_min,int Ji1_int,float Ji1_min,
                                                                 float j1,float j2,float j3,float j4,
                                                                 long double ***trois_j2_part1,long double ***trois_j2_part2,
                                                                 long double ***trois_j1_part1,long double ***trois_j1_part2,
                                                                 long double complex ***I4_inf,long double complex ***I4_sup);

long double  **B4Function (int two_k1, int two_k2, int two_k3, int two_k4,
                                                     float two_rho1, float two_rho2, float two_rho3, float two_rho4,
                                                     int two_j1, int two_j2, int two_j3, int two_j4,
                                                     int two_l1, int two_l2, int two_l3, int two_l4) {

        float j1_1=(float)(two_j1)/2,j1_2=(float)(two_j2)/2,j1_3=(float)(two_j3)/2,j1_4=(float)(two_j4)/2;
        float j2_1=(float)(two_l1)/2,j2_2=(float)(two_l2)/2,j2_3=(float)(two_l3)/2,j2_4=(float)(two_l4)/2;
        float J1=(float)(two_k1)/2,J2=(float)(two_k2)/2,J3=(float)(two_k3)/2,J4=(float)(two_k4)/2;
        float rho1=two_rho1/2 ,rho2=two_rho2/2, rho3=two_rho3/2,rho4=two_rho4/2;

        float j1= minF(j2_1,j1_1);
        float j2= minF(j2_2,j1_2);
        float j3= minF(j2_3,j1_3);
        float j4= minF(j2_4,j1_4);

        int dj1= 2.0*j1+1;
        int dj2= 2.0*j2+1;
        int dj3= 2.0*j3+1;
        int dj4= 2.0*j4+1;

        int Sj1= j1_1+j2_1+1.0;
        int Sj2= j1_2+j2_2+1.0;
        int Sj3= j1_3+j2_3+1.0;
        int Sj4= j1_4+j2_4+1.0;

        float j1_max= maxF(j2_1,j1_1);
        float j2_max= maxF(j2_2,j1_2);
        float j3_max= maxF(j2_3,j1_3);
        float j4_max= maxF(j2_4,j1_4);


        int N= maxF(B4_INTEGRATION_POINTS,(rho1+rho2+rho3+rho4)*10); //points for the integral

        int l1= 128+Sj1*log2(N);   //precision of l bits for D1
        int l2= 128+Sj2*log2(N);   //precision of l bits for D2
        int l3= 128+Sj3*log2(N);   //precision of l bits for D3
        int l4= 128+Sj4*log2(N);   //precision of l bits for D4

        long double complex **D1= malloc(sizeof(long double complex*[dj1]));
        for(int i=0;i<dj1;i++){
                D1[i]= malloc(sizeof(long double complex[N]));
        }
        long double complex **D2= malloc(sizeof(long double complex*[dj2]));
        for(int i=0;i<dj2;i++){
                D2[i]= malloc(sizeof( long double complex[N]));
        }
        long double complex **D3= malloc(sizeof( long double complex*[dj3]));
        for(int i=0;i<dj3;i++){
                D3[i]= malloc(sizeof(long double complex[N]));
        }
        long double complex **D4= malloc(sizeof( long double complex*[dj4]));
        for(int i=0;i<dj4;i++){
                D4[i]= malloc(sizeof( long double complex[N]));
        }

        if((rho1==0)&&(((int)abs(J1))==abs(J1))){
                mpfr_t Zm1[dj1][Sj1],Zn1[dj1][Sj1],W1[dj1][Sj1];
#pragma omp parallel for
                for(int pa=0;pa<dj1;pa++){
                        float p1=pa-j1;
                        int m1_max= j1_1+j2_1-abs(J1-p1);
                        for(int m1=0;m1<=m1_max;m1++){
                                mpfr_init2(Zm1[pa][m1],l1);
                                Z(Zm1[pa][m1],l1,J1,j2_1,j1_1,p1,m1);
                        }
                        int n1_max= j1_1+j2_1-abs(J1+p1);
                        for(int n1=0;n1<=n1_max;n1++){
                                mpfr_init2(Zn1[pa][n1],l1);
                                Z(Zn1[pa][n1],l1,J1,j1_1,j2_1,-p1,n1);
                        }
                        float ms1_min= 0.5*maxF(abs(J1-p1),abs(J1+p1));
                        float ms1_max= j1_1+j2_1-ms1_min;
                        int ms1_bound= ms1_max-ms1_min;
                        for(int ms1=0;ms1<=ms1_bound;ms1++){
                                mpfr_init2(W1[pa][ms1],l1);
                                W(W1[pa][ms1],l1,J1,j2_1,j1_1,-p1,ms1+ms1_min);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pa=0;pa<dj1;pa++){
                                float p1=pa-j1;
                                D1[pa][n]=D_0(Zm1[pa],Zn1[pa],W1[pa],l1,J1,j2_1,j1_1,p1,X);
                        }
                }
        }else{
                mpc_t Ym1[dj1][Sj1],Yn1[dj1][Sj1];
#pragma omp parallel for
                for(int pa=0;pa<dj1;pa++){
                        float p1=pa-j1;
                        int m1_max= j1_1+j2_1-abs(J1-p1);
                        for(int m1=0;m1<=m1_max;m1++){
                                mpc_init2(Ym1[pa][m1],l1);
                                Y(Ym1[pa][m1],l1,J1,rho1,j2_1,j1_1,p1,m1);
                        }
                        int n1_max= j1_1+j2_1-abs(J1+p1);
                        for(int n1=0;n1<=n1_max;n1++){
                                mpc_init2(Yn1[pa][n1],l1);
                                Y(Yn1[pa][n1],l1,J1,-rho1,j1_1,j2_1,-p1,n1);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pa=0;pa<dj1;pa++){
                                float p1=pa-j1;
                                D1[pa][n]=D(Ym1[pa],Yn1[pa],l1,J1,rho1,j2_1,j1_1,p1,X);
                        }
                }
        }

        if((rho2==0)&&(((int)abs(J2))==abs(J2))){
                mpfr_t Zm2[dj2][Sj2],Zn2[dj2][Sj2],W2[dj2][Sj2];
#pragma omp parallel for
                for(int pb=0;pb<dj2;pb++){
                        float p2=pb-j2;
                        int m2_max= j1_2+j2_2-abs(J2-p2);
                        for(int m2=0;m2<=m2_max;m2++){
                                mpfr_init2(Zm2[pb][m2],l2);
                                Z(Zm2[pb][m2],l2,J2,j2_2,j1_2,p2,m2);
                        }
                        int n2_max= j1_2+j2_2-abs(J2+p2);
                        for(int n2=0;n2<=n2_max;n2++){
                                mpfr_init2(Zn2[pb][n2],l2);
                                Z(Zn2[pb][n2],l2,J2,j1_2,j2_2,-p2,n2);
                        }
                        float ms2_min= 0.5*maxF(abs(J2-p2),abs(J2+p2));
                        float ms2_max= j1_2+j2_2-ms2_min;
                        int ms2_bound= ms2_max-ms2_min;
                        for(int ms2=0;ms2<=ms2_bound;ms2++){
                                mpfr_init2(W2[pb][ms2],l2);
                                W(W2[pb][ms2],l2,J2,j2_2,j1_2,-p2,ms2+ms2_min);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pb=0;pb<dj2;pb++){
                                float p2=pb-j2;
                                D2[pb][n]=D_0(Zm2[pb],Zn2[pb],W2[pb],l2,J2,j2_2,j1_2,p2,X);
                        }
                }
        }else{
                mpc_t Ym2[dj2][Sj2],Yn2[dj2][Sj2];
#pragma omp parallel for
                for(int pb=0;pb<dj2;pb++){
                        float p2=pb-j2;
                        int m2_max= j1_2+j2_2-abs(J2-p2);
                        for(int m2=0;m2<=m2_max;m2++){
                                mpc_init2(Ym2[pb][m2],l2);
                                Y(Ym2[pb][m2],l2,J2,rho2,j2_2,j1_2,p2,m2);
                        }
                        int n2_max= j1_2+j2_2-abs(J2+p2);
                        for(int n2=0;n2<=n2_max;n2++){
                                mpc_init2(Yn2[pb][n2],l2);
                                Y(Yn2[pb][n2],l2,J2,-rho2,j1_2,j2_2,-p2,n2);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pb=0;pb<dj2;pb++){
                                float p2=pb-j2;
                                D2[pb][n]=D(Ym2[pb],Yn2[pb],l2,J2,rho2,j2_2,j1_2,p2,X);
                        }
                }
        }

        if((rho3==0)&&(((int)abs(J3))==abs(J3))){
                mpfr_t Zm3[dj3][Sj3],Zn3[dj3][Sj3],W3[dj3][Sj3];
#pragma omp parallel for
                for(int pc=0;pc<dj3;pc++){
                        float p3=pc-j3;
                        int m3_max= j1_3+j2_3-abs(J3-p3);
                        for(int m3=0;m3<=m3_max;m3++){
                                mpfr_init2(Zm3[pc][m3],l3);
                                Z(Zm3[pc][m3],l3,J3,j2_3,j1_3,p3,m3);
                        }
                        int n3_max= j1_3+j2_3-abs(J3+p3);
                        for(int n3=0;n3<=n3_max;n3++){
                                mpfr_init2(Zn3[pc][n3],l3);
                                Z(Zn3[pc][n3],l3,J3,j1_3,j2_3,-p3,n3);
                        }
                        float ms3_min= 0.5*maxF(abs(J3-p3),abs(J3+p3));
                        float ms3_max= j1_3+j2_3-ms3_min;
                        int ms3_bound= ms3_max-ms3_min;
                        for(int ms3=0;ms3<=ms3_bound;ms3++){
                                mpfr_init2(W3[pc][ms3],l3);
                                W(W3[pc][ms3],l3,J3,j2_3,j1_3,-p3,ms3+ms3_min);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pc=0;pc<dj3;pc++){
                                float p3=pc-j3;
                                D3[pc][n]=D_0(Zm3[pc],Zn3[pc],W3[pc],l3,J3,j2_3,j1_3,p3,X);
                        }
                }
        }else{
                mpc_t Ym3[dj3][Sj3],Yn3[dj3][Sj3];
#pragma omp parallel for
                for(int pc=0;pc<dj3;pc++){
                        float p3=pc-j3;
                        int m3_max= j1_3+j2_3-abs(J3-p3);
                        for(int m3=0;m3<=m3_max;m3++){
                                mpc_init2(Ym3[pc][m3],l3);
                                Y(Ym3[pc][m3],l3,J3,rho3,j2_3,j1_3,p3,m3);
                        }
                        int n3_max= j1_3+j2_3-abs(J3+p3);
                        for(int n3=0;n3<=n3_max;n3++){
                                mpc_init2(Yn3[pc][n3],l3);
                                Y(Yn3[pc][n3],l3,J3,-rho3,j1_3,j2_3,-p3,n3);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pc=0;pc<dj3;pc++){
                                float p3=pc-j3;
                                D3[pc][n]=D(Ym3[pc],Yn3[pc],l3,J3,rho3,j2_3,j1_3,p3,X);
                        }
                }
        }

        if((rho4==0)&&(((int)abs(J4))==abs(J4))){
                mpfr_t Zm4[dj4][Sj4],Zn4[dj4][Sj4],W4[dj4][Sj4];
#pragma omp parallel for
                for(int pd=0;pd<dj4;pd++){
                        float p4=pd-j4;
                        int m4_max= j1_4+j2_4-abs(J4-p4);
                        for(int m4=0;m4<=m4_max;m4++){
                                mpfr_init2(Zm4[pd][m4],l4);
                                Z(Zm4[pd][m4],l4,J4,j2_4,j1_4,p4,m4);
                        }
                        int n4_max= j1_4+j2_4-abs(J4+p4);
                        for(int n4=0;n4<=n4_max;n4++){
                                mpfr_init2(Zn4[pd][n4],l4);
                                Z(Zn4[pd][n4],l4,J4,j1_4,j2_4,-p4,n4);
                        }
                        float ms4_min= 0.5*maxF(abs(J4-p4),abs(J4+p4));
                        float ms4_max= j1_4+j2_4-ms4_min;
                        int ms4_bound= ms4_max-ms4_min;
                        for(int ms4=0;ms4<=ms4_bound;ms4++){
                                mpfr_init2(W4[pd][ms4],l4);
                                W(W4[pd][ms4],l4,J4,j2_4,j1_4,-p4,ms4+ms4_min);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pd=0;pd<dj4;pd++){
                                float p4=pd-j4;
                                D4[pd][n]=D_0(Zm4[pd],Zn4[pd],W4[pd],l4,J4,j2_4,j1_4,p4,X);
                        }
                }
        }else{
                mpc_t Ym4[dj4][Sj4],Yn4[dj4][Sj4];
#pragma omp parallel for
                for(int pd=0;pd<dj4;pd++){
                        float p4=pd-j4;
                        int m4_max= j1_4+j2_4-abs(J4-p4);
                        for(int m4=0;m4<=m4_max;m4++){
                                mpc_init2(Ym4[pd][m4],l4);
                                Y(Ym4[pd][m4],l4,J4,rho4,j2_4,j1_4,p4,m4);
                        }
                        int n4_max= j1_4+j2_4-abs(J4+p4);
                        for(int n4=0;n4<=n4_max;n4++){
                                mpc_init2(Yn4[pd][n4],l4);
                                Y(Yn4[pd][n4],l4,J4,-rho4,j1_4,j2_4,-p4,n4);
                        }
                }
#pragma omp parallel for
                for(int n=1;n<N;n++){
                        float X= n/(float)(N);
                        for(int pd=0;pd<dj4;pd++){
                                float p4=pd-j4;
                                D4[pd][n]=D(Ym4[pd],Yn4[pd],l4,J4,rho4,j2_4,j1_4,p4,X);
                        }
                }
        }


        long double complex ***I4_inf= (long double complex ***) malloc(dj1*sizeof( long double complex**));
        for(int i=0;i<dj1;i++){
                I4_inf[i]= (long double complex **) malloc(dj2*sizeof(long double complex*));
                for(int j=0;j<dj2;j++){
                        I4_inf[i][j]= (long double complex *) malloc(dj3*sizeof( long double complex));
                }
        }

        /*for(int k=0;k<dj3;k++){
                for(int i=0;i<dj2;i++){
                        for(int j=0;j<dj1;j++){
                                I4_inf[i][j][k]= 0.+ 0.*I;
                        }
                }
        }*/

        long double complex ***I4_sup= (long double complex ***) malloc(dj1*sizeof( long double complex**));
        for(int i=0;i<dj1;i++){
                I4_sup[i]= (long double complex **) malloc(dj2*sizeof(long double complex*));
                for(int j=0;j<dj2;j++){
                        I4_sup[i][j]= (long double complex *) malloc(dj3*sizeof( long double complex));
                }
        }

        /*for(int k=0;k<dj3;k++){
                for(int i=0;i<dj2;i++){
                        for(int j=0;j<dj1;j++){
                                I4_sup[i][j][k]= 0 +0.*I;
                        }
                }
        }*/

#pragma omp parallel for
        for(int pa=0;pa<dj1;pa++){
                for(int pb=0;pb<dj2;pb++){
                        for(int pc=0;pc<dj3;pc++){
                                int pd= (j1+j2+j3+j4-pa-pb-pc);
                                if((0<=pd)&&(pd<dj4)){
                                        I4_Interval(&I4_inf[pa][pb][pc],&I4_sup[pa][pb][pc],N,D1[pa],D2[pb],D3[pc],D4[pd]);
                                }
                        }
                }
        }

        /*for (int i = 0; i < dj1; i++)
                free(D1[i]);
        free(D1);
        for (int i = 0; i < dj2; i++)
                free(D2[i]);
        free(D2);
        for (int i = 0; i < dj3; i++)
                free(D3[i]);
        free(D3);
        for (int i = 0; i < dj4; i++)
                free(D4[i]);
        free(D4);*/


        float Ji1_min=maxF(fabs(j1_1-j1_2),fabs(j1_3-j1_4));
        float Ji1_max=minF(j1_1+j1_2,j1_3+j1_4);
        int Ji1_bound=Ji1_max-Ji1_min;

        long double ***trois_j1_part1= malloc(sizeof( long double**[Ji1_bound+1]));
        for(int i=0;i<=Ji1_bound;i++){
                trois_j1_part1[i]= malloc(sizeof( long double*[dj1]));
                for(int j=0;j<dj1;j++){
                        trois_j1_part1[i][j]= malloc(sizeof( long double[dj2]));
                }
        }

        for(int Ji1=0;Ji1<=Ji1_bound;Ji1++){
                for(int pa=0;pa<dj1;pa++){
                        float p1=pa-j1;
                        for(int pb=0;pb<dj2;pb++){
                                float p2=pb-j2;
                                trois_j1_part1[Ji1][pa][pb]=wig3jj((int)(2*j1_1),(int)(2*j1_2),(int)(2*(Ji1_min+(float)Ji1)),(int)(2*p1),(int)(2*p2),-(int)(2*(p1+p2)));
                                //printf("%i %i %i %i %i %i\n",(int)(2*j1_1),(int)(2*j1_2),(int)(2*(Ji1_min+(float)Ji1)),(int)(2*p1),(int)(2*p2),-(int)(2*(p1+p2)));

                                //printf("%f\n",trois_j1_part1[Ji1][pa][pb]);
                        }
                }
        }
        long double ***trois_j1_part2= malloc(sizeof( long double**[Ji1_bound+1]));
        for(int i=0;i<=Ji1_bound;i++){
                trois_j1_part2[i]= malloc(sizeof( long double*[dj3]));
                for(int j=0;j<dj3;j++){
                        trois_j1_part2[i][j]= malloc(sizeof( long double[dj4]));
                }
        }

        for(int Ji1=0;Ji1<=Ji1_bound;Ji1++){
                for(int pc=0;pc<dj3;pc++){
                        float p3=pc-j3;
                        for(int pd=0;pd<dj4;pd++){
                                float p4=pd-j4;
                                trois_j1_part2[Ji1][pc][pd]=wig3jj((int)(2*j1_3),(int)(2*j1_4),(int)(2*(Ji1_min+(float)Ji1)),(int)(2*p3),(int)(2*p4),-(int)(2*(p3+p4)));
                        }
                }
        }
        float Ji2_min=maxF(fabs(j2_1-j2_2),fabs(j2_3-j2_4));
        float Ji2_max=minF(j2_1+j2_2,j2_3+j2_4);
        int Ji2_bound=Ji2_max-Ji2_min;
        long double ***trois_j2_part1= malloc(sizeof( long double**[Ji2_bound+1]));
        for(int i=0;i<=Ji2_bound;i++){
                trois_j2_part1[i]= malloc(sizeof( long double*[dj1]));
                for(int j=0;j<dj1;j++){
                        trois_j2_part1[i][j]= malloc(sizeof(long double[dj2]));
                }
        }

        for(int Ji2=0;Ji2<=Ji2_bound;Ji2++){
                for(int pa=0;pa<dj1;pa++){
                        float p1=pa-j1;
                        for(int pb=0;pb<dj2;pb++){
                                float p2=pb-j2;
                                trois_j2_part1[Ji2][pa][pb]=wig3jj((int)(2*j2_1),(int)(2*j2_2),(int)(2*(Ji2_min+(float)Ji2)),(int)(2*p1),(int)(2*p2),-(int)(2*(p1+p2)));

                        }
                }
        }
        long double ***trois_j2_part2= malloc(sizeof( long double**[Ji2_bound+1]));
        for(int i=0;i<=Ji2_bound;i++){
                trois_j2_part2[i]= malloc(sizeof( long double*[dj3]));
                for(int j=0;j<dj3;j++){
                        trois_j2_part2[i][j]= malloc(sizeof( long double[dj4]));
                }
        }
        
        for(int Ji2=0;Ji2<=Ji2_bound;Ji2++){
                for(int pc=0;pc<dj3;pc++){
                        float p3=pc-j3;
                        for(int pd=0;pd<dj4;pd++){
                                float p4=pd-j4;
                                trois_j2_part2[Ji2][pc][pd]=wig3jj((int)(2*j2_3),(int)(2*j2_4),(int)(2*(Ji2_min+(float)Ji2)),(int)(2*p3),(int)(2*p4),-(int)(2*(p3+p4)));

                        }
                }
        }

        long double complex Phase;

        Phase = dPhase (j1_1, j2_1, rho1)*dPhase (j1_2, j2_2, rho2)*dPhase (j1_3, j2_3, rho3)*dPhase (j1_4, j2_4, rho4);

        long double **B4_moy= (long double  **) malloc((Ji2_bound+1)*sizeof( long double  *));
        for(int i=0;i<=Ji2_bound;i++){
                B4_moy[i]= (long double  *) malloc((Ji1_bound+1)*sizeof( long double  ));
        }

        long double complex B4_inf[Ji2_bound+1][Ji1_bound+1];
        long double complex B4_sup[Ji2_bound+1][Ji1_bound+1];

        if( two_j1 == two_j2 && two_l2 == two_l1 ||
             two_j3 == two_j4  && two_l4 == two_l3)
        {

#pragma omp parallel for
                for(int Ji1_int=0;Ji1_int<=Ji1_bound;Ji1_int++){
                        for(int Ji2_int=0;Ji2_int<=Ji2_bound;Ji2_int++){

                                if ( ((Ji1_int + (int)Ji1_min) - (Ji2_int + (int)Ji2_min)) % 2 != 0)
                                {
                                        B4_moy[Ji2_int][Ji1_int] = 0.0;
                                }
                                else
                                {

                                        B4_Interval(&B4_inf[Ji2_int][Ji1_int],&B4_sup[Ji2_int][Ji1_int],
                                                                Ji2_int,Ji2_min,Ji1_int,Ji1_min,
                                                                j1,j2,j3,j4,
                                                                trois_j2_part1,trois_j2_part2,
                                                                trois_j1_part1,trois_j1_part2,
                                                                I4_inf,I4_sup);
                                        B4_moy[Ji2_int][Ji1_int]= creall((long double)0.5*(B4_inf[Ji2_int][Ji1_int]+B4_sup[Ji2_int][Ji1_int])/Phase);
                                }
                        }
                }
        }
        else
        {
#pragma omp parallel for
        for(int Ji1_int=0;Ji1_int<=Ji1_bound;Ji1_int++){
                for(int Ji2_int=0;Ji2_int<=Ji2_bound;Ji2_int++){
                        B4_Interval(&B4_inf[Ji2_int][Ji1_int],&B4_sup[Ji2_int][Ji1_int],
                                                Ji2_int,Ji2_min,Ji1_int,Ji1_min,
                                                j1,j2,j3,j4,
                                                trois_j2_part1,trois_j2_part2,
                                                trois_j1_part1,trois_j1_part2,
                                                I4_inf,I4_sup);

                        B4_moy[Ji2_int][Ji1_int]= creall((long double)0.5*(B4_inf[Ji2_int][Ji1_int]+B4_sup[Ji2_int][Ji1_int])/Phase);
                }
        }
        }

        free(D1);
        free(D2);
        free(D3);
        free(D4);
        free(I4_inf);
        free(I4_sup);
        free(trois_j1_part2);
        free(trois_j1_part1);
        free(trois_j2_part2);
        free(trois_j2_part1);

        return B4_moy;
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

float maxF(float e,float f){
        if(e>=f){
                return e;
        }else{
                return f;
        }
}
float minF(float e,float f){
        if(e<=f){
                return e;
        }else{
                return f;
        }
}
long double fact(float j){
    int i=j;
    if ( (j<0)||((j-i)!=0) ) {
        printf("Problem!");
    }
    if ((i==1)||(i==0)) {
        return 1;
    } else {
        long double f=1.0;
        for(int n=2;n<=i;n++){
            f*=n;
        }
        return f;
    }
}
void f(mpfr_t f_res,int l,float J,float j,int s,float p){

    int s_term1= s+fabs(J-p);
    int s_term2= j-s+0.5*(abs(J+p)-abs(J-p));
    int s_term3= j-s-0.5*(abs(J+p)+abs(J-p));
    mpz_t s_fact,s_term1_fact,s_term2_fact,s_term3_fact;
        mpz_init(s_fact);
        mpz_fac_ui(s_fact,s);
        mpz_init(s_term1_fact);
        mpz_fac_ui(s_term1_fact,s_term1);
        mpz_init(s_term2_fact);
        mpz_fac_ui(s_term2_fact,s_term2);
        mpz_init(s_term3_fact);
        mpz_fac_ui(s_term3_fact,s_term3);
    mpfr_t norm_s;
        mpfr_init2(norm_s,l);
        mpfr_set_z(norm_s,s_fact,MPFR_RNDN);
            mpz_clear(s_fact);
        mpfr_mul_z(norm_s,norm_s,s_term1_fact,MPFR_RNDN);
            mpz_clear(s_term1_fact);
        mpfr_mul_z(norm_s,norm_s,s_term2_fact,MPFR_RNDN);
            mpz_clear(s_term2_fact);
        mpfr_mul_z(norm_s,norm_s,s_term3_fact,MPFR_RNDN);
            mpz_clear(s_term3_fact);
    mpfr_ui_div(f_res,1,norm_s,MPFR_RNDN);
        mpfr_clear(norm_s);

}
void alpha(mpc_t alpha_res,int l,float J,float rho,float j2,float j1,int s2,int s1,float p,int m){

    mpc_t z,z_temp;
        mpc_init2(z,l);
        mpc_init2(z_temp,l);
        mpc_set_ui(z,1,MPC_RNDNN);
    for(int k=0;k<=abs(J-p)+s1+s2;k++){
        mpc_set_d_d(z_temp,j1-k-m,-rho,MPC_RNDNN);
        mpc_mul(z,z,z_temp,MPC_RNDNN);
    }
        mpc_clear(z_temp);
    mpc_ui_div(z,1,z,MPC_RNDNN);

    int Lm= j1+j2-abs(J-p)-s1-s2-m;
    mpz_t m_fact,Lm_fact,mfact_Lmfact;
        mpz_init(m_fact);
        mpz_fac_ui(m_fact,m);
        mpz_init(Lm_fact);
        mpz_fac_ui(Lm_fact,Lm);
        mpz_init(mfact_Lmfact);
        mpz_mul(mfact_Lmfact,m_fact,Lm_fact);
            mpz_clear(m_fact);
            mpz_clear(Lm_fact);
    mpfr_t z_norm;
        mpfr_init2(z_norm,l);
        mpfr_set_z(z_norm,mfact_Lmfact,MPFR_RNDN);
            mpz_clear(mfact_Lmfact);
        mpfr_si_div(z_norm,pow(-1,m),z_norm,MPFR_RNDN);
    mpc_mul_fr(alpha_res,z,z_norm,MPC_RNDNN);
        mpc_clear(z);
        mpfr_clear(z_norm);
}
void alpha_0(mpfr_t alpha_res,int l,float J,float j2,float j1,int s2,int s1,float p,int m){

    int L1= abs(J-p)+s1+s2;

    int j_m= j1-m;    //in this peculiar case, j1 (and j2) is integer ! As m.
    int jLm= L1-j_m;

    mpfr_t z_fr;
        mpfr_init2(z_fr,l);
    if((0<=j_m)&&(j_m<=L1)){
        mpfr_t temp_fr;
            mpfr_init2(temp_fr,l);

        mpfr_t sum_jLm;
            mpfr_init2(sum_jLm,l);
            mpfr_set_ui(sum_jLm,0,MPFR_RNDN);
            for(int k=1;k<=jLm;k++){
                mpfr_set_ui(temp_fr,k,MPFR_RNDN);
                mpfr_ui_div(temp_fr,1,temp_fr,MPFR_RNDN);
                mpfr_add(sum_jLm,temp_fr,sum_jLm,MPFR_RNDN);
            }
        mpfr_t sum_j_m;
            mpfr_init2(sum_j_m,l);
            mpfr_set_ui(sum_j_m,0,MPFR_RNDN);
            for(int k=1;k<=j_m;k++){
                mpfr_set_ui(temp_fr,k,MPFR_RNDN);
                mpfr_ui_div(temp_fr,1,temp_fr,MPFR_RNDN);
                mpfr_add(sum_j_m,temp_fr,sum_j_m,MPFR_RNDN);
            }
        mpfr_sub(temp_fr,sum_jLm,sum_j_m,MPFR_RNDN);
            mpfr_clear(sum_jLm);
            mpfr_clear(sum_j_m);

        mpz_t j_m_fact,jLm_fact,jLmfact_jmfact;
            mpz_init(j_m_fact);
            mpz_fac_ui(j_m_fact,j_m);
            mpz_init(jLm_fact);
            mpz_fac_ui(jLm_fact,jLm);
            mpz_init(jLmfact_jmfact);
            mpz_mul(jLmfact_jmfact,j_m_fact,jLm_fact);
                mpz_clear(j_m_fact);
                mpz_clear(jLm_fact);

        mpz_mul_si(jLmfact_jmfact,jLmfact_jmfact,pow(-1,jLm));

        mpfr_div_z(z_fr,temp_fr,jLmfact_jmfact,MPFR_RNDN);
            mpfr_clear(temp_fr);
            mpz_clear(jLmfact_jmfact);
    }
    if(j_m<0){
        mpfr_t temp_fr;
            mpfr_init2(temp_fr,l);

        mpz_t j_m_fact,jLm_fact;
            mpz_init(j_m_fact);
            mpz_fac_ui(j_m_fact,-(j_m+1));
            mpz_init(jLm_fact);
            mpz_fac_ui(jLm_fact,jLm);

        mpfr_set_z(temp_fr,j_m_fact,MPFR_RNDN);
            mpz_clear(j_m_fact);
        mpfr_mul_si(temp_fr,temp_fr,-pow(-1,L1),MPFR_RNDN);

        mpfr_div_z(z_fr,temp_fr,jLm_fact,MPFR_RNDN);
            mpfr_clear(temp_fr);
            mpz_clear(jLm_fact);
    }
    if(j_m>L1){
        mpfr_t temp_fr;
            mpfr_init2(temp_fr,l);

        mpz_t j_m_fact,jLm_fact;
            mpz_init(j_m_fact);
            mpz_fac_ui(j_m_fact,j_m);
            mpz_init(jLm_fact);
            mpz_fac_ui(jLm_fact,-(jLm+1));

        mpfr_set_z(temp_fr,jLm_fact,MPFR_RNDN);
            mpz_clear(jLm_fact);

        mpfr_div_z(z_fr,temp_fr,j_m_fact,MPFR_RNDN);
            mpfr_clear(temp_fr);
            mpz_clear(j_m_fact);
    }

    int Lm= j1+j2-L1-m;
    mpz_t m_fact,Lm_fact,mfact_Lmfact;
        mpz_init(m_fact);
        mpz_fac_ui(m_fact,m);
        mpz_init(Lm_fact);
        mpz_fac_ui(Lm_fact,Lm);
        mpz_init(mfact_Lmfact);
        mpz_mul(mfact_Lmfact,m_fact,Lm_fact);
            mpz_clear(m_fact);
            mpz_clear(Lm_fact);
    mpfr_t z_norm;
        mpfr_init2(z_norm,l);
        mpfr_set_z(z_norm,mfact_Lmfact,MPFR_RNDN);
            mpz_clear(mfact_Lmfact);
        mpfr_si_div(z_norm,pow(-1,m),z_norm,MPFR_RNDN);
    mpfr_mul(alpha_res,z_fr,z_norm,MPFR_RNDN);
        mpfr_clear(z_fr);
        mpfr_clear(z_norm);
}
void sigma(mpfr_t sigma_res,int l,float J,float j2,float j1,int s2,int s1,float p,float ms){

    int s_term1= ms-s1-0.5*abs(J+p);
    int s_term2= ms-s2-0.5*abs(J-p);
    int s_term3= j2+s1+0.5*abs(J+p)-ms;
    int s_term4= j1+s2+0.5*abs(J-p)-ms;

    mpz_t sterm1_fact,sterm2_fact,sterm3_fact,sterm4_fact,sterm_prod;
        mpz_init(sterm1_fact);
        mpz_fac_ui(sterm1_fact,s_term1);
        mpz_init(sterm2_fact);
        mpz_fac_ui(sterm2_fact,s_term2);
        mpz_init(sterm3_fact);
        mpz_fac_ui(sterm3_fact,s_term3);
        mpz_init(sterm4_fact);
        mpz_fac_ui(sterm4_fact,s_term4);
        mpz_init(sterm_prod);
        mpz_mul(sterm_prod,sterm1_fact,sterm2_fact);
        mpz_mul(sterm_prod,sterm_prod,sterm3_fact);
        mpz_mul(sterm_prod,sterm_prod,sterm4_fact);
        mpz_mul_si(sterm_prod,sterm_prod,pow(-1,s1+s2));
            mpz_clear(sterm1_fact);
            mpz_clear(sterm2_fact);
            mpz_clear(sterm3_fact);
            mpz_clear(sterm4_fact);

    mpfr_t sigma_temp;
        mpfr_init2(sigma_temp,l);
        mpfr_set_z(sigma_temp,sterm_prod,MPFR_RNDN);
        mpfr_ui_div(sigma_res,1,sigma_temp,MPFR_RNDN);
            mpz_clear(sterm_prod);
}
void sqrtj(mpfr_t sqrtj_res,int l,float J,float j,float p){

    int dj= 2.0*j+1.0;
    int jpJ= j+J,jmJ= j-J;
    int jpp= j+p,jmp= j-p;
    mpz_t jpJ_fact,jmJ_fact;
        mpz_init(jpJ_fact);
        mpz_fac_ui(jpJ_fact,jpJ);
        mpz_init(jmJ_fact);
        mpz_fac_ui(jmJ_fact,jmJ);
    mpz_t jpp_fact,jmp_fact;
        mpz_init(jpp_fact);
        mpz_fac_ui(jpp_fact,jpp);
        mpz_init(jmp_fact);
        mpz_fac_ui(jmp_fact,jmp);
    mpfr_t norm_j;
        mpfr_init2(norm_j,l);
        mpfr_set_ui(norm_j,dj,MPFR_RNDN);
        mpfr_mul_z(norm_j,norm_j,jpJ_fact,MPFR_RNDN);
            mpz_clear(jpJ_fact);
        mpfr_mul_z(norm_j,norm_j,jmJ_fact,MPFR_RNDN);
            mpz_clear(jmJ_fact);
        mpfr_mul_z(norm_j,norm_j,jpp_fact,MPFR_RNDN);
            mpz_clear(jpp_fact);
        mpfr_mul_z(norm_j,norm_j,jmp_fact,MPFR_RNDN);
            mpz_clear(jmp_fact);
    mpfr_sqrt(sqrtj_res,norm_j,MPFR_RNDN);
        mpfr_clear(norm_j);

}
void Y(mpc_t Y_res,int l,float J,float rho,float j2,float j1,float p,int m){

    int s1_max= j1-0.5*(abs(J+p)+abs(J-p));
    int s2_max= j2-0.5*(abs(J+p)+abs(J-p));

    mpc_t sum_Y,sum_Y_temp;
        mpc_init2(sum_Y,l);
        mpc_init2(sum_Y_temp,l);
        mpc_set_ui_ui(sum_Y,0,0,MPC_RNDNN);
    mpfr_t f1,f2;
        mpfr_init2(f1,l);
        mpfr_init2(f2,l);
    mpz_t L1_fact,L2_fact;
        mpz_init(L1_fact);
        mpz_init(L2_fact);
    mpfr_t L1fact_L2fact;
        mpfr_init2(L1fact_L2fact,l);
    mpc_t alpha_res;
        mpc_init2(alpha_res,l);
    for(int s1=0;s1<=s1_max;s1++){
        f(f1,l,J,j1,s1,p);
        for(int s2=0;s2<=s2_max;s2++){
            f(f2,l,J,j2,s2,p);

            int L1= abs(J-p)+s1+s2;
            int L2= j1+j2-abs(J-p)-s1-s2;
            if((s2<=m)&&(m<=L2+s2)){
                mpz_fac_ui(L1_fact,L1);
                mpz_fac_ui(L2_fact,L2);
                mpfr_set_si(L1fact_L2fact,pow(-1,s1+s2),MPFR_RNDN);
                mpfr_mul_z(L1fact_L2fact,L1fact_L2fact,L1_fact,MPFR_RNDN);
                mpfr_mul_z(L1fact_L2fact,L1fact_L2fact,L2_fact,MPFR_RNDN);
                mpfr_mul(L1fact_L2fact,L1fact_L2fact,f1,MPFR_RNDN);
                mpfr_mul(L1fact_L2fact,L1fact_L2fact,f2,MPFR_RNDN);
                alpha(alpha_res,l,J,rho,j2,j1,s2,s1,p,m-s2);
                mpc_mul_fr(sum_Y_temp,alpha_res,L1fact_L2fact,MPC_RNDNN);
                mpc_add(sum_Y,sum_Y,sum_Y_temp,MPC_RNDNN);
            }
        }
    }
        mpz_clear(L1_fact);
        mpz_clear(L2_fact);
        mpfr_clear(f1);
        mpfr_clear(f2);
        mpfr_clear(L1fact_L2fact);
        mpc_clear(alpha_res);
        mpc_clear(sum_Y_temp);

    mpfr_t sqrtj1,sqrtj2,normj;
        mpfr_init2(sqrtj1,l);
        sqrtj(sqrtj1,l,J,j1,p);
        mpfr_init2(sqrtj2,l);
        sqrtj(sqrtj2,l,J,j2,p);
        mpfr_init2(normj,l);
        mpfr_mul(normj,sqrtj1,sqrtj2,MPFR_RNDN);
            mpfr_clear(sqrtj1);
            mpfr_clear(sqrtj2);

    mpc_mul_fr(Y_res,sum_Y,normj,MPC_RNDNN);
        mpc_clear(sum_Y);
        mpfr_clear(normj);
}
void Z(mpfr_t Z_res,int l,float J,float j2,float j1,float p,int m){

    int s1_max= j1-0.5*(abs(J+p)+abs(J-p));
    int s2_max= j2-0.5*(abs(J+p)+abs(J-p));

    mpfr_t sum_Z,sum_Z_temp;
        mpfr_init2(sum_Z,l);
        mpfr_init2(sum_Z_temp,l);
        mpfr_set_ui(sum_Z,0,MPFR_RNDN);
    mpfr_t f1,f2;
        mpfr_init2(f1,l);
        mpfr_init2(f2,l);
    mpz_t L1_fact,L2_fact;
        mpz_init(L1_fact);
        mpz_init(L2_fact);
    mpfr_t L1fact_L2fact;
        mpfr_init2(L1fact_L2fact,l);
    mpfr_t alpha_res;
        mpfr_init2(alpha_res,l);
    for(int s1=0;s1<=s1_max;s1++){
        f(f1,l,J,j1,s1,p);
        for(int s2=0;s2<=s2_max;s2++){
            f(f2,l,J,j2,s2,p);

            int L1= abs(J-p)+s1+s2;
            int L2= j1+j2-abs(J-p)-s1-s2;
            if((s2<=m)&&(m<=L2+s2)){
                mpz_fac_ui(L1_fact,L1);
                mpz_fac_ui(L2_fact,L2);
                mpfr_set_si(L1fact_L2fact,pow(-1,s1+s2),MPFR_RNDN);
                mpfr_mul_z(L1fact_L2fact,L1fact_L2fact,L1_fact,MPFR_RNDN);
                mpfr_mul_z(L1fact_L2fact,L1fact_L2fact,L2_fact,MPFR_RNDN);
                mpfr_mul(L1fact_L2fact,L1fact_L2fact,f1,MPFR_RNDN);
                mpfr_mul(L1fact_L2fact,L1fact_L2fact,f2,MPFR_RNDN);
                alpha_0(alpha_res,l,J,j2,j1,s2,s1,p,m-s2);
                mpfr_mul(sum_Z_temp,alpha_res,L1fact_L2fact,MPFR_RNDN);
                mpfr_add(sum_Z,sum_Z,sum_Z_temp,MPFR_RNDN);
            }
        }
    }
        mpz_clear(L1_fact);
        mpz_clear(L2_fact);
        mpfr_clear(f1);
        mpfr_clear(f2);
        mpfr_clear(L1fact_L2fact);
        mpfr_clear(alpha_res);
        mpfr_clear(sum_Z_temp);

    mpfr_t sqrtj1,sqrtj2,normj;
        mpfr_init2(sqrtj1,l);
        sqrtj(sqrtj1,l,J,j1,p);
        mpfr_init2(sqrtj2,l);
        sqrtj(sqrtj2,l,J,j2,p);
        mpfr_init2(normj,l);
        mpfr_mul(normj,sqrtj1,sqrtj2,MPFR_RNDN);
            mpfr_clear(sqrtj1);
            mpfr_clear(sqrtj2);

    mpfr_mul(Z_res,sum_Z,normj,MPFR_RNDN);
        mpfr_clear(sum_Z);
        mpfr_clear(normj);
}
void W(mpfr_t W_res,int l,float J,float j2,float j1,float p,float ms){

    int s1_max= j1-0.5*(abs(J+p)+abs(J-p));
    int s2_max= j2-0.5*(abs(J+p)+abs(J-p));

    mpfr_t sum_W,sum_W_temp;
        mpfr_init2(sum_W,l);
        mpfr_init2(sum_W_temp,l);
        mpfr_set_ui(sum_W,0,MPFR_RNDN);
    mpfr_t f1,f2;
        mpfr_init2(f1,l);
        mpfr_init2(f2,l);
    mpz_t L1_fact,L2_fact;
        mpz_init(L1_fact);
        mpz_init(L2_fact);
    mpfr_t L1fact_L2fact;
        mpfr_init2(L1fact_L2fact,l);
    mpfr_t sigma_res;
        mpfr_init2(sigma_res,l);
    for(int s1=0;s1<=s1_max;s1++){
        f(f1,l,J,j1,s1,-p);
        for(int s2=0;s2<=s2_max;s2++){
            f(f2,l,J,j2,s2,p);

            float ms_min1= s1+0.5*abs(J+p);
            float ms_min2= s2+0.5*abs(J-p);
            float ms_min= maxF(ms_min1,ms_min2);

            float ms_max1= j1+s2+0.5*abs(J-p);
            float ms_max2= j2+s1+0.5*abs(J+p);
            float ms_max= minF(ms_max1,ms_max2);

            int L1s= ms_max1-ms_min1;
            int L2s= ms_max2-ms_min2;

            if((ms_min<=ms)&&(ms<=ms_max)){
                mpz_fac_ui(L1_fact,L1s);
                mpz_fac_ui(L2_fact,L2s);
                mpfr_set_si(L1fact_L2fact,pow(-1,s1+s2),MPFR_RNDN);
                mpfr_mul_z(L1fact_L2fact,L1fact_L2fact,L1_fact,MPFR_RNDN);
                mpfr_mul_z(L1fact_L2fact,L1fact_L2fact,L2_fact,MPFR_RNDN);
                mpfr_mul(L1fact_L2fact,L1fact_L2fact,f1,MPFR_RNDN);
                mpfr_mul(L1fact_L2fact,L1fact_L2fact,f2,MPFR_RNDN);
                sigma(sigma_res,l,J,j2,j1,s2,s1,p,ms);
                mpfr_mul(sum_W_temp,sigma_res,L1fact_L2fact,MPFR_RNDN);
                mpfr_add(sum_W,sum_W,sum_W_temp,MPFR_RNDN);
            }
        }
    }
        mpz_clear(L1_fact);
        mpz_clear(L2_fact);
        mpfr_clear(f1);
        mpfr_clear(f2);
        mpfr_clear(L1fact_L2fact);
        mpfr_clear(sigma_res);
        mpfr_clear(sum_W_temp);

    mpfr_t sqrtj1,sqrtj2,normj;
        mpfr_init2(sqrtj1,l);
        sqrtj(sqrtj1,l,J,j1,p);
        mpfr_init2(sqrtj2,l);
        sqrtj(sqrtj2,l,J,j2,p);
        mpfr_init2(normj,l);
        mpfr_mul(normj,sqrtj1,sqrtj2,MPFR_RNDN);
            mpfr_clear(sqrtj1);
            mpfr_clear(sqrtj2);

    mpfr_mul(W_res,sum_W,normj,MPFR_RNDN);
        mpfr_clear(sum_W);
        mpfr_clear(normj);
}
long double complex D(mpc_t Ym[],mpc_t Yn[],int l,float J,float rho,float j2,float j1,float p,long double X){

    mpfr_t X_fr;
        mpfr_init2(X_fr,l);
        mpfr_set_ld(X_fr,X,MPFR_RNDN);

    mpc_t X_c;
        mpc_init2(X_c,l);
        mpc_set_fr(X_c,X_fr,MPC_RNDNN);

    int m_max= j1+j2-abs(J-p);
    int n_max= j1+j2-abs(J+p);

    mpc_t sum_m,km;
        mpc_init2(sum_m,l);
            mpc_set_ui_ui(sum_m,0,0,MPC_RNDNN);
        mpc_init2(km,l);
            mpc_set_d_d(km,abs(J-p)-1.0,rho,MPC_RNDNN);
    for(int m=0;m<=m_max;m++){
        mpc_add_ui(km,km,2,MPC_RNDNN);
        mpc_t sum_m_temp;
            mpc_init2(sum_m_temp,l);
            mpc_pow(sum_m_temp,X_c,km,MPC_RNDNN);
            mpc_mul(sum_m_temp,sum_m_temp,Ym[m],MPC_RNDNN);
        mpc_add(sum_m,sum_m,sum_m_temp,MPC_RNDNN);
            mpc_clear(sum_m_temp);
    }
        mpc_clear(km);

    mpc_t sum_n,kn;
        mpc_init2(sum_n,l);
            mpc_set_ui_ui(sum_n,0,0,MPC_RNDNN);
        mpc_init2(kn,l);
            mpc_set_d_d(kn,abs(J+p)-1.0,-rho,MPC_RNDNN);
    for(int n=0;n<=n_max;n++){
        mpc_add_ui(kn,kn,2,MPC_RNDNN);
        mpc_t sum_n_temp;
            mpc_init2(sum_n_temp,l);
            mpc_pow(sum_n_temp,X_c,kn,MPC_RNDNN);
            mpc_mul(sum_n_temp,sum_n_temp,Yn[n],MPC_RNDNN);
        mpc_add(sum_n,sum_n,sum_n_temp,MPC_RNDNN);
            mpc_clear(sum_n_temp);
    }
        mpc_clear(kn);
    mpc_mul_si(sum_n,sum_n,pow(-1,(int)(j2-j1)),MPC_RNDNN);

    mpc_clear(X_c);

    mpc_t sum_mn;
        mpc_init2(sum_mn,l);
        mpc_add(sum_mn,sum_m,sum_n,MPC_RNDNN);
            mpc_clear(sum_m);
            mpc_clear(sum_n);

    mpfr_t exp_term;
        mpfr_init2(exp_term,l);
        mpfr_mul(exp_term,X_fr,X_fr,MPFR_RNDN);
            mpfr_clear(X_fr);
        mpfr_ui_sub(exp_term,1,exp_term,MPFR_RNDN);
        mpfr_pow_ui(exp_term,exp_term,(int)(j2+j1)+1,MPFR_RNDN);

    mpc_div_fr(sum_mn,sum_mn,exp_term,MPC_RNDNN);
        mpfr_clear(exp_term);

    long double complex D_res=mpc_get_ldc(sum_mn,MPC_RNDNN);
        mpc_clear(sum_mn);
    return D_res;
}
long double complex D_0(mpfr_t Zm[],mpfr_t Zn[],mpfr_t W[],int l,float J,float j2,float j1,float p,long double X){

    /*if(X==0){
        return 0;
    }
    if(X==1){
        if(j1==j2){return 1;}else{return 0;}
    }*/ //not necessary for the integral

    mpfr_t X_fr;
        mpfr_init2(X_fr,l);
        mpfr_set_ld(X_fr,X,MPFR_RNDN);

    int m_max= j1+j2-abs(J-p);
    int n_max= j1+j2-abs(J+p);

    float ms_min= 0.5*maxF(abs(J-p),abs(J+p));
    float ms_max= j1+j2-ms_min;
    int ms_bound= ms_max-ms_min;

    mpfr_t sum_m,km;
        mpfr_init2(sum_m,l);
            mpfr_set_ui(sum_m,0,MPFR_RNDN);
        mpfr_init2(km,l);
            mpfr_set_d(km,abs(J-p)-1.0,MPFR_RNDN);
    for(int m=0;m<=m_max;m++){
        mpfr_add_ui(km,km,2,MPFR_RNDN);
        mpfr_t sum_m_temp;
            mpfr_init2(sum_m_temp,l);
            mpfr_pow(sum_m_temp,X_fr,km,MPFR_RNDN);
            mpfr_mul(sum_m_temp,sum_m_temp,Zm[m],MPFR_RNDN);
        mpfr_add(sum_m,sum_m,sum_m_temp,MPFR_RNDN);
            mpfr_clear(sum_m_temp);
    }
        mpfr_clear(km);

    mpfr_t sum_n,kn;
        mpfr_init2(sum_n,l);
            mpfr_set_ui(sum_n,0,MPFR_RNDN);
        mpfr_init2(kn,l);
            mpfr_set_d(kn,abs(J+p)-1.0,MPFR_RNDN);
    for(int n=0;n<=n_max;n++){
        mpfr_add_ui(kn,kn,2,MPFR_RNDN);
        mpfr_t sum_n_temp;
            mpfr_init2(sum_n_temp,l);
            mpfr_pow(sum_n_temp,X_fr,kn,MPFR_RNDN);
            mpfr_mul(sum_n_temp,sum_n_temp,Zn[n],MPFR_RNDN);
        mpfr_add(sum_n,sum_n,sum_n_temp,MPFR_RNDN);
            mpfr_clear(sum_n_temp);
    }
        mpfr_clear(kn);
    mpfr_mul_si(sum_n,sum_n,pow(-1,(int)(j2-j1)),MPFR_RNDN);

    mpfr_t sum_ms,kms;
        mpfr_init2(sum_ms,l);
            mpfr_set_ui(sum_ms,0,MPFR_RNDN);
        mpfr_init2(kms,l);
            mpfr_set_d(kms,2.0*ms_min-1.0,MPFR_RNDN);
    for(int ms=0;ms<=ms_bound;ms++){
        mpfr_add_ui(kms,kms,2,MPFR_RNDN);
        mpfr_t sum_ms_temp;
            mpfr_init2(sum_ms_temp,l);
            mpfr_pow(sum_ms_temp,X_fr,kms,MPFR_RNDN);
            mpfr_mul(sum_ms_temp,sum_ms_temp,W[ms],MPFR_RNDN);
        mpfr_add(sum_ms,sum_ms,sum_ms_temp,MPFR_RNDN);
            mpfr_clear(sum_ms_temp);
    }
        mpfr_clear(kms);
    mpfr_mul_si(sum_ms,sum_ms,2*pow(-1,(int)(j1+J+p)),MPFR_RNDN);
        mpfr_t r_fr;
            mpfr_init2(r_fr,l);
            mpfr_log(r_fr,X_fr,MPFR_RNDN);
            mpfr_mul_si(r_fr,r_fr,-1,MPFR_RNDN);
    mpfr_mul(sum_ms,sum_ms,r_fr,MPFR_RNDN);
        mpfr_clear(r_fr);

    mpfr_t sum_mn_ms;
        mpfr_init2(sum_mn_ms,l);
        mpfr_add(sum_mn_ms,sum_m,sum_n,MPFR_RNDN);
        mpfr_add(sum_mn_ms,sum_mn_ms,sum_ms,MPFR_RNDN);
            mpfr_clear(sum_m);
            mpfr_clear(sum_n);
            mpfr_clear(sum_ms);

    mpfr_t exp_term;
        mpfr_init2(exp_term,l);
        mpfr_mul(exp_term,X_fr,X_fr,MPFR_RNDN);
            mpfr_clear(X_fr);
        mpfr_ui_sub(exp_term,1,exp_term,MPFR_RNDN);
        mpfr_pow_ui(exp_term,exp_term,(int)(j2+j1)+1,MPFR_RNDN);

    mpfr_div(sum_mn_ms,sum_mn_ms,exp_term,MPFR_RNDN);
        mpfr_clear(exp_term);

    long double D_res=mpfr_get_ld(sum_mn_ms,MPFR_RNDN);
        mpfr_clear(sum_mn_ms);
    return D_res;
}
long double mesure_x(long double X){
    return powl(1.0-X*X,2)/(X*X*X)/(16.0L*M_PI);
    //return powl(1.0-X*X,2)/(X*X*X);
}
void I4_Interval(long double complex *I4_inf,long double complex *I4_sup,int N,long double complex *D1,long double complex *D2,long double complex *D3,long double complex *D4){
    long double R_min=0,R_max=0;
    long double I_min=0,I_max=0;

    long double dN=1.0L/N;

    //for n=0
        long double R_term1= 0;
        long double R_term2= mesure_x(dN)*creal(D1[1]*D2[1]*D3[1]*D4[1]);
        if(R_term1<R_term2){
            R_min+= R_term1;
            R_max+= R_term2;
        }else{
            R_min+= R_term2;
            R_max+= R_term1;
        }
        long double I_term1= 0;
        long double I_term2= mesure_x(dN)*cimag(D1[1]*D2[1]*D3[1]*D4[1]);
        if(I_term1<I_term2){
            I_min+= I_term1;
            I_max+= I_term2;
        }else{
            I_min+= I_term2;
            I_max+= I_term1;
        }

    //for n=1 to n=N-2
    for(int n=1;n<N-1;n++){
        int n2=n+1;
        long double X=n*dN;
        long double  X2=n2*dN;
        R_term1= mesure_x(X)*creal(D1[n]*D2[n]*D3[n]*D4[n]);
        R_term2= mesure_x(X2)*creal(D1[n2]*D2[n2]*D3[n2]*D4[n2]);
        if(R_term1<R_term2){
            R_min+= R_term1;
            R_max+= R_term2;
        }else{
            R_min+= R_term2;
            R_max+= R_term1;
        }
        I_term1= mesure_x(X)*cimag(D1[n]*D2[n]*D3[n]*D4[n]);
        I_term2= mesure_x(X2)*cimag(D1[n2]*D2[n2]*D3[n2]*D4[n2]);
        if(I_term1<I_term2){
            I_min+= I_term1;
            I_max+= I_term2;
        }else{
            I_min+= I_term2;
            I_max+= I_term1;
        }
    }

    //for n=N-1
        R_term1= mesure_x(1-dN)*creal(D1[N-1]*D2[N-1]*D3[N-1]*D4[N-1]);
        R_term2= 0;
        if(R_term1<R_term2){
            R_min+= R_term1;
            R_max+= R_term2;
        }else{
            R_min+= R_term2;
            R_max+= R_term1;
        }
        I_term1= mesure_x(1-dN)*cimag(D1[N-1]*D2[N-1]*D3[N-1]*D4[N-1]);
        I_term2= 0;
        if(I_term1<I_term2){
            I_min+= I_term1;
            I_max+= I_term2;
        }else{
            I_min+= I_term2;
            I_max+= I_term1;
        }

        *I4_inf= dN*(R_min+I*I_min);
        *I4_sup= dN*(R_max+I*I_max);

}
long double intertwiner_j(float Ji,long double trois_j_part1,long double trois_j_part2,float m1,float m2,float m3,float m4){
    long double result=0;

    float M12=-(m1+m2);
    float M34= m3+m4;
    if( (fabs(M12)<=Ji) && (M12==M34) ){
        //result= pow(-1.0,(int)(Ji-M12))*sqrtl(2.0*Ji+1.0)*trois_j_part1*trois_j_part2;   // <--- with normalization sqrt(2*K+1)
        result= pow(-1.0,(int)(Ji-M12))*trois_j_part1*trois_j_part2;   // <--- without the normalization sqrt(2*K+1)
    }
    return result;
}
long double complex B4_Interval(long double complex *B4_inf,long double complex *B4_sup,
                                int Ji2_int,float Ji2_min,int Ji1_int,float Ji1_min,
                                float j1,float j2,float j3,float j4,
                                long double ***trois_j2_part1,long double ***trois_j2_part2,
                                long double ***trois_j1_part1,long double ***trois_j1_part2,
                                long double complex ***I4_inf,long double complex ***I4_sup){

        long double complex sum_B4_inf = 0. + 0.*I;
        long double complex sum_B4_sup = 0. + 0.*I;

        float Ji2 = Ji2_int+Ji2_min;
        float Ji1 = Ji1_int+Ji1_min;
        float Ji= minF(Ji2,Ji1);

        int dj1= 2.0*j1+1;
        int dj2= 2.0*j2+1;
        int dj3= 2.0*j3+1;
        int dj4= 2.0*j4+1;

        for(int pa=0;pa<dj1;pa++){
                float p1=pa-j1;
                for(int pb=0;pb<dj2;pb++){
                        float p2=pb-j2;
                        for(int pc=0;pc<dj3;pc++){
                                float p3=pc-j3;
                                float p4=-p1-p2-p3;
                                if(fabs(p4)<=j4){
                                        int pd=p4+j4;

                                        float M12=-(p1+p2);
                                        float M34= p3+p4;

                                        if(fabs(M12)<=Ji){
                                                long double inter1= pow(-1.0,(int)(Ji1-M12))
                                                *trois_j1_part1[Ji1_int][pa][pb]*trois_j1_part2[Ji1_int][pc][pd];
                                                long double inter2= pow(-1.0,(int)(Ji2-M12))
                                                *trois_j2_part1[Ji2_int][pa][pb]*trois_j2_part2[Ji2_int][pc][pd];

                                                long double inter_temp= inter2*inter1;
                                                if(inter_temp>=0){
                                                        sum_B4_sup+= inter_temp*I4_sup[pa][pb][pc];
                                                        sum_B4_inf+= inter_temp*I4_inf[pa][pb][pc];
                                                }else{
                                                        sum_B4_sup+= inter_temp*I4_inf[pa][pb][pc];
                                                        sum_B4_inf+= inter_temp*I4_sup[pa][pb][pc];
                                                }
                                        }
                                }
                        }
                }
        }

        *B4_inf= sum_B4_inf;
        *B4_sup= sum_B4_sup;
}
long double complex ComplexGamma(double complex z){

        long double zr = creal(z);
        long double zi = cimag(z);

        gsl_sf_result *lnr = malloc(2*sizeof(long double));
        gsl_sf_result *arg = malloc(2*sizeof(long double));
        gsl_sf_lngamma_complex_e(zr, zi, lnr, arg);

        long double complex ComplexGammaValue;

        long double module = exp(lnr ->val);
        long double argument = arg->val;

        free(lnr);
        free(arg);

        ComplexGammaValue = module*cexp(I*argument);

        return ComplexGammaValue;

}
long double complex dPhase (float j, float l, float rho){

        long double complex Phase1 = ComplexGamma(j+rho*I+1)/cabs(ComplexGamma(j+rho*I+1));
        long double complex Phase2 = ComplexGamma(l-rho*I+1)/cabs(ComplexGamma(l+rho*I+1));

        long double complex Phase = cpow(-1,-(j-l)/2)*Phase1*Phase2;

        return Phase;


}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// End of Collet's code.
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


void B4_Hash (dspin two_j1,  dspin two_j2, dspin two_j3,
                            dspin two_j4,  dspin two_j5, dspin two_j6,
                            dspin two_j7,  dspin two_j8, dspin two_j9,
                            dspin two_j10, dspin two_Dl,
                            float Immirzi) {

    //////////////////// Parameters for Boosters ////////////////////

    dspin  two_k1, two_k2, two_k3, two_k4, two_k5,
                 two_k6, two_k7, two_k8, two_k9, two_k10;

    float two_rho1, two_rho2, two_rho3, two_rho4, two_rho5,
                two_rho6, two_rho7, two_rho8, two_rho9, two_rho10;


    two_k1 = two_j1; two_k2 = two_j2; two_k3 = two_j3; two_k4 = two_j4; two_k5 = two_j5;
    two_k6 = two_j6; two_k7 = two_j7; two_k8 = two_j8; two_k9 = two_j9; two_k10 = two_j10;

    two_rho1 = Immirzi*(float)two_j1; two_rho2 = Immirzi*(float)two_j2; two_rho3 = Immirzi*(float)two_j3;
    two_rho4 = Immirzi*(float)two_j4; two_rho5 = Immirzi*(float)two_j5; two_rho6 = Immirzi*(float)two_j6;
    two_rho7 = Immirzi*(float)two_j7; two_rho8 = Immirzi*(float)two_j8; two_rho9 = Immirzi*(float)two_j9;
    two_rho10 = Immirzi*(float)two_j10;

    // Folder check.
    check_data_4simplex(Immirzi);

    //////////////////// Hash Table initialization ///////////////////
    khash_t(HashTableBooster) *h = NULL;

    //////////////////// build paths ///////////////////

    struct data_folders* fd = get_data_folders(Immirzi);

    char filename[1024];
    char path_boost[1024];
    char path_boost_dli[1024];
    char path_key_boost[1024];
    char path_key_boost_dli[1024];

    // build path for boosts
    sprintf(filename, "%i.%i.%i.%i_%i.boost", two_j5, two_j6, two_j9, two_j10, two_Dl);
    strcpy(path_boost, fd->foursimp_imm_boost);
    strcat(path_boost, filename);

    // build path for boost keys
    sprintf(filename, "%i.%i.%i.%i_%i.kboost", two_j5, two_j6, two_j9, two_j10, two_Dl);
    strcpy(path_key_boost, fd->foursimp_aux_imm_boost);
    strcat(path_key_boost, filename);

    //////////////////// Check if data already exists and load previously computed tables ////////////////////
    int data_yes = 0;
    int data_check = 0;
    for (int two_i = two_Dl; two_i >= 0; two_i -= 2) {

        sprintf(filename, "%i.%i.%i.%i_%i.boost", two_j5, two_j6, two_j9, two_j10, two_i);
        strcpy(path_boost_dli, fd->foursimp_imm_boost);
        strcat(path_boost_dli, filename);

        if(file_exist (path_boost_dli) != 0) {

            data_yes = 1;
            data_check = two_i;
            h = kh_load(HashTableBooster, path_boost_dli);
            break;

        }

    }

    if (h == NULL) {
        h = kh_init(HashTableBooster);
    }

    // initialize a table with keys to know
    // what boosters have been already computed

    khash_t(HashTableKeyBooster) *h2 = NULL;

    for (int two_i = two_Dl; two_i >= 0; two_i -= 2) {

        sprintf(filename, "%i.%i.%i.%i_%i.kboost", two_j5, two_j6, two_j9, two_j10, two_i);
        strcpy(path_key_boost_dli, fd->foursimp_aux_imm_boost);
        strcat(path_key_boost_dli, filename);

        if (file_exist(path_key_boost_dli) != 0) {
            h2 = kh_load(HashTableKeyBooster, path_key_boost_dli);
            break;
        }

    }

    if (h2 == NULL) {
        h2 = kh_init(HashTableKeyBooster);
    }

    //////////////////// Check if I have the boosters. If not compute ////////////////////

    HashTableKeyBooster_key_t key = {two_j1, two_j2, two_j3, two_j4, two_j7, two_j8};

    if (kh_get(HashTableKeyBooster, h2, key) != kh_end(h2) &&
            data_yes == 1 && data_check == two_Dl) {

        kh_destroy(HashTableKeyBooster, h2);
        kh_destroy(HashTableBooster, h);
        return;
        
    }

    int ret;

    // put auxiliary key for boosters
    khint_t k = kh_put(HashTableKeyBooster, h2, key, &ret);

    // check for errors
    if (ret == -1) {
        error("error inserting key into aux Booster hash table");
    }

    //////////////////////////////////////////////
    // Boosters Hashing
    //////////////////////////////////////////////

    // Intertwiners range
    int two_Ji1_minA = (int)max(abs(two_j2-two_j3),abs(two_j5-two_j4));
    int two_Ji1_maxA = (int)min(two_j2+two_j3,two_j5+two_j4);
    int two_Ji1_boundA = two_Ji1_maxA-two_Ji1_minA;

    // All needed l's combinations
    for (dspin two_l2 = two_j2; two_l2 <= two_j2+two_Dl; two_l2 += 2) {
    for (dspin two_l3 = two_j3; two_l3 <= two_j3+two_Dl; two_l3 += 2) {
    for (dspin two_l4 = two_j4; two_l4 <= two_j4+two_Dl; two_l4 += 2) {

        int two_Ji2_minA = (int)max(abs(two_l2-two_l3),abs(two_j5-two_l4));
        int two_Ji2_maxA = (int)min(two_l2+two_l3,two_j5+two_l4);
        int two_Ji2_boundA = two_Ji2_maxA-two_Ji2_minA;

        HashTableBooster_key_t keyCheckA = {two_j2, two_j3, two_j5, two_j4,
                                            two_l2,two_l3,two_j5,two_l4,
                                            two_Ji1_minA,two_Ji2_minA};

        // Check we haven't already the function
        if (kh_get(HashTableBooster, h, keyCheckA) == kh_end(h)){

            // Compute the function for every intertwiner
            long double  **B4_moyA = B4Function(two_k2, two_k3, two_k5, two_k4,
                                                two_rho2, two_rho3, two_rho5, two_rho4,
                                                two_j2, two_j3, two_j5, two_j4,
                                                two_l2, two_l3, two_j5, two_l4);

            // Save all intertwiners combinations
            for(int two_Ji1A = 0; two_Ji1A <= two_Ji1_boundA; two_Ji1A += 2){
                for(int two_Ji2A = 0; two_Ji2A <= two_Ji2_boundA; two_Ji2A += 2){

                    int retA;
                    HashTableBooster_key_t keyA = {two_j2, two_j3, two_j5, two_j4,
                                                   two_l2,two_l3,two_j5,two_l4,
                                                   two_Ji1A+two_Ji1_minA,two_Ji2A+two_Ji2_minA};

                    khint_t kA = kh_put(HashTableBooster, h, keyA, &retA);

                    // check for error in hashing the ket
                    if (retA == -1) {
                        error("error inserting key into first Booster hash table");
                    }

                    // if the key is correct put the value
                    if (retA == 1) {

                        int Ji1A = two_Ji1A/2;
                        int Ji2A = two_Ji2A/2;
                        double valB4A = (double)B4_moyA[Ji2A][Ji1A];

                        kh_val(h,kA) = valB4A;
                        kA = kh_get(HashTableBooster, h, keyA);
                    
                        // check the inserted value 

                        if (kh_val(h,kA) == 0 && valB4A != 0) {

                            // retry
                            kh_val(h,kA) = valB4A;

                            // check again and now if not fail
                            if (kh_val(h,kA) == 0 && valB4A != 0) {
                                error("error inserting value into first Booster hash table");
                            }
                        }

                    }
                }
            }

        free (B4_moyA);

        }

        } // l4
        } // l3
        } // l2

        if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

            //////////////////////////////////////////////
            // hashing REDUCIBLE case
            //////////////////////////////////////////////


            int two_Ji1_minB = (int)max(abs(two_j10-two_j7),abs(two_j1-two_j2));
            int two_Ji1_maxB = (int)min(two_j10+two_j7,two_j1+two_j2);
            int two_Ji1_boundB = two_Ji1_maxB-two_Ji1_minB;

            for (dspin two_l7 = two_j7 ; two_l7 <= two_j7+two_Dl; two_l7 += 2) {
            for (dspin two_l1 = two_j1; two_l1 <= two_j1+two_Dl; two_l1 += 2) {
            for (dspin two_l2 = two_j2; two_l2 <= two_j2+two_Dl; two_l2 += 2) {

                        int two_Ji2_minB = (int)max(abs(two_j10-two_l7),abs(two_l1-two_l2));
                        int two_Ji2_maxB = (int)min(two_j10+two_l7,two_l1+two_l2);
                        int two_Ji2_boundB = two_Ji2_maxB-two_Ji2_minB;

                        HashTableBooster_key_t keyCheckB = {two_j10, two_j7, two_j1, two_j2,
                                                            two_j10,two_l7,two_l1,two_l2,
                                                            two_Ji1_minB,two_Ji2_minB};

                        if (kh_get(HashTableBooster, h, keyCheckB) == kh_end(h)){

                            long double **B4_moyB = B4Function(two_k10, two_k7, two_k1, two_k2,
                                                               two_rho10, two_rho7, two_rho1, two_rho2,
                                                               two_j10, two_j7, two_j1, two_j2,
                                                               two_j10, two_l7, two_l1, two_l2);

                            for (int two_Ji1B = 0; two_Ji1B <= two_Ji1_boundB; two_Ji1B += 2) {
                                for (int two_Ji2B = 0; two_Ji2B <= two_Ji2_boundB; two_Ji2B += 2) {

                                    int retB;
                                    HashTableBooster_key_t keyB = {two_j10, two_j7, two_j1, two_j2,
                                                                                    two_j10,two_l7,two_l1,two_l2,
                                                                                    two_Ji1B+two_Ji1_minB,two_Ji2B+two_Ji2_minB};
                                    khint_t kB = kh_put(HashTableBooster, h, keyB, &retB);
                                    
                                    // check for error in hashing the key
                                    if (retB == -1) {
                                        error("error inserting key into second Booster hash table");
                                    }

                                    // if the key is correct put the value
                                    if (retB == 1) {

                                        int Ji1B = two_Ji1B/2;
                                        int Ji2B = two_Ji2B/2;
                                        double valB4B = (double)B4_moyB[Ji2B][Ji1B];

                                        kh_val(h,kB) = valB4B;
                                        kB = kh_get(HashTableBooster, h, keyB);
                                        
                                        // check the inserted value 
                                        if (kh_val(h,kB) == 0 && valB4B != 0) {

                                            // retry
                                            kh_val(h,kB) = valB4B;

                                            // check again and now if not fail
                                            if (kh_val(h,kB) == 0 && valB4B != 0) {
                                                error("error inserting value into second Booster hash table");
                                            }

                                        }

                                    }

                                }
                            }

                    free (B4_moyB);

                    }

            } // l2
            } // l1
            } // l7

        } else {

            //////////////////////////////////////////////
            // hashing IRREDUCIBLE case
            //////////////////////////////////////////////

            int two_Ji1_minB = (int)max(abs(two_j1-two_j10),abs(two_j7-two_j2));
            int two_Ji1_maxB = (int)min(two_j1+two_j10,two_j7+two_j2);
            int two_Ji1_boundB = two_Ji1_maxB-two_Ji1_minB;

            for (dspin two_l7 = two_j7; two_l7 <= two_j7+two_Dl; two_l7 += 2) {
            for (dspin two_l1 = two_j1; two_l1 <= two_j1+two_Dl; two_l1 += 2) {
            for (dspin two_l2 = two_j2; two_l2 <= two_j2+two_Dl; two_l2 += 2) {

                int two_Ji2_minB=(int)max(abs(two_l1-two_j10),abs(two_l7-two_l2));
                int two_Ji2_maxB=(int)min(two_l1+two_j10,two_l7+two_l2);
                int two_Ji2_boundB=two_Ji2_maxB-two_Ji2_minB;

                HashTableBooster_key_t keyCheckB = {two_j1, two_j10, two_j7, two_j2,
                                                    two_l1, two_j10, two_l7, two_l2,
                                                    two_Ji1_minB,two_Ji2_minB};

                if (kh_get(HashTableBooster, h, keyCheckB) == kh_end(h)) {

                    long double **B4_moyB = B4Function(two_k1, two_k10, two_k7, two_k2,
                                                       two_rho1, two_rho10, two_rho7, two_rho2,
                                                       two_j1, two_j10, two_j7, two_j2,
                                                       two_l1, two_j10, two_l7, two_l2);

                    for (int two_Ji1B = 0; two_Ji1B <= two_Ji1_boundB; two_Ji1B += 2) {
                        for (int two_Ji2B = 0; two_Ji2B <= two_Ji2_boundB; two_Ji2B += 2) {

                            int retB;
                            HashTableBooster_key_t keyB = {two_j1, two_j10, two_j7, two_j2,
                                                           two_l1, two_j10, two_l7, two_l2,
                                                           two_Ji1B+two_Ji1_minB,two_Ji2B+two_Ji2_minB};
                            khint_t kB = kh_put(HashTableBooster, h, keyB, &retB);
                            
                            // check for error in hashing the key
                            if (retB == -1) {
                                error("error inserting key into second Booster hash table");
                            }

                            // if the key is correct put the value
                            if (retB == 1){

                                int Ji1B = two_Ji1B/2;
                                int Ji2B = two_Ji2B/2;
                                double valB4B = (double)B4_moyB[Ji2B][Ji1B];

                                kh_val(h,kB) = valB4B;
                                kB = kh_get(HashTableBooster, h, keyB);
                                
                                // check the inserted value 
                                if (kh_val(h,kB) == 0 && valB4B != 0) {

                                    // retry
                                    kh_val(h,kB) = valB4B;

                                    // check again and now if not fail
                                    if (kh_val(h,kB) == 0 && valB4B != 0) {
                                        error("error inserting value into second Booster hash table");
                                    }
                                    
                                }              
                            }
                        }
                    }

                free (B4_moyB);

                }

            } // l2
            } // l1
            } // l7

        }

    int two_Ji1_minC = (int)max(abs(two_j9-two_j8),abs(two_j3-two_j1));
    int two_Ji1_maxC = (int)min(two_j9+two_j8,two_j3+two_j1);
    int two_Ji1_boundC = two_Ji1_maxC-two_Ji1_minC;

    for (dspin two_l8 = two_j8 ; two_l8 <= two_j8+two_Dl; two_l8 += 2) {
        for (dspin two_l3 = two_j3; two_l3 <= two_j3+two_Dl; two_l3 += 2) {
            for (dspin two_l1 = two_j1; two_l1 <= two_j1+two_Dl; two_l1 += 2) {

                int two_Ji2_minC = (int)max(abs(two_j9-two_l8),abs(two_l3-two_l1));
                int two_Ji2_maxC = (int)min(two_j9+two_l8,two_l3+two_l1);
                int two_Ji2_boundC = two_Ji2_maxC-two_Ji2_minC;

                HashTableBooster_key_t keyCheckC = {two_j9,two_j8,two_j3,two_j1,
                                                    two_j9,two_l8,two_l3,two_l1,
                                                    two_Ji1_minC,two_Ji2_minC};

                if (kh_get(HashTableBooster, h, keyCheckC) == kh_end(h)) {

                    long double  **B4_moyC = B4Function(two_k9,  two_k8,  two_k3,  two_k1,
                                                        two_rho9,  two_rho8,  two_rho3,  two_rho1,
                                                        two_j9,  two_j8,  two_j3,  two_j1,
                                                        two_j9,  two_l8,  two_l3,  two_l1);

                    for(int two_Ji1C = 0; two_Ji1C <= two_Ji1_boundC; two_Ji1C += 2) {
                        for(int two_Ji2C = 0; two_Ji2C <= two_Ji2_boundC; two_Ji2C += 2) {

                            int retC;
                            HashTableBooster_key_t keyC = {two_j9,two_j8,two_j3,two_j1,
                                                                                         two_j9,two_l8,two_l3,two_l1,
                                                                                         two_Ji1C+two_Ji1_minC,two_Ji2C+two_Ji2_minC};
                            khint_t kC = kh_put(HashTableBooster, h, keyC, &retC);

                            //check for error in hashing the key
                            if (retC == -1) {
                                error("error inserting key into third Booster hash table");
                            }

                            //if the key is correct put the value
                            if (retC == 1){

                                int Ji1C = two_Ji1C/2;
                                int Ji2C = two_Ji2C/2;
                                double valB4C = (double)B4_moyC[Ji2C][Ji1C];

                                kh_val(h,kC) = valB4C;
                                kC = kh_get(HashTableBooster, h, keyC);
                                
                                //check the inserted value 
                                if (kh_val(h,kC) == 0 && valB4C != 0) {

                                    // retry
                                    kh_val(h,kC) = valB4C;

                                    // check again and now if not fail
                                    if (kh_val(h,kC) == 0 && valB4C != 0) {
                                        error("error inserting value into third Booster hash table");
                                    }
                                }          
                            }            
                        }
                    }
                free (B4_moyC);
                }
            }
        }
    }

    if (get_coupling() == SL2CFOAM_COUPLING_REDUCIBLE) {

        //////////////////////////////////////////////
        // hashing REDUCIBLE case
        //////////////////////////////////////////////

        int two_Ji1_minD = (int)max(abs(two_j4-two_j6),abs(two_j7-two_j8));
        int two_Ji1_maxD = (int)min(two_j4+two_j6,two_j7+two_j8);
        int two_Ji1_boundD = two_Ji1_maxD-two_Ji1_minD;

        for (dspin two_l4 = two_j4 ; two_l4 <= two_j4+two_Dl; two_l4 += 2) {
            for (dspin two_l7 = two_j7; two_l7 <= two_j7+two_Dl; two_l7 += 2) {
                for (dspin two_l8 = two_j8; two_l8 <= two_j8+two_Dl; two_l8 += 2) {

                    int two_Ji2_minD = (int)max(abs(two_l4-two_j6),abs(two_l7-two_l8));
                    int two_Ji2_maxD = (int)min(two_l4+two_j6,two_l7+two_l8);
                    int two_Ji2_boundD = two_Ji2_maxD-two_Ji2_minD;

                    HashTableBooster_key_t keyCheckD = {two_j4,two_j6,two_j7,two_j8,
                                                        two_l4,two_j6,two_l7,two_l8,
                                                        two_Ji1_minD,two_Ji2_minD};

                    if (kh_get(HashTableBooster, h, keyCheckD) == kh_end(h)) {

                        long double **B4_moyD = B4Function(two_k4, two_k6, two_k7, two_k8,
                                                           two_rho4, two_rho6, two_rho7, two_rho8,
                                                           two_j4, two_j6, two_j7, two_j8,
                                                           two_l4, two_j6, two_l7, two_l8);

                        for(int two_Ji1D = 0; two_Ji1D <= two_Ji1_boundD; two_Ji1D += 2) {
                            for(int two_Ji2D = 0; two_Ji2D <= two_Ji2_boundD; two_Ji2D += 2) {

                                int retD;
                                HashTableBooster_key_t keyD = {two_j4,two_j6,two_j7,two_j8,
                                                                                two_l4,two_j6,two_l7,two_l8,
                                                                                two_Ji1D+two_Ji1_minD,two_Ji2D+two_Ji2_minD};
                                khint_t kD = kh_put(HashTableBooster, h, keyD, &retD);

                                // check for error in hashing the key
                                if (retD == -1) {
                                    error("error inserting key into fourth Booster hash table");
                                }

                                // if the key is correct put the value
                                if (retD == 1){

                                    int Ji1D = two_Ji1D/2;
                                    int Ji2D = two_Ji2D/2;
                                    double valB4D = (double)B4_moyD[Ji2D][Ji1D];

                                    kh_val(h,kD) = valB4D;
                                    kD = kh_get(HashTableBooster, h, keyD);
                                    
                                    //check the inserted value 
                                    if (kh_val(h,kD) == 0 && valB4D != 0) {

                                        // retry
                                        kh_val(h,kD) = valB4D;

                                        // check again and now if not fail
                                        if (kh_val(h,kD) == 0 && valB4D != 0) {
                                            error("error inserting value into fourth Booster hash table");
                                        }
                                    }          
                                }
                            }
                        }

                    free (B4_moyD);

                    }
                }
            }
        }

    } else {

        //////////////////////////////////////////////
        // hashing IRREDUCIBLE case
        //////////////////////////////////////////////

        int two_Ji1_minD = (int)max(abs(two_j4-two_j7),abs(two_j8-two_j6));
        int two_Ji1_maxD = (int)min(two_j4+two_j7,two_j8+two_j6);
        int two_Ji1_boundD = two_Ji1_maxD-two_Ji1_minD;

        for (dspin two_l4 = two_j4; two_l4 <= two_j4+two_Dl; two_l4 += 2) {
        for (dspin two_l7 = two_j7; two_l7 <= two_j7+two_Dl; two_l7 += 2) {
        for (dspin two_l8 = two_j8; two_l8 <= two_j8+two_Dl; two_l8 += 2) {

            int two_Ji2_minD = (int)max(abs(two_l4-two_l7),abs(two_l8-two_j6));
            int two_Ji2_maxD = (int)min(two_l4+two_l7,two_l8+two_j6);
            int two_Ji2_boundD = two_Ji2_maxD-two_Ji2_minD;

            HashTableBooster_key_t keyCheckD = {two_j4,two_j7,two_j8,two_j6,
                                                two_l4,two_l7,two_l8,two_j6,
                                                two_Ji1_minD,two_Ji2_minD};

            if (kh_get(HashTableBooster, h, keyCheckD) == kh_end(h)){

                long double **B4_moyD = B4Function(two_k4, two_k7, two_k8, two_k6,
                                                two_rho4, two_rho7, two_rho8, two_rho6,
                                                two_j4,two_j7,two_j8,two_j6,
                                                two_l4,two_l7,two_l8,two_j6);

                for(int two_Ji1D = 0; two_Ji1D <= two_Ji1_boundD; two_Ji1D += 2) {
                    for(int two_Ji2D = 0; two_Ji2D <= two_Ji2_boundD; two_Ji2D += 2) {

                        int retD;
                        HashTableBooster_key_t keyD = {two_j4,two_j7,two_j8,two_j6,
                                                    two_l4,two_l7,two_l8,two_j6,
                                                    two_Ji1D+two_Ji1_minD,two_Ji2D+two_Ji2_minD};
                        khint_t kD = kh_put(HashTableBooster, h, keyD, &retD);

                        //check for error in hashing the key
                        if (retD == -1) {
                            error("error inserting key into fourth Booster hash table");
                        }

                        //if the key is correct put the value
                        if (retD == 1){

                            int Ji1D = two_Ji1D/2;
                            int Ji2D = two_Ji2D/2;
                            double valB4D = (double)B4_moyD[Ji2D][Ji1D];

                            kh_val(h,kD) = valB4D;
                            kD = kh_get(HashTableBooster, h, keyD);
                            
                            //check the inserted value 
                            if (kh_val(h,kD) == 0 && valB4D != 0) {

                                // retry
                                kh_val(h,kD) = valB4D;

                                // check again and now if not fail
                                if (kh_val(h,kD) == 0 && valB4D != 0) {
                                    error("error inserting value into fourth Booster hash table");
                                }

                            }      

                        }
                    }
                }

            free (B4_moyD);

            }

        } // l8
        } // l7
        } // l4

    }

    // write boosts
    kh_write(HashTableKeyBooster, h2, path_key_boost);
    kh_write(HashTableBooster, h, path_boost);

    // free memory
    kh_destroy(HashTableKeyBooster, h2);
    kh_destroy(HashTableBooster, h);

}
