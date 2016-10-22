/* 
 * 
 * https://opensource.org/licenses/BSD-3-Clause
* 
* Copyright (c) 2016, JodeInt developers
* All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
* 
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 
* 2. Redistributions in binary form must reproduce the above copyright notice, this list
* of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived
* from this software without specific prior written permission.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
* BUT NOT LIMITED TO,THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
* THE COPYRIGHT HOLDER OR CONTRIBUTORS BELIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
* THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* EOM
*/ 
#ifndef H__CRADUA5__H
#define H__CRADUA5__H

typedef void (*radau5_ode_func)(int* n, double* x, double* y, double* f, double* rpar, int* ipar);
typedef void (*radau5_jacobian)(int* n, double* x, double* y, double* dfy, int* ldfy, double* rpar, int* ipar);
typedef void (*radau5_mass)(int* n, double* am,int* lmas,double* rpar, int* ipar);
typedef void (*radau5_solout)(int* nr, double* xold, double* x, double* y, double* cont, int* lrc, int* n, double* rpar, int* ipar, int* irtrn); 

typedef enum radau5_itol_e { ALL_SCALAR = 0, ALL_ARRAY=1} RADAU5_ITOLERANCE;
typedef enum radau5_jacobian_e {INTERNAL=0, JAC_USER_PROVIDED=1 } RADAU5_JACOBIAN;
typedef enum radau5_mass_matrix_e {IDENTITY_MATRIX=0, MASS_USER_PROVIDED=1} RADAU5_MASS_MATRIX;
typedef enum radau5_iout_e { NEVER_CALLED=0, OUTPUT=1} RADAU5_IOUT;

typedef enum radau5_idid_e {SUCCESS=1,
                            SUCCESS_INTR=2,
                            INPUT_INCONSISTENT=-1,
                            NMAX_TOO_SMALL=-2,
                            STEP_TOO_SMALL=-3,
                            STIFF_PROBLEM=-4} RADAU5_RETVAL;

void radau5(
            int* n,
            void (*fcn)(int* , double* , double* , double* , double* , int* ),
            double* x,
            double* y,
            double* xend,
            double* h,
            double* rtol,
            double* atol,
            int* itol,
            void (*jac)(int* , double* , double* , double* , int* , double* , int* ),
            int* ijac,
            int* mljac,
            int* mujac,
            void (*mas)(int* , double* ,int* ,double* , int* ),
            int* imas,
            int* mlmas,
            int* mumas,
            void (*solout)(int* , double* , double* , double* , double* , int*, int*, double*, int* , int* ),
            int* iout,
            double* work,
            int* lwork,
            int* iwork,
            int* liwork,
            double* rpar,
            int* ipar,
            int* idid
            );

RADAU5_RETVAL radau5_basic(double* stack, double* y, radau5_ode_func f_func, int n, double xstart, double xfinal, double deltax);

#endif
