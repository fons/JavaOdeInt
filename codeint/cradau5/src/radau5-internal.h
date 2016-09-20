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
#ifndef H__RADAU5_INTERNAL__H
#define H__RADAU5_INTERNAL__H

#include "../../codeintdeps/include/stack.h"
#include "../include/cradau5.h"

typedef struct radau5_params_s {
    radau5_ode_func f_func;
    radau5_jacobian jac_func;
    radau5_mass     mass_func;
    radau5_solout   solout;
    int neq;
    RADAU5_ITOLERANCE itol;
    union {
        double *rtol_vec;
        double rtol_val;
    } rtol;

    union {
        double *atol_vec;
        double atol_val;
    } atol;
    int mljac;
    int mujac;
    double h;
    double* rwork;
    int*    iwork;
    int lrw;
    int liw;
    RADAU5_JACOBIAN ijac;
    RADAU5_MASS_MATRIX imas;
    int mlmas;
    int mumas;
    RADAU5_IOUT iout;
    RADAU5_RETVAL retval;
    double *rpar;
    int *ipar;
} radau5_params;

void radau5_(
    int* n,
    void (*fcn)(int* , double* , double* , double* , double* , int* ),
    double* x,
    double* y,
    double* xend,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int* , double* , double* , double** , int* , double* , int* ),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int* , double** ,int* ,double* , int* ),
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


#endif
