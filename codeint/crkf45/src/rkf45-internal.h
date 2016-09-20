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
#ifndef H__RKF45_INTERNAL__H
#define H__RKF45_INTERNAL__H
#include "../../codeintdeps/include/stack.h" 
#include "../include/crkf45.h"

typedef struct rkf45_params_s {

    rkf45_ode_func f_func;
    int neq;
    double rtol;
    double atol;
    union {
        RKF45_IFLAG  iflag;
        RKF45_RETVAL retval;
    } iflag;
    double* rwork;
    int*    iwork;
    int lrw;
    int liw;
} rkf45_params;

double* write_to_stack(double* stack, int neq, int* index, double t_new, double* q);
extern void r8_rkf45_ ( void (*f)(double*, double* , double*),
                        int* neqn,
                        double* y,
                        double* t,
                        double* tout,
                        double* relerr,
                        double* abserr, 
                        int* iflag,
                        double* work,
                        int* iwork );


#endif
