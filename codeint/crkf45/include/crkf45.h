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
#ifndef H__CRKF45__H
#define H__CRKF45__H

typedef void (*rkf45_ode_func)(double* t, double* y, double* yp);

typedef enum rkf45_iflag_in_e  {ONE_STEP = -1, NORMAL=1} RKF45_IFLAG;

typedef enum rkf45_retval_e {ONE_STEP_SUCCES = -2,
                             SUCCESS=2,
                             RELATIVE_ERROR_TOO_SMALL=3,
                             TOO_MANY_EVAKIUATIONS=4,
                             SOLUTION_VANISHED=5,
                             REQUESTED_ACCURACY_TOO_SMALL=6,
                             USE_ONE_STEP_MODE=7,
                             INVALID_INPUT=8} RKF45_RETVAL;

RKF45_RETVAL rkf45_basic(double* stack, double* q,  rkf45_ode_func f_func,int neq, double t0, double tf, double dt);

void rkf45 ( void (*f)(double*, double* , double*),
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
