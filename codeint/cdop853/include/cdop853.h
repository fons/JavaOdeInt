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
#ifndef H__DOP853__H
#define H__DOP853__H

/*
  n : diremsion
  x : independent variable
  y : vector
  f : derivative of y
 */

typedef void (*dop853_ode_func)(const int* n, const double *x, const double *y, double *f, double* rpar, int* ipar);

typedef void (*dop853_solout)(int* nr, double* xold, double* x, double* y, int* n, double* con, int* icomp, int* nd, double* rpar, int* ipar, int* irtrn);


typedef enum dop853_idid_e {SUCCESS=1,
                            SUCCESS_INTR=2,
                            INPUT_INCONSISTENT=-1,
                            NMAX_TOO_SMALL=-2,
                            STEP_TOO_SMALL=-3,
                            STIFF_PROBLEM=-4} DOP853_RETVAL;

typedef enum dop853_itol_e { ALL_SCALAR = 1, ALL_ARRAY=2} DOP853_ITOLERANCE;
typedef enum dop853_iout_e { NEVER_CALLED=0, OUTPUT=2, DENSE_OUTPUT=3 } DOP853_IOUT; 

DOP853_RETVAL dop853_basic(double* stack, double* y, dop853_ode_func f_func, int n, double xstart, double xend, double deltax);

void dop853(
    int* n,
    void (*fcn) (const int* n, const double *x, const double *y, double *f, double* rpar, int* ipar),
    double* x,
    double* y,
    double* xend,
    double* rtol,
    double* atol,
    const int* itol,
    void (*solout)(int* nr, double* xold, double* x, double* y, int* n, double* con, int* icomp, int* nd, double* rpar, int* ipar, int* irtrn),
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
