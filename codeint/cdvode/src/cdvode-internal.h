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
#ifndef H__CDVODE_INTERNAL__H
#define H__CDVODE_INTERNAL__H
#include "../../../codeint/codeintdeps/include/stack.h"
#include "../include/cdvode.h"

/*
  declaring the odepack fortran interfaces
 */


extern void dvode_(void (*f)(const int *neq, const double *t, const double *y, double *ydot, double* rpar, int* ipar),
                   const int *neq,
                   double *y,
                   double *t,
                   const double *tout,
                   const int *itol,
                   const double *rtol,
                   const double *atol,
                   const int *itask,
                   int *istate,
                   const int *iopt,
                   double *rwork,
                   const int *lrw,
                   int *iwork,
                   const int *liw,
                   void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd, double* rpar, int* ipar),
                   const int *mf,
                   double* rpar,
                   int* ipar);





#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#define MIN(a, b) ((a) <= (b) ? (a) : (b))

typedef enum cvode_function_e {DVODE} CDVODE_FUNCTION;
typedef int CDVODE_METHOD_FLAG;
/*
 * this captures all the input params 
 */
typedef struct dvode_params_s {
    CDVODE_FUNCTION ode_func;
    cdvode_ode_func f_func;
    int neq;
    int itol;
    union {
        double *rtol_vec;
        double rtol_val;
    } rtol;

    union {
        double *atol_vec;
        double atol_val;
    } atol;

    CDVODE_ITASK itask;
    
    union {
        CDVODE_ISTATE_IN  istate_in;
        CDVODE_ISTATE_OUT istate_out;
        int istate;
    } istate;
    
    CDVODE_OPTIONAL_INPUT_FLAG iopt;
    double* rwork;
    int* iwork;
    int lrw;
    int liw;
    cdvode_jac_func jac;
    struct {
        CDVODE_JSV jsv;
        CDVODE_METH meth;
        CDVODE_MITER miter;
    } method_params;

    CDVODE_METHOD_FLAG mf;

    double *rpar;
    int* ipar;
} cdvode_params;


double* write_to_stack(double* stack, int neq, int* index, double t_new, double* q);


#endif


