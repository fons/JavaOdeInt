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
#ifndef H__CZVODE__H
#define H__CZVODE__H

/*
 * function pointers compatible with the odepack 
 */

typedef struct zvode_complex_s {
    double dr;
    double di;
} ZVODE_COMPLEX;

/*
 * ode function 
 */
typedef void (*czvode_ode_func)(const int *neq, const double *t_, const ZVODE_COMPLEX *y, ZVODE_COMPLEX *ydot, ZVODE_COMPLEX* rpar, int* ipar);

/*
 * jacobian 
 */
typedef void (*czvode_jac_func) (const int *neq, const double *t_, const ZVODE_COMPLEX *y, const int *ml, const int *mu, ZVODE_COMPLEX *pd, const int *nrowpd, ZVODE_COMPLEX* rpar, int* ipar);

typedef enum czvode_itol_e { ALL_SCALAR = 1, ATOL_ARRAY=2, RTOL_ARRAY=3, ALL_ARRAY=4} CZVODE_ITOLERANCE;
typedef enum code_itask_e {NORMAL = 1, ONESTEP_ONLY = 2, STOP_AT_FIRST_MESH = 3, NORMAL_TCRIT = 4, ONESTEP_TCRIT = 5} CZVODE_ITASK;
typedef enum czvode_istate_in_e {FIRST_CALL = 1, NEXT_CALL = 2, NEXT_CALL_WITH_CHANGES = 3} CZVODE_ISTATE_IN;
typedef enum czvode_istate_out_e { NOTHING_DONE           = 1,
                                  SUCCESS_DONE           = 2,
                                  MXSTEPS_EXCEEDED       = -1,
                                  TO_MUCH_ACCURACY       = -2,
                                  ILLEGAL_INPUT          = -3 ,
                                  ERROR_TEST_FAILURES    = -4,
                                  CONVERGENCE_FAILURES   = -5,
                                  ZERO_ERR_TOLERANCE     = -6} CZVODE_ISTATE_OUT;

typedef enum czvode_iopt_e {NO_OPTIONAL_INPUTS = 0, OPTIONAL_INPUTS = 1} CZVODE_OPTIONAL_INPUT_FLAG;

typedef enum czvode_jacobian_saving_strategy_e {SAVE_COPY=1 , NO_COPY=-1} CZVODE_JSV;

typedef enum czvode_method_e {ADAMS=1,
                             BDF=2 } CZVODE_METH;

typedef enum czvode_iteration_method_e {NO_JACOBIAN=0,
                                       EXTERNAL_FULL_JACOBIAN = 1,
                                       INTERNAL_FULL_JACOBIAN = 2,
                                       INTERNAL_DIAG_JACOBIAN = 3,
                                       EXTERNAL_BAND_JACOBIAN = 4,
                                       INTERNAL_BAND_JACOBIAN = 5} CZVODE_MITER;


typedef enum czvode_ode_err_e {UNCHANGED = 2, SUCCESS = 0, ERROR = -1, UNKNOWN_ERROR = -20} CZVODE_ODE_RETVAL;



CZVODE_ODE_RETVAL zvode_basic(double* stack, ZVODE_COMPLEX* q, czvode_ode_func f_func,int neq, double t0, double tf, double dt, CZVODE_METH meth);

void zvode(void (*f)(const int *neq, const double *t, const ZVODE_COMPLEX *y, ZVODE_COMPLEX *ydot, ZVODE_COMPLEX* rpar, int* ipar),
           const int *neq,
           ZVODE_COMPLEX *y,
           double *t,
           const double *tout,
           const int *itol,
           const double *rtol,
           const double *atol,
           const int *itask,
           int *istate,
           const int *iopt,
           ZVODE_COMPLEX* zwork,
           const int *lzw,
           double *rwork,
           const int *lrw,
           int *iwork,
           const int *liw,
           void (*jac)(const int *neq, const double *t, const ZVODE_COMPLEX *y, const int *ml, const int *mu, ZVODE_COMPLEX *pd, const int *nrowpd, ZVODE_COMPLEX* rpar, int* ipar),
           const int *mf,
           ZVODE_COMPLEX* rpar,
           int* ipar);


#endif
