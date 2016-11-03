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

#ifndef H__DASPK__H
#define H__DASPK__H


typedef enum daspk_info_setting {DEFAULT = 0, CUSTOM = 1} DASPK_INFO_SETTING;
typedef enum daspk_constraint_check {INITIAL_CONDITION_ONLY = 1,
                                     ENFORCE_NON_NEGATIVE_Y = 2,
                                     BOTH_CONSTRAINTS       = 3 } DASPK_CONSTRAINT_CHECK;
typedef enum daspk_consistency_calc {CALC_ALG = 1, CALC_Y = 2} DASPK_CONSISTENCY_CALC;
/*
  elements of the INFO array
 */
typedef enum daspk_info {
    FIRST_CALL        =  0,
    TOLERANCE_TYPE    =  1,
    OVERSHOOT_TOUT    =  2,
    HAS_TCRIT         =  3,
    JACOBIAN          =  4,
    MATRIX            =  5,
    MAX_STEP_SIZE     =  6,
    INITIAL_STEP_SIZE =  7,
    MAX_ORDER         =  8,
    CONSTRAINTS       =  9,
    CONSISTENCY       = 10,
    DIRECT_OR_KRYLOV  = 11,
    KRYLOV_VALUES     = 12,
    CONTINUE_AFTER_IC_CALC = 13,
    KRYLOV_PRECOND    = 14,
    LOCAL_ERROR_CONTROL = 15,
    IC_HEURISTICS       = 16,
    EXTRA_PRINTING      = 17
} DASPK_INFO;
typedef enum daspk_idid_e {SUCCESS_NOT_TOUT    = 1,
                           SUCCESS_TSTOP       = 2,
                           SUCCESS_TOUT        = 3,
                           SUCCESS_INIT_CALC   = 4,
                           TOO_MANY_STEPS      = -1,
                           TOLERANCE_TOO_SMALL = -2,
                           LOCAL_ERROR_FAILED  = -3,
                           PRECOND_FAILURE     = -5,
                           REPEATED_FAILURE    = -6,
                           CONVERGENCE_FAILED  = -7,
                           SINGULAR_MATRIX     = -8,
                           REPEATED_FAILURE_ERROR    = -9,
                           REPEATED_FAILURE_IRES    = -10,
                           IRES_SIGNALLED_FAILURE   = -11,
                           INITIAL_VALUE_FAILED     = -12,
                           PSOL_FAILED              = -13,
                           KRYLOV_FAILED            = -14,
                           TASK_TERMINATED          = -33} DASPK_RETVAL;



void daspk(
    void(*res)(double*, double*, double*, double*, double* , int* , double*, int*),
    int* neq,
    double* t,
    double* y,
    double* yprime,
    double* tout,
    int* info,
    double* rtol,
    double* atol,
    int* idid,
    double* rwork,
    int* lrw,
    int* iwork,
    int* liw,
    double* rpar,
    int* ipar,
    void (*jac)(double*, double*, double*, double*, double*),
    void (*psol)(int*, double*, double*, double*, double*, double* , double*, double*, double*, int*, double* ));

#endif
