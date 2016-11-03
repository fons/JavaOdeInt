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

#ifndef H__CMEBDFI__H
#define H__CMEBDFI__H

typedef enum mebdfi_method_e {
    BDF_USER_FULL_JAC         = 21,
    BDF_INTERNAL_FULL_JAC     = 22,
    BDF_USER_BAND_JAC         = 23,
    BDF_INTERNAL_DIAG_JAC     = 24
} MEBDFI__METHOD_FLAG;

typedef enum mebdfi_itol_e { ALL_SCALAR = 2, ATOL_ARRAY=3, RTOL_ARRAY=4, ALL_ARRAY=5} MEBDFI_ITOLERANCE;
typedef enum mebdfi_idid_in_e {CONTINUE   = 0,
                               FIRST_CALL = 1,
                               CONTINUE_WITH_CHANGES = -1,
                               STOP_AT_FIRST_MESH = 2,
                               ONESTEP_ONLY = 3 } MEBDFI_IDID_IN;

typedef enum mebdfi_idid_out_e { SUCCESS_DONE                = 0,
                                 FAILED_ERROR_TEST           = -1,
                                 REPEATED_FAILURES           = -2,
                                 CONVERGENCE_FAIlURES        = -3 ,
                                 ILLEGAL_INPUT               = -4,
                                 PARAM_CHANGE_NOT_APPLIED    = -5,
                                 INTEGRATION_STEPS_EXCEEDED  = -6,
                                 STEPSIZE_TOO_SMALL          = -7,
                                 INSUFFICIENT_WORK_SPACE     = -11,
                                 INSUFFICIENT_IWORK_SPACE    = -12} MEBDFI_IDID_OUT;

void mebdfi(
    int* n,
    double* t0,
    double* h0,
    double* y0,
    double* yprime,
    double* tout,
    double* tend,
    int* mf,
    int* idid,
    int* lout,
    int* lwork,
    double* work,
    int* liwork,
    int* iwork,
    int* mbnd,
    int* maxder,
    int* itol,    
    double* rtol,
    double* atol,
    double* rpar,
    int* ipar,
            void (*pderv)(double* , /*t*/
                          double* , /*y*/
                          double* , /*pd*/
                          int*,    /*n*/
                          double* , /*yprime*/
                          int*,     /*mbnd*/
                          double* , /*con*/
                          int*    , /*ipar*/
                          double* , /*rpar*/
                          int* /*ier*/
                          ),

    void (*resid)(int*, double*, double* , double*, double*, int*, double* , int* ),
    int *ierr
    );

#endif
