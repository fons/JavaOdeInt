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

#ifndef H__CBIMD__H
#define H__CBIMD__H
typedef enum bimd_itol_e { ALL_SCALAR = 0, ALL_ARRAY=1} BIMD_ITOLERANCE;
typedef enum bimd_jac_type_e {NUMERICAL_JACOBIAN = 0, USER_PROVIDED = 1} BIMD_JAC_TYPE;
typedef enum bimd_mass_matrix_e {IDENTITY_MATRIX=0, MASS_USER_PROVIDED=1} BIMD_MASS_MATRIX;
typedef enum bimd_iout_e { NEVER_CALLED=0, OUTPUT=1} BIMD_IOUT;

typedef enum bimd_idid_e {SUCCESS=1,
                          SUCCESS_INTR=2,
                          INPUT_INCONSISTENT=-1,
                          NMAX_TOO_SMALL=-2,
                          STEP_TOO_SMALL=-3,
                          SINGULAR_MATRIX=-4,
                          TOO_MANY_NEWTON_FAILURES,
                          FCN_JAC_ERROR_RETURNED} BIMD_RETVAL;

void bimd(
    int* m,
    void (*fcn)(int* , double* , double* , double* , int* , double* , int* ),
    double* t0,
    double* tend,
    double* y0,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int* , double* , double* , double** , int* , int*, double* , int* ),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int* , double* ,int* ,int*, double* , int* ),
    int* imas,
    int* mlmas,
    int* mumas,
    void (*solout)(int* , int*, int*, double* , double* , double* , double* , double**, double*, int* , int* ),
    int* iout,
    double* work,
    int* lwork,
    int* iwork,
    int* liwork,
    double* rpar,
    int* ipar,
    int *idid
    );

#endif
