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
#ifndef H__GNI_IRK2__H
#define H__GNI_IRK2__H

typedef void (*gnicodes_ode_func)(int* n, double* x, double* q, double* f, double* rpar, int* ipar);

typedef void (*gnicodes_solfix)(int* nr, double* xold, double* x, double* p, double* q, int* n, int* irtrn, double* rpar, int* ipar);

typedef enum gni_irk2_method_e {IRK2_METH_2=2, IRK2_METH_4 = 4, IRK2_METH_6 = 6} GNI_IRK2_METH;
typedef enum gni_lmm2_method_e {LMM2_METH_201=201, LMM2_METH_401 = 401, LMM2_METH_801 = 801,LMM2_METH_802 = 802,LMM2_METH_803 = 803} GNI_LMM2_METH;

typedef enum gnicodes_iout_e { NEVER_CALLED=0, OUTPUT=1} GNICODES_IOUT;


void gni_irk2_basic(double* stack, double* p, double* q, gnicodes_ode_func f_func, int n, double xstart, double xfinal, double deltax, GNI_IRK2_METH meth);
void gni_lmm2_basic(double* stack, double* p, double* q, gnicodes_ode_func f_func, int n, double xstart, double xfinal, double deltax, GNI_LMM2_METH meth);

void gni_irk2(
    int* n,
    void(*f)(int* n, double* x, double* q, double* f, double* rpar, int* ipar),
    int* nstep,
    double* x,
    double* p,
    double* q,
    double* xend,
    int* meth,
    void (*s) (int* nr, double* xold, double* x, double* p, double* q, int* n, int* irtrn, double* rpar, int* ipar),
    int* iout,
    double* rpar,
    int* ipar
    );

void gni_lmm2(
    int* n,
    void(*f)(int* n, double* x, double* q, double* f, double* rpar, int* ipar),
    int* nstep,
    double* x,
    double* p,
    double* q,
    double* xend,
    int* meth,
    void (*s) (int* nr, double* xold, double* x, double* p, double* q, int* n, int* irtrn, double* rpar, int* ipar),
    int* iout,
    double* rpar,
    int* ipar
    );

#endif
