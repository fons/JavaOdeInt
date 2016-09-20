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
#ifndef H__CQUADPACK_INTERNAL__H
#define H__CQUADPACK_INTERNAL__H

extern void dqagse_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqagie_(double (*f)(double *y),
                    double* a,
                    int* inf,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqawoe_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* omega,
                    int* integr,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    int* icall,
                    int* maxpl,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    int* last,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* nnlog,
                    int* momcom,
                    double* chebmo
    );

extern void dqawce_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* c,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqawfe_(double (*f)(double *y),
                    double* a,
                    double* omega,
                    int* integr,
                    double* epsabs,
                    int* limlst,
                    int* limit,
                    int* maxpl,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* rslst,
                    double* erlst,
                    int*    ierlst,
                    int* lst,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* nnlog,
                    double* chebmo
    );

extern void dqawse_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* alfa,
                    double* beta,
                    int* integr,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqagpe_(double (*f)(double *y),
                    double* a,
                    double* b,
                    int* npts2,
                    double* points,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    double* pts,
                    int* level,
                    int*  ndin,
                    int* iord,
                    int* last);

typedef struct quadpack_input_s {
    double epsabs;
    double epsrel;
    int    limlst;
    int    limit;
    int    icall;
    int    maxpl;
} QUADPACK_INPUT;

typedef struct quadpack_output_s {
    double  abserr;
    int     neval;
    int     ier;
    double* rslst;
    double* erlst;
    int*    ierlst;
    int     lst;
    double* alist;
    double* blist;
    double* rlist;
    double* elist;
    int*    iord;
    int     last;
    int*    nnlog;
    int     momcom;
    double* chebmo;
    double* pts;
    int*    level;
    int*    ndin;
} QUADPACK_OUTPUT;

typedef double QUADPACK_RESULT;


void default_quadpack_input(QUADPACK_INPUT* qpi);
void default_quadpack_output(QUADPACK_INPUT qpi, QUADPACK_OUTPUT* qpo);
void points_quadpack_output(int npts2, QUADPACK_OUTPUT* qpo);
void free_quadpack_output(QUADPACK_OUTPUT* qpo);
#endif
