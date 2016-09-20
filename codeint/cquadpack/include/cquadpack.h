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
#ifndef H__CQUADPACK__H
#define H__CQUADPACK__H

typedef double (*cquadpack_ode_func)(double*); 

typedef enum quadpack_errno_e {SUCCESS=0,
                               MAX_SUBDIV=1,
                               ROUNDOFF_ER=2,
                               BAD_INTEGRAND=3,
                               NO_CONVERGANCE=4,
                               SLOW_CONVERGANCE=5,
                               INVALID_INPUT=6 ,
                               BAD_INTEGRAND_BEHAVIOUR=7} QUADPACK_ERRNO;

typedef enum quadpack_infinity_e {PLUS_INFINITY=1, MINUS_INFINITY=-1, MINUS_TO_PLUS_INFINITY=2} QUADPACK_INFINITY; 

typedef enum quadpack_trig_weight_function_e {COS=1, SIN=2} QUADPACK_TRIG_WEIGHT_FUNCTION;
typedef enum quadpack_log_weight_function_e  {LOGW_1=1,LOGW_2=2,LOGW_3=3,LOGW_4=4} QUADPACK_LOG_WEIGHT_FUNCTION; 

static char* LOG_WEIGHT[20] = {"1  (x-a)**alfa*(b-x)**beta",
                               "2  (x-a)**alfa*(b-x)**beta*log(x-a)",
                               "3  (x-a)**alfa*(b-x)**beta*log(b-x)",
                               "4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)"
};

char* log_weight_to_string(QUADPACK_LOG_WEIGHT_FUNCTION w);

QUADPACK_ERRNO qagse_basic(cquadpack_ode_func func, double start, double end, double* result);
QUADPACK_ERRNO qagie_basic(cquadpack_ode_func func, double bound, QUADPACK_INFINITY infinity, double* result);
QUADPACK_ERRNO qawoe_basic(cquadpack_ode_func func, double a, double b, double omega, QUADPACK_TRIG_WEIGHT_FUNCTION integr, double* result);
QUADPACK_ERRNO qawce_basic(cquadpack_ode_func func, double start, double end, double weight, double* result);
QUADPACK_ERRNO qawfe_basic(cquadpack_ode_func func, double a, double omega, QUADPACK_TRIG_WEIGHT_FUNCTION integr, double* result);
QUADPACK_ERRNO qawse_basic(cquadpack_ode_func func, double a, double b, double alfa, double beta, QUADPACK_LOG_WEIGHT_FUNCTION integr, double* result);
QUADPACK_ERRNO qagpe_basic(cquadpack_ode_func func, double start, double end, int npts2, double* points, double* result);

void dqagse(double (*f)(double *y),
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

void dqagie(double (*f)(double *y),
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

void dqawoe(double (*f)(double *y),
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

void dqawce(double (*f)(double *y),
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

void dqawfe(double (*f)(double *y),
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

void dqawse(double (*f)(double *y),
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


void dqagpe(double (*f)(double *y),
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



#endif
