#ifndef H__CODEPACK_INTERNAL__H
#define H__CODEPACK_INTERNAL__H

#include "../include/codepack.h"

/*
  declaring the odepack fortran interfaces
 */


extern void dlsode_(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
                    void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd),
                    const int *mf);

extern void dlsodes_(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
                     void (*jac)(const int *neq, const double *t, const double *y, const int *j, const double *ian, double *jan, double *pdj),
                    const int *mf);

extern void dlsoda_(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
                    void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd),
                    const int *jt);


extern void dlsodar_(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
                     void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd), const int *jt,\
                     void (*g)(const int *neq, const double *t, const double *y, const int *ng, double *gout), const int *ng,
                     int *jroot);



#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#define MIN(a, b) ((a) <= (b) ? (a) : (b))

typedef enum odepack_function_e {DLSODE, DLSODES, DLSODA, DLSODAR} ODEPACK_FUNCTION;

/*
 * this captures all the input params 
 */
typedef struct lsod_params_s {
    ODEPACK_FUNCTION ode_func;
    codepack_ode_func f_func;
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

    CODEPACK_ITASK itask;
    
    union {
        CODEPACK_ISTATE_IN  istate_in;
        CODEPACK_ISTATE_OUT istate_out;
        int istate;
    } istate;
    
    CODEPACK_OPTIONAL_INPUT_FLAG iopt;
    double* rwork;
    int* iwork;
    int lrw;
    int liw;
    union {
        codepack_jac_func_1 jac1;
        codepack_jac_func_2 jac2;
    } jac;
    CODEPACK_JAC_TYPE jt;
    CODEPACK_METHOD_FLAG mf;
    codepack_g_func g;
    int ng;
    int* jroot;
} lsod_params;

int istate(CODEPACK_ISTATE_OUT istat);
double* write_to_stack(double* stack, int neq, int* index, double t_new, double* q);


lsod_params* create_basic_lsode_params(int neq, codepack_ode_func f_func, CODEPACK_METHOD_FLAG mf);
CODEPACK_ISTATE_OUT lsode(double t, double *t0, double *q, lsod_params *sodap);


lsod_params* create_basic_lsodes_params(int neq, codepack_ode_func f_func, CODEPACK_METHOD_FLAG mf);
CODEPACK_ISTATE_OUT lsodes(double t, double *t0, double *q, lsod_params *sodap);

lsod_params* create_basic_lsoda_params(int neq, codepack_ode_func f_func);
CODEPACK_ISTATE_OUT lsoda(double t, double *t0, double *q, lsod_params *sodap);

lsod_params* create_basic_lsodar_params(int neq, codepack_ode_func f_func);
CODEPACK_ISTATE_OUT lsodar(double t, double *t0, double *q, lsod_params *sodap);


#endif


