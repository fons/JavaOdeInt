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


