#ifndef H__GNI_IRK2__INTERNAL__H
#define H__GNI_IRK2__INTERNAL__H
#include "../../codeintdeps/include/stack.h"
#include "../include/cgnicodes.h"

typedef struct gnicodes_params_s
{
    gnicodes_ode_func f_func;
    gnicodes_solfix   solfix;
    int neq;
    int nstep;
    GNICODES_IOUT iout;
    union {
        GNI_IRK2_METH irk2_meth;
        GNI_LMM2_METH lmm2_meth;
    } meth;
    int lr;
    double *rpar;
    int li;
    int *ipar;
} gnicodes_params;

extern void gni_irk2_(
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

extern void gni_lmm2_(
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
