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
