#include <stdio.h>
#include <stdlib.h>

#include "codepack-internal.h"

static int lwm(CODEPACK_METHOD_FLAG mf, int neq, int nnz)
{
    switch (mf) {
        /*21,22,121,122*/
    case BDF_USER_FULL_JAC:
    case BDF_USER_JAC_NO_IA_JA:
        return 2*nnz + 2*neq + (nnz+9*neq)/2;
    case BDF_INTERNAL_FULL_JAC:
    case BDF_INTERNAL_JAC_IA_JA:
        return 2*nnz + 2*neq + (nnz+10*neq)/2;
    default :
        return 0;
    };

}
static int lrw(CODEPACK_METHOD_FLAG mf, int neq, int lwm)
{
    switch (mf) {
        /*10*/
    case ADAMS_BASIC :
        return 20 + 16 * neq;
        /*21,22,121,122*/
    case BDF_USER_FULL_JAC:
    case BDF_INTERNAL_FULL_JAC:
    case BDF_USER_JAC_NO_IA_JA:
    case BDF_INTERNAL_JAC_IA_JA:
        return 20 + 9 * neq + lwm;
    default :
        return 0;
    };
}

static int liw(CODEPACK_METHOD_FLAG mf, int neq, int nnz)
{

    switch (mf) {
        /*10;121;122  => moss =0; miter=0*/
    case ADAMS_BASIC :
    case BDF_USER_JAC_NO_IA_JA:
    case BDF_INTERNAL_JAC_IA_JA:
        return 30;
        /*21,22 MOSS=0; MTER=1,2*/
    case BDF_USER_FULL_JAC:
    case BDF_INTERNAL_FULL_JAC:
        return 31 + neq + nnz;

    default :
        return 0;
    
    }
}

static int supported_methods_bool(CODEPACK_METHOD_FLAG mf)
{
    switch (mf) {
    case ADAMS_BASIC:
    case BDF_USER_FULL_JAC:
    case BDF_INTERNAL_FULL_JAC:
    case  BDF_USER_JAC_NO_IA_JA:
    case  BDF_INTERNAL_JAC_IA_JA:
        return 1;
    default :
        return 0;
    }
}

static lsod_params* create_basic_lsodes_params(int neq, codepack_ode_func f_func, CODEPACK_METHOD_FLAG mf)
{
    if (! supported_methods_bool(mf)) {
        return NULL;
    }
    const double rtol = 0.0;
    const double atol = 1e-12;
    lsod_params* dsp = calloc(1, sizeof(lsod_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->ode_func = DLSODES;
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->itask = NORMAL;
    dsp->istate.istate_in = FIRST_CALL;
    dsp->iopt     = NO_OPTIONAL_INPUTS;
    dsp->lrw      = lrw(mf, neq, lwm(mf, neq, neq*neq));
    dsp->rwork    = (double*) calloc(dsp->lrw, sizeof(double));
    dsp->liw      = liw(mf,neq, neq*neq);
    dsp->iwork    = (int *)  calloc(dsp->liw, sizeof(int));
    dsp->jac.jac2 = NULL;
    dsp->mf       = mf;

    return dsp;
}

static void free_params(lsod_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}

static CODEPACK_ISTATE_OUT lsodes(lsod_params *dlsodap, double tnext, double *t, double *y)
{

    //fprintf(stderr, "calling dlsoda_\n");
    dlsodes_(dlsodap->f_func,
            &dlsodap->neq,
            y,
            t,
            &tnext,
            &dlsodap->itol,
            &dlsodap->rtol.rtol_val,
            &dlsodap->atol.atol_val,
            (const int*) &dlsodap->itask,
            &dlsodap->istate.istate,
            (const int*) &dlsodap->iopt,
            dlsodap->rwork,
            &dlsodap->lrw,
            dlsodap->iwork,
            &dlsodap->liw,
            dlsodap->jac.jac2,
            (const int*)&dlsodap->mf);

    return dlsodap->istate.istate_out;
}

CODEPACK_ODE_RETVAL lsodes_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt, CODEPACK_METHOD_FLAG mf)
{
    double t = 0.0;
    double tnext = 0;
    int max_retries = 5;
    int retry = 0;
    CODEPACK_ODE_RETVAL ode_ret = SUCCESS;
    int index = 0;
    lsod_params* dlsop = create_basic_lsodes_params(neq, f_func, mf);

    if (dlsop == NULL) {
        return ERROR;
    }

    stack = write_to_stack(stack, neq, &index, t, q);    
    t = t0;

    while(t < tf){
        retry = 0;
        tnext = t + dt;
        CODEPACK_ISTATE_OUT return_code = SUCCESS_DONE;
        do {
            return_code = lsodes(dlsop, tnext, &t, q);
            if (return_code == MAX_STEPS_EXCEEDED) {
                retry++;
                if (retry >= max_retries) {
                    break;
                }
                dlsop->iopt      = OPTIONAL_INPUTS;
                dlsop->iwork[5] += 2000;
                dlsop->istate.istate_in = NEXT_CALL_WITH_CHANGES;
                fprintf(stderr, "increased max steps to %d for retry %d \n", dlsop->iwork[5], retry);
                t = tnext - dt;
            }
            else {
                retry = 0;
                break;
            }
        } while (retry > 0);
        ode_ret = istate(return_code);
        if (ode_ret < 0) {
            break;
        }
        stack = write_to_stack(stack, neq, &index, t, q);
    }

    free_params(dlsop);
    return ode_ret;   
}

void dlsodes(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
             const int *mf)
{
    dlsodes_(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf);
}
