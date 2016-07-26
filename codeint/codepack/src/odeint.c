
#include <stdio.h>
#include <stdlib.h>


#include "codepack-internal.h"



CODEPACK_ODE_RETVAL odeint(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt)
{
    double t;
    CODEPACK_ODE_RETVAL ode_ret;
    int index = 0;
    lsod_params* dlsop = create_basic_lsoda_params(neq, f_func);
    t = t0;
    while(t < tf){    
        CODEPACK_ISTATE_OUT ret = lsoda(t + dt, &t, q, dlsop);
        ode_ret = istate(ret);
        if (ode_ret < 0) {
            return ode_ret;
        }
        stack = write_to_stack(stack, neq, &index, (t+dt), q);
    }
    
    return ode_ret;
}


