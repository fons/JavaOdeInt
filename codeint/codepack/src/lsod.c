#include <stdio.h>
#include <stdlib.h>
#include "codepack-internal.h"

CODEPACK_ODE_RETVAL istate(CODEPACK_ISTATE_OUT istat)
{
    switch (istat) {
    case NOTHING_DONE:
        return 2;
    case SUCCESS_DONE:
        return 0;
    case MAX_STEPS_EXCEEDED:
    case TO_MUCH_ACCURACY:
    case ILLEGAL_INPUT:
    case ERROR_TEST_FAILURES:
    case CONVERGENCE_FAILURES:
    case ZERO_ERR_TOLERANCE:
    case TOO_SMALL_WORK_ARRAY:
        return -1;
    default :
        return -20;
    }

}

