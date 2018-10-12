/**
 *
 * @file parsec/runtime_context.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC context routines
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2017-01-12
 *
 */
#include <stdlib.h>
#include "chameleon_parsec.h"

/**
 *  Create new context
 */
void RUNTIME_context_create( CHAM_context_t *chamctxt )
{
    /* In case of PaRSEC, this is done in init */
    chamctxt->scheduler = RUNTIME_SCHED_PARSEC;
    return;
}

/**
 *  Clean the context
 */
void RUNTIME_context_destroy( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

/**
 *
 */
void RUNTIME_enable(int lever)
{
    switch (lever)
    {
    case CHAMELEON_PROFILING_MODE:
        break;
    default:
        return;
    }
    return;
}

/**
 *
 */
void RUNTIME_disable(int lever)
{
    switch (lever)
    {
    case CHAMELEON_PROFILING_MODE:
        break;
    default:
        return;
    }
    return;
}
