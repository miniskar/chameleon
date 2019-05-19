/**
 *
 * @file starpu/runtime_context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU context routines
 *
 * @version 0.9.2
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2015-05-22
 *
 */
#include <stdlib.h>
#include "chameleon_starpu.h"

#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 3))
/* Defined by StarPU as external function */
#else
#if ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 2))
int _starpu_is_initialized(void);
#define starpu_is_initialized() _starpu_is_initialized()
#else
#define starpu_is_initialized() 0
#endif
#endif

/**
 *  Create new context
 */
void RUNTIME_context_create( CHAM_context_t *chamctxt )
{
    starpu_conf_t *conf;

    chamctxt->scheduler = RUNTIME_SCHED_STARPU;

    if (! starpu_is_initialized() ) {
        chamctxt->schedopt = (void*) malloc (sizeof(struct starpu_conf));
        conf = chamctxt->schedopt;

        starpu_conf_init( conf );
    }
    else {
        chamctxt->schedopt = NULL;
    }

    return;
}

/**
 *  Clean the context
 */
void RUNTIME_context_destroy( CHAM_context_t *chamctxt )
{
    /* StarPU was already initialized by an external library */
    if (chamctxt->schedopt) {
        free(chamctxt->schedopt);
    }
    return;
}

/**
 *
 */
void RUNTIME_enable( void *runtime_ctxt, int lever )
{
    switch (lever)
    {
    case CHAMELEON_DAG:
        fprintf(stderr, "StarPU is providing DAG generation through tracing support (CHAMELEON_PROFILING_MODE)\n");
        break;
    case CHAMELEON_PROFILING_MODE:
        starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
        break;
    case CHAMELEON_BOUND:
        starpu_bound_start(0, 0);
        break;
    default:
        return;
    }

    (void)runtime_ctxt;
    return;
}

/**
 *
 */
void RUNTIME_disable( void *runtime_ctxt, int lever )
{
    switch (lever)
    {
    case CHAMELEON_DAG:
        fprintf(stderr, "StarPU is providing DAG generation through tracing support (CHAMELEON_PROFILING_MODE)\n");
        break;
    case CHAMELEON_PROFILING_MODE:
        starpu_profiling_status_set(STARPU_PROFILING_DISABLE);
        break;
    case CHAMELEON_BOUND:
        starpu_bound_stop();
        break;
    default:
        return;
    }

    (void)runtime_ctxt;
    return;
}
