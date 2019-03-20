/**
 *
 * @file openmp/runtime_options.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU options routines
 *
 * @version 0.9.2
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2018-06-15
 *
 */
#include <stdlib.h>
#include "chameleon_openmp.h"

void RUNTIME_options_init( RUNTIME_option_t *option, CHAM_context_t *chamctxt,
                           RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    option->sequence   = sequence;
    option->request    = request;
    option->profiling  = CHAMELEON_PROFILING == CHAMELEON_TRUE;
    option->parallel   = CHAMELEON_PARALLEL == CHAMELEON_TRUE;
    option->priority   = RUNTIME_PRIORITY_MIN;
    option->ws_wsize   = 0;
    option->ws_hsize   = 0;
    option->ws_worker  = NULL;
    option->ws_host    = NULL;
    return;
}

void RUNTIME_options_finalize( RUNTIME_option_t *option, CHAM_context_t *chamctxt )
{
    (void)option;
    (void)chamctxt;
    return;
}

int RUNTIME_options_ws_alloc( RUNTIME_option_t *options, size_t worker_size, size_t host_size )
{
    if (worker_size > 0) {
        /*
         * NOTE: we set the size, but instead of doing a malloc shared by multiple workers,
         * we just create a VLA in the relevant codelets, within the task's body.
         * This way we ensure the "scratch" is thread local and not shared by multiple threads.
         */
        options->ws_wsize = worker_size;
    }
    return CHAMELEON_SUCCESS;
}

int RUNTIME_options_ws_free( RUNTIME_option_t *options )
{
    if (options->ws_wsize) {
        options->ws_wsize = 0;
    }
    return CHAMELEON_SUCCESS;
}
