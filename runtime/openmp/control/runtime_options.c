/**
 *
 * @file openmp/runtime_options.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP options routines
 *
 * @version 1.2.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Philippe Virouleau
 * @date 2022-02-22
 *
 */
#include "chameleon_openmp.h"

void RUNTIME_options_init( RUNTIME_option_t *options, CHAM_context_t *chamctxt,
                           RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    options->sequence  = sequence;
    options->request   = request;
    options->profiling = CHAMELEON_STATISTICS == CHAMELEON_TRUE;
    options->parallel  = CHAMELEON_PARALLEL == CHAMELEON_TRUE;
    options->priority  = RUNTIME_PRIORITY_MIN;
    options->workerid  = -1;
    options->forcesub  = 0;
    options->ws_wsize  = 0;
    options->ws_hsize  = 0;
    options->ws_worker = NULL;
    options->ws_host   = NULL;
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
    if ( worker_size > 0 ) {
        /*
         * NOTE: we set the size, but instead of doing a malloc shared by multiple workers,
         * we just create a VLA in the relevant codelets, within the task's body.
         * This way we ensure the "scratch" is thread local and not shared by multiple threads.
         */
        options->ws_wsize = worker_size;
    }
    (void)host_size;
    return CHAMELEON_SUCCESS;
}

int RUNTIME_options_ws_free( RUNTIME_option_t *options )
{
    if (options->ws_wsize) {
        options->ws_wsize = 0;
    }
    return CHAMELEON_SUCCESS;
}
