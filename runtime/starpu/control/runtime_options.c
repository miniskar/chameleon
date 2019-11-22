/**
 *
 * @file starpu/runtime_options.c
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
 * @date 2015-05-22
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "chameleon_starpu.h"

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
    int ret = 0;
    if ( worker_size > 0 ) {
        CHAM_tile_t tile = {
            .format = CHAMELEON_TILE_FULLRANK,
            .m      = worker_size,
            .n      = 1,
            .ld     = worker_size,
            .mat    = NULL,
        };
        options->ws_wsize = worker_size;
        starpu_cham_tile_register( (starpu_data_handle_t*)(&(options->ws_worker)),
                                   -1, &tile, sizeof(char) );
    }
    if ( host_size > 0 ) {
        options->ws_hsize = host_size;
        ret = RUNTIME_starpu_ws_alloc((CHAMELEON_starpu_ws_t**)&(options->ws_host),
                                      host_size, CHAMELEON_CUDA, CHAMELEON_HOST_MEM);
    }
    return ret;
}

int RUNTIME_options_ws_free( RUNTIME_option_t *options )
{
    int ret = 0;
    if ( options->ws_worker != NULL ) {
        starpu_data_unregister_submit((starpu_data_handle_t)(options->ws_worker));
        options->ws_worker = NULL;
    }
    if ( options->ws_host != NULL ) {
        starpu_task_wait_for_all();
        ret = RUNTIME_starpu_ws_free( (CHAMELEON_starpu_ws_t*)(options->ws_host) );
        options->ws_host = NULL;
    }
    return ret;
}
