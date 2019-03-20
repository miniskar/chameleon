/**
 *
 * @file quark/runtime_options.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark options routines
 *
 * @version 0.9.2
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @date 2015-05-22
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "chameleon_quark.h"

void RUNTIME_options_init( RUNTIME_option_t *options, CHAM_context_t *chamctxt,
                           RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    /* Create the task flag */
    quark_option_t *qopt = (quark_option_t*) malloc(sizeof(quark_option_t));

    /* Initialize task_flags */
    memset( qopt, 0, sizeof(Quark_Task_Flags) );
    qopt->flags.task_lock_to_thread = -1;
    qopt->flags.task_thread_count   =  1;
    qopt->flags.thread_set_to_manual_scheduling = -1;

    /* Initialize options */
    options->sequence   = sequence;
    options->request    = request;
    options->profiling  = CHAMELEON_PROFILING == CHAMELEON_TRUE;
    options->parallel   = CHAMELEON_PARALLEL == CHAMELEON_TRUE;
    options->priority   = RUNTIME_PRIORITY_MIN;

    options->ws_wsize   = 0;
    options->ws_hsize   = 0;
    options->ws_worker  = NULL;
    options->ws_host    = NULL;

    /* quark in options */
    qopt->quark = (Quark*)(chamctxt->schedopt);
    options->schedopt = qopt;

    QUARK_Task_Flag_Set((Quark_Task_Flags*)qopt, TASK_SEQUENCE, (intptr_t)(sequence->schedopt));

    return;
}

void RUNTIME_options_finalize( RUNTIME_option_t *options, CHAM_context_t *chamctxt )
{
    /* we can free the task_flags without waiting for quark
       because they should have been copied for every task */
    (void)chamctxt;
    free( options->schedopt );
    return;
}

int RUNTIME_options_ws_alloc( RUNTIME_option_t *options, size_t worker_size, size_t host_size )
{
    options->ws_wsize = worker_size;
    options->ws_hsize = host_size;
    return CHAMELEON_SUCCESS;
}

int RUNTIME_options_ws_free( RUNTIME_option_t *options )
{
    options->ws_wsize = 0;
    options->ws_hsize = 0;
    return CHAMELEON_SUCCESS;
}
