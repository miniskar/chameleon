/**
 *
 * @file starpu/codelet_zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset StarPU codelet
 *
 * @version 1.1.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @date 2021-03-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

struct cl_zlaset_args_s {
    cham_uplo_t uplo;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileA;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaset_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaset_args_s clargs;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    starpu_codelet_unpack_args( cl_arg, &clargs );
    TCORE_zlaset( clargs.uplo, clargs.m, clargs.n, clargs.alpha, clargs.beta, tileA );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlaset, cl_zlaset_cpu_func )

void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                         const CHAM_desc_t *A, int Am, int An )
{
    struct cl_zlaset_args_s clargs = {
        .uplo  = uplo,
        .m     = m,
        .n     = n,
        .alpha = alpha,
        .beta  = beta,
        .tileA = A->get_blktile( A, Am, An ),
    };
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    starpu_option_request_t *schedopt = (starpu_option_request_t *)(request->schedopt);
    int                      workerid;
    char                    *cl_name = "zlaset";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_W(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlaset_callback : NULL;

    /* Fix the worker id */
    workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlaset,
        /* Task codelet arguments */
        STARPU_VALUE, &clargs, sizeof(struct cl_zlaset_args_s),
        STARPU_W,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );
}
