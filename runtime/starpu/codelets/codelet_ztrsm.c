/**
 *
 * @file starpu/codelet_ztrsm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm StarPU codelet
 *
 * @version 1.1.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
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

struct cl_ztrsm_args_s {
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t transA;
    cham_diag_t diag;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_ztrsm_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_ztrsm_args_s clargs;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &clargs );
    TCORE_ztrsm( clargs.side, clargs.uplo, clargs.transA, clargs.diag,
                 clargs.m, clargs.n, clargs.alpha, tileA, tileB );
}

#ifdef CHAMELEON_USE_CUDA
static void
cl_ztrsm_cuda_func(void *descr[], void *cl_arg)
{
    struct cl_ztrsm_args_s clargs;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &clargs );

    RUNTIME_getStream(stream);

    CUDA_ztrsm(
        clargs.side, clargs.uplo, clargs.transA, clargs.diag,
        clargs.m, clargs.n,
        (cuDoubleComplex*)&(clargs.alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( ztrsm, cl_ztrsm_cpu_func, cl_ztrsm_cuda_func, STARPU_CUDA_ASYNC )

void INSERT_TASK_ztrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn )
{
    struct cl_ztrsm_args_s clargs = {
        .side   = side,
        .uplo   = uplo,
        .transA = transA,
        .diag   = diag,
        .m      = m,
        .n      = n,
        .alpha  = alpha,
        .tileA  = A->get_blktile( A, Am, An ),
        .tileB  = B->get_blktile( B, Bm, Bn ),
    };
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    starpu_option_request_t *schedopt = (starpu_option_request_t *)(request->schedopt);
    int                      workerid;
    char                    *cl_name = "ztrsm";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Callback fro profiling information */
    callback = options->profiling ? cl_ztrsm_callback : NULL;

    /* Fix the worker id */
    workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ztrsm,
        /* Task codelet arguments */
        STARPU_VALUE, &clargs, sizeof(struct cl_ztrsm_args_s),
        STARPU_R,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_RW,     RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );

    (void)nb;
}
