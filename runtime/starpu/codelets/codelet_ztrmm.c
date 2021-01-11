/**
 *
 * @file starpu/codelet_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm StarPU codelet
 *
 * @version 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztrmm_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t transA;
    cham_diag_t diag;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &M, &N, &alpha);
    TCORE_ztrmm(side, uplo,
        transA, diag,
        M, N,
        alpha, tileA,
        tileB);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_ztrmm_cuda_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t transA;
    cham_diag_t diag;
    int M;
    int N;
    cuDoubleComplex alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &M, &N, &alpha);

    RUNTIME_getStream(stream);

    CUDA_ztrmm(
        side, uplo, transA, diag, M, N, &alpha,
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif /* CHAMELEON_USE_CUDA */
#endif /* !defined(CHAMELEON_SIMULATION) */


/*
 * Codelet definition
 */
CODELETS(ztrmm, cl_ztrmm_cpu_func, cl_ztrmm_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_ztrmm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                      const CHAM_desc_t *B, int Bm, int Bn)
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlaset( options, ChamUpperLower, m, n,
                                   alpha, alpha, B, Bm, Bn );
    }

    (void)nb;
    struct starpu_codelet *codelet = &cl_ztrmm;
    void (*callback)(void*) = options->profiling ? cl_ztrmm_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,      &side,                sizeof(int),
        STARPU_VALUE,      &uplo,                sizeof(int),
        STARPU_VALUE,    &transA,                sizeof(int),
        STARPU_VALUE,      &diag,                sizeof(int),
        STARPU_VALUE,         &m,                        sizeof(int),
        STARPU_VALUE,         &n,                        sizeof(int),
        STARPU_VALUE,     &alpha,         sizeof(CHAMELEON_Complex64_t),
        STARPU_R,                  RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_RW,                 RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztrmm",
#endif
        0);
}
