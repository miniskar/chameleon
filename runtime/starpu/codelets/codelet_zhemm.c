/**
 *
 * @file starpu/codelet_zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm StarPU codelet
 *
 * @version 1.2.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @date 2022-02-22
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zhemm_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &M, &N, &alpha, &beta);
    TCORE_zhemm(side, uplo,
        M, N,
        alpha, tileA,
        tileB,
        beta, tileC);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zhemm_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    cham_side_t side;
    cham_uplo_t uplo;
    int M;
    int N;
    cuDoubleComplex alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    cuDoubleComplex beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &M, &N, &alpha, &beta);

    CUDA_zhemm(
        side, uplo,
        M, N,
        &alpha, tileA->mat, tileA->ld,
                tileB->mat, tileB->ld,
        &beta,  tileC->mat, tileC->ld,
        handle );
}
#endif /* CHAMELEON_USE_CUDA */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zhemm, cl_zhemm_cpu_func, cl_zhemm_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zhemm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                      const CHAM_desc_t *B, int Bm, int Bn,
                      CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn)
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                                    beta, C, Cm, Cn );
    }

    (void)nb;
    struct starpu_codelet *codelet = &cl_zhemm;
    void (*callback)(void*) = options->profiling ? cl_zhemm_callback : NULL;
    int accessC = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &side,                sizeof(int),
        STARPU_VALUE,    &uplo,                sizeof(int),
        STARPU_VALUE,       &m,                        sizeof(int),
        STARPU_VALUE,       &n,                        sizeof(int),
        STARPU_VALUE,   &alpha,         sizeof(CHAMELEON_Complex64_t),
        STARPU_R,               RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,               RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_VALUE,    &beta,         sizeof(CHAMELEON_Complex64_t),
        accessC,                RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zhemm",
#endif
        0);
}
