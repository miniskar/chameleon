/**
 *
 * @file starpu/codelet_zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @date 2014-11-16
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
    CHAMELEON_Complex64_t *A;
    int ldA;
    CHAMELEON_Complex64_t *B;
    int ldB;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *C;
    int ldC;

    A = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldB = STARPU_MATRIX_GET_LD( descr[1] );
    ldC = STARPU_MATRIX_GET_LD( descr[2] );

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &M, &N, &alpha, &beta);
    CORE_zhemm(side, uplo,
        M, N,
        alpha, A, ldA,
        B, ldB,
        beta, C, ldC);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zhemm_cuda_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_uplo_t uplo;
    int M;
    int N;
    cuDoubleComplex alpha;
    const cuDoubleComplex *A;
    int ldA;
    const cuDoubleComplex *B;
    int ldB;
    cuDoubleComplex beta;
    cuDoubleComplex *C;
    int ldC;

    A = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);

    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldB = STARPU_MATRIX_GET_LD( descr[1] );
    ldC = STARPU_MATRIX_GET_LD( descr[2] );

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &M, &N, &alpha, &beta);

    RUNTIME_getStream(stream);

    CUDA_zhemm(
        side, uplo,
        M, N,
        &alpha, A, ldA,
        B, ldB,
        &beta, C, ldC,
        stream);

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
CODELETS(zhemm, 3, cl_zhemm_cpu_func, cl_zhemm_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zhemm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int ldA,
                      const CHAM_desc_t *B, int Bm, int Bn, int ldB,
                      CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldC)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zhemm;
    void (*callback)(void*) = options->profiling ? cl_zhemm_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &side,                sizeof(int),
        STARPU_VALUE,    &uplo,                sizeof(int),
        STARPU_VALUE,       &m,                        sizeof(int),
        STARPU_VALUE,       &n,                        sizeof(int),
        STARPU_VALUE,   &alpha,         sizeof(CHAMELEON_Complex64_t),
        STARPU_R,               RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,               RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_VALUE,    &beta,         sizeof(CHAMELEON_Complex64_t),
        STARPU_RW,               RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zhemm",
#endif
        0);
    (void)ldC;
    (void)ldB;
    (void)ldA;
}
