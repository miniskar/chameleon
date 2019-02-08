/**
 *
 * @file starpu/codelet_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zgemm(const RUNTIME_option_t *options,
                      cham_trans_t transA, cham_trans_t transB,
                      int m, int n, int k, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                                                   const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                      CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgemm;
    void (*callback)(void*) = options->profiling ? cl_zgemm_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &transA,            sizeof(int),
        STARPU_VALUE,    &transB,            sizeof(int),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &alpha,             sizeof(CHAMELEON_Complex64_t),
        STARPU_R,         RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_R,         RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_VALUE,    &ldb,               sizeof(int),
        STARPU_VALUE,    &beta,              sizeof(CHAMELEON_Complex64_t),
        STARPU_RW,        RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_VALUE,    &ldc,               sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgemm",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgemm_cpu_func(void *descr[], void *cl_arg)
{
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *C;
    int ldc;

    A = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);
    CORE_zgemm(transA, transB,
        m, n, k,
        alpha, A, lda,
        B, ldb,
        beta, C, ldc);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zgemm_cuda_func(void *descr[], void *cl_arg)
{
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    cuDoubleComplex alpha;
    const cuDoubleComplex *A;
    int lda;
    const cuDoubleComplex *B;
    int ldb;
    cuDoubleComplex beta;
    cuDoubleComplex *C;
    int ldc;

    A = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);

    RUNTIME_getStream( stream );

    CUDA_zgemm(
        transA, transB,
        m, n, k,
        &alpha, A, lda,
                B, ldb,
        &beta,  C, ldc,
        stream);

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
CODELETS(zgemm, 3, cl_zgemm_cpu_func, cl_zgemm_cuda_func, STARPU_CUDA_ASYNC)
