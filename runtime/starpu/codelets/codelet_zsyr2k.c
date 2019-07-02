/**
 *
 * @file starpu/codelet_zsyr2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyr2k StarPU codelet
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
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zsyr2k_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
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
    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);
    CORE_zsyr2k(uplo, trans,
                 n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zsyr2k_cuda_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
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
    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);

    RUNTIME_getStream(stream);

    CUDA_zsyr2k( uplo, trans,
                 n, k, &alpha, A, ldA, B, ldB, &beta, C, ldC,
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
CODELETS(zsyr2k, 3, cl_zsyr2k_cpu_func, cl_zsyr2k_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zsyr2k(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_trans_t trans,
                       int n, int k, int nb,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int ldA,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldB,
                       CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldC)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zsyr2k;
    void (*callback)(void*) = options->profiling ? cl_zsyr2k_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,      &uplo,                sizeof(int),
        STARPU_VALUE,     &trans,                sizeof(int),
        STARPU_VALUE,         &n,                        sizeof(int),
        STARPU_VALUE,         &k,                        sizeof(int),
        STARPU_VALUE,     &alpha,         sizeof(CHAMELEON_Complex64_t),
        STARPU_R,                 RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,                 RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_VALUE,      &beta,         sizeof(CHAMELEON_Complex64_t),
        STARPU_RW,                 RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zsyr2k",
#endif
        0);
    (void)ldC;
    (void)ldB;
    (void)ldA;
}
