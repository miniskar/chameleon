/**
 *
 * @file starpu/codelet_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk StarPU codelet
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
static void cl_zherk_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    double alpha;
    CHAMELEON_Complex64_t *A;
    int ldA;
    double beta;
    CHAMELEON_Complex64_t *C;
    int ldC;

    A = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    C = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldC = STARPU_MATRIX_GET_LD( descr[1] );

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);
    CORE_zherk(uplo, trans,
        n, k,
        alpha, A, ldA,
        beta, C, ldC);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zherk_cuda_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    double alpha;
    const cuDoubleComplex *A;
    int ldA;
    double beta;
    cuDoubleComplex *C;
    int ldC;

    A = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);

    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldC = STARPU_MATRIX_GET_LD( descr[1] );

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);

    RUNTIME_getStream(stream);

    CUDA_zherk(
        uplo, trans,
        n, k,
        &alpha, A, ldA,
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
CODELETS(zherk, 2, cl_zherk_cpu_func, cl_zherk_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zherk(const RUNTIME_option_t *options,
                      cham_uplo_t uplo, cham_trans_t trans,
                      int n, int k, int nb,
                      double alpha, const CHAM_desc_t *A, int Am, int An, int ldA,
                      double beta, const CHAM_desc_t *C, int Cm, int Cn, int ldC)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zherk;
    void (*callback)(void*) = options->profiling ? cl_zherk_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,              sizeof(int),
        STARPU_VALUE,    &trans,             sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &alpha,             sizeof(double),
        STARPU_R,         RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,    &beta,              sizeof(double),
        STARPU_RW,        RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zherk",
#endif
        0);
    (void)ldC;
    (void)ldA;
}
