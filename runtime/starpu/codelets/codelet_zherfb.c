/**
 *
 * @file starpu/codelet_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherfb StarPU codelet
 *
 * @version 0.9.2
 * @author Hatem Ltaief
 * @author Lucas Barros de Assis
 * @date 2016-12-09
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zherfb_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int n;
    int k;
    int ib;
    int nb;
    const CHAMELEON_Complex64_t *A;
    int ldA;
    const CHAMELEON_Complex64_t *T;
    int ldT;
    CHAMELEON_Complex64_t *C;
    int ldC;
    CHAMELEON_Complex64_t *WORK;
    int ldWORK;

    A    = (const CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    T    = (const CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    WORK = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]); /* ib * nb */

    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldT = STARPU_MATRIX_GET_LD( descr[1] );
    ldC = STARPU_MATRIX_GET_LD( descr[2] );

    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &k, &ib, &nb, &ldWORK);

    CORE_zherfb(uplo, n, k, ib, nb, A, ldA, T, ldT, C, ldC, WORK, ldWORK);
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zherfb_cuda_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int n;
    int k;
    int ib;
    int nb;
    const cuDoubleComplex *A;
    int ldA;
    const cuDoubleComplex *T;
    int ldT;
    cuDoubleComplex *C;
    int ldC;
    cuDoubleComplex *WORK;
    int ldWORK;

    RUNTIME_getStream(stream);

    A    = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    T    = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C    = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    WORK = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[3]); /* ib * nb */

    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldT = STARPU_MATRIX_GET_LD( descr[1] );
    ldC = STARPU_MATRIX_GET_LD( descr[2] );

    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &k, &ib, &nb, &ldWORK);

    CUDA_zherfb( uplo, n, k, ib, nb, A, ldA, T, ldT, C, ldC, WORK, ldWORK, stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zherfb, 4, cl_zherfb_cpu_func, cl_zherfb_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zherfb(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An, int ldA,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldT,
                       const CHAM_desc_t *C, int Cm, int Cn, int ldC)
{
    struct starpu_codelet *codelet = &cl_zherfb;
    void (*callback)(void*) = options->profiling ? cl_zherfb_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,              sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_VALUE,    &nb,                sizeof(int),
        STARPU_R,         RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,         RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        STARPU_RW,        RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &nb,                sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zherfb",
#endif
        0);
    (void)ldC;
    (void)ldT;
    (void)ldA;
}
