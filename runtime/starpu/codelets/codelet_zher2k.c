/**
 *
 * @file starpu/codelet_zher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k StarPU codelet
 *
 * @version 1.3.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Loris Lucido
 * @date 2023-07-06
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zher2k_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    double beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);
    TCORE_zher2k(uplo, trans,
                n, k, alpha, tileA, tileB, beta, tileC);
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zher2k_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    cuDoubleComplex alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    double beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);

    CUDA_zher2k( uplo, trans,
                 n, k,
                 &alpha, tileA->mat, tileA->ld,
                         tileB->mat, tileB->ld,
                 &beta,  tileC->mat, tileC->ld,
                 handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void cl_zher2k_hip_func(void *descr[], void *cl_arg)
{
    hipblasHandle_t handle = starpu_hipblas_get_local_handle();
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    hipblasDoubleComplex alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    double beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);

    HIP_zher2k( uplo, trans,
                 n, k,
                 &alpha, tileA->mat, tileA->ld,
                         tileB->mat, tileB->ld,
                 &beta,  tileC->mat, tileC->ld,
                 handle );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( zher2k, cl_zher2k_cpu_func, cl_zher2k_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( zher2k, cl_zher2k_cpu_func, cl_zher2k_cuda_func, STARPU_CUDA_ASYNC )
#endif

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void
INSERT_TASK_zher2k( const RUNTIME_option_t *options,
                    cham_uplo_t uplo, cham_trans_t trans,
                    int n, int k, int nb,
                    CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                 const CHAM_desc_t *B, int Bm, int Bn,
                    double beta,                 const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    (void)nb;
    struct starpu_codelet *codelet = &cl_zher2k;
    void (*callback)(void*) = options->profiling ? cl_zher2k_callback : NULL;
    int accessC = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,      &uplo,                sizeof(int),
        STARPU_VALUE,     &trans,                sizeof(int),
        STARPU_VALUE,         &n,                        sizeof(int),
        STARPU_VALUE,         &k,                        sizeof(int),
        STARPU_VALUE,     &alpha,         sizeof(CHAMELEON_Complex64_t),
        STARPU_R,                 RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R,                 RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_VALUE,      &beta,                     sizeof(double),
        accessC,                  RTBLKADDR(C, ChamComplexDouble, Cm, Cn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zher2k",
#endif
        0);
}
