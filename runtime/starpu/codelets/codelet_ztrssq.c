/**
 *
 * @file starpu/codelet_ztrssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrssq StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztrssq_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_diag_t diag;
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int ldA;
    double *SCALESUMSQ;

    A          = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    SCALESUMSQ = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    starpu_codelet_unpack_args(cl_arg, &uplo, &diag, &m, &n);
    CORE_ztrssq( uplo, diag, m, n, A, ldA, &SCALESUMSQ[0], &SCALESUMSQ[1]);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztrssq, 2, cl_ztrssq_cpu_func)

void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An, int ldA,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    struct starpu_codelet *codelet = &cl_ztrssq;
    void (*callback)(void*) = options->profiling ? cl_ztrasm_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(SCALESUMSQ, SCALESUMSQm, SCALESUMSQn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,                      sizeof(int),
        STARPU_VALUE,    &diag,                      sizeof(int),
        STARPU_VALUE,    &m,                         sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_R,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_RW,       RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztrssq",
#endif
        0);
    (void)ldA;
}
