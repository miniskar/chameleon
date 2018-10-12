/**
 *
 * @file starpu/codelet_zhessq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhessq StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void INSERT_TASK_zhessq( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, int n,
                        const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    struct starpu_codelet *codelet = &cl_zhessq;
    void (*callback)(void*) = options->profiling ? cl_zgessq_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(SCALESUMSQ, SCALESUMSQm, SCALESUMSQn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,                       sizeof(int),
        STARPU_VALUE,    &n,                          sizeof(int),
        STARPU_R,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                        sizeof(int),
        STARPU_RW,       RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zhessq",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zhessq_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *SCALESUMSQ;

    A          = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    SCALESUMSQ = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &lda);
    CORE_zhessq( uplo, n, A, lda, &SCALESUMSQ[0], &SCALESUMSQ[1] );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhessq, 2, cl_zhessq_cpu_func)
