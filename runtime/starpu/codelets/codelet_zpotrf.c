/**
 *
 * @file starpu/codelet_zpotrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf StarPU codelet
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
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 */
void INSERT_TASK_zpotrf(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zpotrf;
    void (*callback)(void*) = options->profiling ? cl_zpotrf_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,                      sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_RW,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        STARPU_VALUE,    &(options->sequence),       sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,    &(options->request),        sizeof(RUNTIME_request_t*),
        /* STARPU_SCRATCH,   options->ws_worker, */
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zpotrf",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zpotrf_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info = 0;

    A = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &lda, &iinfo, &sequence, &request);
    CORE_zpotrf(uplo, n, A, lda, &info);

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zpotrf, 1, cl_zpotrf_cpu_func)

