/**
 *
 * @file starpu/codelet_zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
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
static void cl_zgetrf_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int ldA;
    int *IPIV;
    cham_bool_t check_info;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info = 0;

    A = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    ldA = STARPU_MATRIX_GET_LD( descr[0] );

    starpu_codelet_unpack_args(cl_arg, &m, &n, &IPIV, &check_info, &iinfo, &sequence, &request);
    CORE_zgetrf( m, n, A, ldA, IPIV, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf, 1, cl_zgetrf_cpu_func)

void INSERT_TASK_zgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An, int ldA,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo )
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,             &m,                        sizeof(int),
        STARPU_VALUE,             &n,                        sizeof(int),
        STARPU_RW,                     RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,                  &IPIV,                      sizeof(int*),
        STARPU_VALUE,    &check_info,                sizeof(cham_bool_t),
        STARPU_VALUE,         &iinfo,                        sizeof(int),
        STARPU_VALUE,    &(options->sequence),       sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,    &(options->request),        sizeof(RUNTIME_request_t*),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf",
#endif
        0);
    (void)ldA;
}
