/**
 *
 * @file starpu/codelet_ztrtri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
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
static void cl_ztrtri_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_diag_t diag;
    int N;
    CHAM_tile_t *tileA;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info = 0;

    tileA = cti_interface_get(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &diag, &N, &iinfo, &sequence, &request);
    TCORE_ztrtri(uplo, diag, N, tileA, &info);

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztrtri, 1, cl_ztrtri_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo )
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_ztrtri;
    void (*callback)(void*) = options->profiling ? cl_ztrtri_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,              sizeof(int),
        STARPU_VALUE,    &diag,              sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_RW,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,    &iinfo,             sizeof(int),
        STARPU_VALUE,    &(options->sequence),       sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,    &(options->request),        sizeof(RUNTIME_request_t*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztrtri",
#endif
        0);
}
