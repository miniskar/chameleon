/**
 *
 * @file parsec/codelet_zgessq.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessq PaRSEC codelet
 *
 * @version 0.9.2
 * @author Reazul Hoque
 * @date 2015-11-04
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgessq_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_store_t storev;
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *SCALESUMSQ;

    parsec_dtd_unpack_args(
        this_task, &storev, &m, &n, &A, &lda, &SCALESUMSQ );

    CORE_zgessq( storev, m, n, A, lda, SCALESUMSQ );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zgessq( const RUNTIME_option_t *options,
                        cham_store_t storev, int m, int n,
                        const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgessq_parsec, options->priority, "gessq",
        sizeof(cham_store_t), &storev,       VALUE,
        sizeof(int),          &m,            VALUE,
        sizeof(int),          &n,            VALUE,
        PASSED_BY_REF,   RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),          &lda,          VALUE,
        PASSED_BY_REF,   RTBLKADDR( SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn ), INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
