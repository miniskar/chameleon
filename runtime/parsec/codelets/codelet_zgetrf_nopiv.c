/**
 *
 * @file parsec/codelet_zgetrf_nopiv.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv PaRSEC codelet
 *
 * @version 1.2.0
 * @author Reazul Hoque
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgetrf_nopiv_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    int ib;
    CHAMELEON_Complex64_t *A;
    int lda;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &ib, &A, &lda, &iinfo, &sequence, &request );

    CORE_zgetrf_nopiv( m, n, ib, A, lda, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zgetrf_nopiv(const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgetrf_nopiv_parsec, options->priority,  "getrf_nopiv",
        sizeof(int),                 &m,                          VALUE,
        sizeof(int),                 &n,                          VALUE,
        sizeof(int),                 &ib,                         VALUE,
        PASSED_BY_REF,               RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(int),                 &iinfo,                      VALUE,
        sizeof(RUNTIME_sequence_t*), &(options->sequence),        VALUE,
        sizeof(RUNTIME_request_t*),  &(options->request),         VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
