/**
 *
 * @file codelet_ztrasm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrasm PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztrasm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_store_t storev;
    cham_uplo_t uplo;
    cham_diag_t diag;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *work;

    parsec_dtd_unpack_args(
        this_task, &storev, &uplo, &diag, &M, &N, &A, &lda, &work );

    CORE_ztrasm( storev, uplo, diag, M, N, A, lda, work );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztrasm(const RUNTIME_option_t *options,
                       cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztrasm_parsec, options->priority, "trasm",
        sizeof(int),     &storev,                VALUE,
        sizeof(int),     &uplo,                  VALUE,
        sizeof(int),     &diag,                  VALUE,
        sizeof(int),            &M,                     VALUE,
        sizeof(int),            &N,                     VALUE,
        PASSED_BY_REF,          RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),            &lda,                   VALUE,
        PASSED_BY_REF,          RTBLKADDR( B, double, Bm, Bn ),     INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
