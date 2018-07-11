/**
 *
 * @file codelet_zsymm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm PaRSEC codelet
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
CORE_zsymm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *C;
    int LDC;

    parsec_dtd_unpack_args(
        this_task, &side, &uplo, &M, &N, &alpha, &A, &LDA, &B, &LDB, &beta, &C, &LDC );

    CORE_zsymm( side, uplo, M, N,
                alpha, A, LDA,
                       B, LDB,
                beta,  C, LDC);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zsymm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                      const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                      CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zsymm_parsec, options->priority, "symm",
        sizeof(int),     &side,              VALUE,
        sizeof(int),     &uplo,              VALUE,
        sizeof(int),            &m,                 VALUE,
        sizeof(int),            &n,                 VALUE,
        sizeof(CHAMELEON_Complex64_t),  &alpha,         VALUE,
        PASSED_BY_REF,          RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),            &lda,               VALUE,
        PASSED_BY_REF,          RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), morse_parsec_get_arena_index( B ) | INPUT,
        sizeof(int),            &ldb,               VALUE,
        sizeof(CHAMELEON_Complex64_t),  &beta,          VALUE,
        PASSED_BY_REF,          RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), morse_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int),            &ldc,               VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
