/**
 *
 * @file codelet_ztpmqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpmqrt PaRSEC codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztpmqrt_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    const CHAMELEON_Complex64_t *V;
    int ldv;
    const CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;
    CHAMELEON_Complex64_t *WORK;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &M, &N, &K, &L, &ib, &V, &ldv, &T, &ldt, &A, &lda, &B, &ldb, &WORK );

    CORE_ztpmqrt( side, trans, M, N, K, L, ib,
                  V, ldv, T, ldt, A, lda, B, ldb, WORK );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int M, int N, int K, int L, int ib, int nb,
                         const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztpmqrt_parsec, options->priority, "tpmqrt",
        sizeof(int), &side,  VALUE,
        sizeof(int), &trans, VALUE,
        sizeof(int),        &M,     VALUE,
        sizeof(int),        &N,     VALUE,
        sizeof(int),        &K,     VALUE,
        sizeof(int),        &L,     VALUE,
        sizeof(int),        &ib,    VALUE,
        PASSED_BY_REF,       RTBLKADDR( V, CHAMELEON_Complex64_t, Vm, Vn ), chameleon_parsec_get_arena_index( V ) | INPUT,
        sizeof(int),        &ldv,   VALUE,
        PASSED_BY_REF,       RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | INPUT,
        sizeof(int),        &ldt,   VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT,
        sizeof(int),        &lda,   VALUE,
        PASSED_BY_REF,       RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int),        &ldb,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb, NULL, SCRATCH,
        PARSEC_DTD_ARG_END );
}
