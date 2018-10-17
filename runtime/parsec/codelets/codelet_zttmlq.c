/**
 *
 * @file parsec/codelet_zttmlq.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zttmlq PaRSEC codelet
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
CORE_zttmlq_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_trans_t trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A1;
    int lda1;
    CHAMELEON_Complex64_t *A2;
    int lda2;
    CHAMELEON_Complex64_t *V;
    int ldv;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib, &A1, &lda1, &A2, &lda2, &V, &ldv, &T, &ldt, &WORK, &ldwork );

    CORE_zttmlq( side, trans, m1, n1, m2, n2, k, ib, A1, lda1,
                 A2, lda2, V, ldv, T, ldt, WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zttmlq(const RUNTIME_option_t *options,
                       cham_side_t side, cham_trans_t trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                       const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                       const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldt)
{
    int ldwork = side == ChamLeft ? ib : nb;

    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zttmlq_parsec, options->priority, "ttmlq",
        sizeof(int),     &side,                      VALUE,
        sizeof(int),     &trans,                     VALUE,
        sizeof(int),            &m1,                        VALUE,
        sizeof(int),            &n1,                        VALUE,
        sizeof(int),            &m2,                        VALUE,
        sizeof(int),            &n2,                        VALUE,
        sizeof(int),            &k,                         VALUE,
        sizeof(int),            &ib,                        VALUE,
        PASSED_BY_REF,          RTBLKADDR( A1, CHAMELEON_Complex64_t, A1m, A1n ), chameleon_parsec_get_arena_index( A1 ) | INOUT,
        sizeof(int),            &lda1,                      VALUE,
        PASSED_BY_REF,          RTBLKADDR( A2, CHAMELEON_Complex64_t, A2m, A2n ), chameleon_parsec_get_arena_index( A2 ) | INOUT | AFFINITY,
        sizeof(int),            &lda2,                      VALUE,
        PASSED_BY_REF,          RTBLKADDR( V, CHAMELEON_Complex64_t, Vm, Vn ), chameleon_parsec_get_arena_index( V ) | INPUT,
        sizeof(int),            &ldv,                       VALUE,
        PASSED_BY_REF,          RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | INPUT,
        sizeof(int),            &ldt,                       VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,    NULL,           SCRATCH,
        sizeof(int),            &ldwork,                    VALUE,
        PARSEC_DTD_ARG_END );
}
