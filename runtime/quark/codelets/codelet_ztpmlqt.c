/**
 *
 * @file quark/codelet_ztpmlqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpmlqt Quark codelet
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @date 2018-01-31
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static void
CORE_ztpmlqt_quark( Quark *quark )
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

    quark_unpack_args_16( quark, side, trans, M, N, K, L, ib,
                          V, ldv, T, ldt, A, lda, B, ldb, WORK );

    CORE_ztpmlqt( side, trans, M, N, K, L, ib,
                  V, ldv, T, ldt, A, lda, B, ldb, WORK );
}

void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int M, int N, int K, int L, int ib, int nb,
                         const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TPMLQT;

    int shapeV = ( L == 0 ) ? 0 : (QUARK_REGION_L | QUARK_REGION_D);

    QUARK_Insert_Task(
        opt->quark, CORE_ztpmlqt_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &side,  VALUE,
        sizeof(int),              &trans, VALUE,
        sizeof(int),                     &M,     VALUE,
        sizeof(int),                     &N,     VALUE,
        sizeof(int),                     &K,     VALUE,
        sizeof(int),                     &L,     VALUE,
        sizeof(int),                     &ib,    VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR( V, CHAMELEON_Complex64_t, Vm, Vn ), INPUT | shapeV,
        sizeof(int),                     &ldv,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), INPUT,
        sizeof(int),                     &ldt,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), INOUT,
        sizeof(int),                     &lda,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), INOUT | LOCALITY,
        sizeof(int),                     &ldb,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  NULL, SCRATCH,
        0);
}
