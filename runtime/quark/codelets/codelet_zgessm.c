/**
 *
 * @file codelet_zgessm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessm Quark codelet
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
#include "coreblas/cblas.h"
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zgessm_quark(Quark *quark)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    CHAMELEON_Complex64_t *L;
    int ldl;
    CHAMELEON_Complex64_t *D;
    int ldd;
    CHAMELEON_Complex64_t *A;
    int lda;

    quark_unpack_args_11(quark, m, n, k, ib, IPIV, L, ldl, D, ldd, A, lda);
    CORE_zgessm(m, n, k, ib, IPIV, D, ldd, A, lda);
}

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zgessm applies the factors L computed by CORE_zgetrf_incpiv to
 *  a complex M-by-N tile A.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] K
 *         The number of columns of the tile L. K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_zgetrf_incpiv.
 *
 * @param[in] L
 *         The M-by-K lower triangular tile.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,M).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the application of L.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *         \retval CHAMELEON_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */
void INSERT_TASK_zgessm(const RUNTIME_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       const CHAM_desc_t *L, int Lm, int Ln, int ldl,
                       const CHAM_desc_t *D, int Dm, int Dn, int ldd,
                       const CHAM_desc_t *A, int Am, int An, int lda)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GESSM;
    QUARK_Insert_Task(opt->quark, CORE_zgessm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(L, CHAMELEON_Complex64_t, Lm, Ln),             INPUT | QUARK_REGION_L,
        sizeof(int),                        &ldl,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(D, CHAMELEON_Complex64_t, Dm, Dn),             INPUT | QUARK_REGION_L,
        sizeof(int),                        &ldd,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}
