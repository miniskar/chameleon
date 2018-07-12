/**
 *
 * @file codelet_zlaset2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset2 Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlaset2_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;

    quark_unpack_args_6(quark, uplo, M, N, alpha, A, LDA);
    CORE_zlaset2(uplo, M, N, alpha, A, LDA);
}

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zlaset2 - Sets the elements of the matrix A to alpha.
 *  Not LAPACK compliant! Read below.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = ChamUpper: STRICT Upper part of A is set to alpha;
 *          = ChamLower: STRICT Lower part of A is set to alpha;
 *          = ChamUpperLower: ALL elements of A are set to alpha.
 *          Not LAPACK Compliant.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set to alpha accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 */
void INSERT_TASK_zlaset2(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int M, int N,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int LDA)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LASET;
    QUARK_Insert_Task(opt->quark, CORE_zlaset2_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha, VALUE,
        sizeof(CHAMELEON_Complex64_t)*M*N,     RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}
