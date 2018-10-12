/**
 *
 * @file quark/codelet_zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
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

void CORE_zgeadd_quark(Quark *quark)
{
    cham_trans_t trans;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *B;
    int LDB;

    quark_unpack_args_9(quark, trans, M, N, alpha, A, LDA, beta, B, LDB);
    CORE_zgeadd(trans, M, N, alpha, A, LDA, beta, B, LDB);
    return;
}

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  INSERT_TASK_zgeadd adds two general matrices together as in PBLAS pzgeadd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = ChamNoTrans:   op(A) = A
 *          = ChamTrans:     op(A) = A'
 *          = ChamConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices op(A) and B.
 *
 * @param[in] N
 *          Number of columns of the matrices op(A) and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans = ChamNoTrans, LDA-by-M
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,k), with k=M, if
 *          trans = ChamNoTrans, and k=N otherwise.
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *          On exit, B = alpha * op(A) + beta * B
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
void INSERT_TASK_zgeadd(const RUNTIME_option_t *options,
                       cham_trans_t trans, int m, int n, int nb,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                       CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GEADD;
    QUARK_Insert_Task(opt->quark, CORE_zgeadd_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                 &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha, VALUE,
        sizeof(CHAMELEON_Complex64_t)*lda*n,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(CHAMELEON_Complex64_t),         &beta,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*ldb*n,    RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);

    (void)nb;
}
