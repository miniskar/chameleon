/**
 *
 * @file pzgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
#define B(m, n) B,  m,  n
#define C(m, n) C,  m,  n
/**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 */
void chameleon_pzgemm(cham_trans_t transA, cham_trans_t transB,
                         CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                         CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int m, n, k;
    int ldam, ldak, ldbn, ldbk, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for (m = 0; m < C->mt; m++) {
        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
        ldcm = BLKLDD(C, m);
        for (n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
            /*
             *  A: ChamNoTrans / B: ChamNoTrans
             */
            if (transA == ChamNoTrans) {
                ldam = BLKLDD(A, m);
                if (transB == ChamNoTrans) {
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k), ldam,  /* lda * Z */
                                   B(k, n), ldbk,  /* ldb * Y */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
                /*
                 *  A: ChamNoTrans / B: Cham[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(B, n);
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k), ldam,  /* lda * Z */
                                   B(n, k), ldbn,  /* ldb * Z */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
            }
            /*
             *  A: Cham[Conj]Trans / B: ChamNoTrans
             */
            else {
                if (transB == ChamNoTrans) {
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m), ldak,  /* lda * X */
                                   B(k, n), ldbk,  /* ldb * Y */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
                /*
                 *  A: Cham[Conj]Trans / B: Cham[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(B, n);
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m), ldak,  /* lda * X */
                                   B(n, k), ldbn,  /* ldb * Z */
                            zbeta, C(m, n), ldcm); /* ldc * Y */
                    }
                }
            }
            RUNTIME_data_flush( sequence, C(m, n) );
        }
        if (transA == ChamNoTrans) {
            for (k = 0; k < A->nt; k++) {
                RUNTIME_data_flush( sequence, A(m, k) );
            }
        } else {
            for (k = 0; k < A->mt; k++) {
                RUNTIME_data_flush( sequence, A(k, m) );
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
