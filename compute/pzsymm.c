/**
 *
 * @file pzsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm parallel algorithm
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

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
/**
 *  Parallel tile symmetric matrix-matrix multiplication - dynamic scheduling
 */
void chameleon_pzsymm(cham_side_t side, cham_uplo_t uplo,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                          CHAMELEON_Complex64_t beta, CHAM_desc_t *C,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int ldak, ldam, ldan, ldbk, ldbm, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    /*
     *  ChamLeft
     */
    if (side == ChamLeft) {
        for (k = 0; k < C->mt; k++) {
            tempkm = k == C->mt-1 ? C->m-k*C->mb : C->mb;
            ldak = BLKLDD(A, k);
            ldbk = BLKLDD(B, k);
            zbeta = k == 0 ? beta : zone;

            for (n = 0; n < C->nt; n++) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                for (m = 0; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    ldam = BLKLDD(A, m);
                    ldcm = BLKLDD(C, m);

                    /*
                     *  ChamLeft / ChamLower
                     */
                    if (uplo == ChamLower) {
                        if (k < m) {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(m, k), ldam,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm);
                        }
                        else {
                            if (k == m) {
                                INSERT_TASK_zsymm(
                                    &options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm);
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    &options,
                                    ChamTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(k, m), ldak,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm);
                            }
                        }
                    }
                    /*
                     *  ChamLeft / ChamUpper
                     */
                    else {
                        if (k < m) {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm);
                        }
                        else {
                            if (k == m) {
                                INSERT_TASK_zsymm(
                                    &options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm);
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    &options,
                                    ChamNoTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(m, k), ldam,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm);
                            }
                        }
                    }
                }
                RUNTIME_data_flush( sequence, B(k, n) );
            }
            if (uplo == ChamLower) {
                for (n = 0; n <= k; n++) {
                    RUNTIME_data_flush( sequence, A(k, n) );
                }
            }
            else {
                for (m = 0; m <= k; m++) {
                    RUNTIME_data_flush( sequence, A(m, k) );
                }
            }
        }
    }
    /*
     *  ChamRight
     */
    else {
        for (k = 0; k < C->nt; k++) {
            tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
            ldak = BLKLDD(A, k);
            zbeta = k == 0 ? beta : zone;

            for (m = 0; m < C->mt; m++) {
                tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                ldbm = BLKLDD(B, m);
                ldcm = BLKLDD(C, m);

                for (n = 0; n < C->nt; n++) {
                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                    ldan = BLKLDD(A, n);

                    /*
                     *  ChamRight / ChamLower
                     */
                    if (uplo == ChamLower) {
                        if (k < n) {
                           INSERT_TASK_zgemm(
                               &options,
                               ChamNoTrans, ChamTrans,
                               tempmm, tempnn, tempkn, A->mb,
                               alpha, B(m, k), ldbm,
                                      A(n, k), ldan,
                               zbeta, C(m, n), ldcm);
                        }
                        else {
                            if (k == n) {
                               INSERT_TASK_zsymm(
                                   &options,
                                   ChamRight, uplo,
                                   tempmm, tempnn, A->mb,
                                   alpha, A(k, k), ldak,
                                          B(m, k), ldbm,
                                   zbeta, C(m, n), ldcm);
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    &options,
                                    ChamNoTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k), ldbm,
                                           A(k, n), ldak,
                                    zbeta, C(m, n), ldcm);
                            }
                        }
                    }
                    /*
                     *  ChamRight / ChamUpper
                     */
                    else {
                        if (k < n) {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,
                                       A(k, n), ldak,
                                zbeta, C(m, n), ldcm);
                        }
                        else {
                            if (k == n) {
                                INSERT_TASK_zsymm(
                                    &options,
                                    ChamRight, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,
                                           B(m, k), ldbm,
                                    zbeta, C(m, n), ldcm);
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    &options,
                                    ChamNoTrans, ChamTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k), ldbm,
                                           A(n, k), ldan,
                                    zbeta, C(m, n), ldcm);
                            }
                        }
                    }
                }
                RUNTIME_data_flush( sequence, B(m, k) );
            }
            if (uplo == ChamLower) {
                for (n = 0; n <= k; n++) {
                    RUNTIME_data_flush( sequence, A(k, n) );
                }
            }
            else {
                for (m = 0; m <= k; m++) {
                    RUNTIME_data_flush( sequence, A(m, k) );
                }
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
