/**
 *
 * @file pzpotrimm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrimm parallel algorithm
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Ali M Charara
 * @date 2014-11-16
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 */
void chameleon_pzpotrimm(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int ldbm, ldcm;
    int ldak, ldam, ldan;
    int tempkm, tempmm, tempnn, tempkn;

    CHAMELEON_Complex64_t alpha = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t beta  = (CHAMELEON_Complex64_t) 0.0;
    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;


    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    /*
     *  ChamLower
     */
    if (uplo == ChamLower) {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);

            INSERT_TASK_zpotrf(
                &options,
                ChamLower, tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, ChamLower, ChamConjTrans, ChamNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                ldan = BLKLDD(A, n);
                INSERT_TASK_zherk(
                    &options,
                    ChamLower, ChamNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k), ldan,
                     1.0, A(n, n), ldan);

                for (m = n+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k), ldam,
                               A(n, k), ldan,
                        zone,  A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(n, k) );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, A->nt + k);

            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, uplo, ChamNoTrans, ChamNonUnit,
                    tempmm, tempkn, A->mb,
                    mzone, A(k, k), ldak,
                           A(m, k), ldam);
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                for (n = 0; n < k; n++) {
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamNoTrans,
                        tempmm, A->nb, tempkn, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(m, k) );
            }
            for (n = 0; n < k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, uplo, ChamNoTrans, ChamNonUnit,
                    tempkn, A->nb, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            INSERT_TASK_ztrtri(
                &options,
                uplo, ChamNonUnit,
                tempkn, A->mb,
                A(k, k), ldak, A->nb*k);

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, 2*A->nt + k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for(n = 0; n < k; n++) {
                ldan = BLKLDD(A, n);
                INSERT_TASK_zherk(
                    &options,
                    uplo, ChamConjTrans,
                    A->mb, tempkm, A->mb,
                    1.0, A(k, n), ldak,
                    1.0, A(n, n), ldan);

                for(m = n+1; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        A->mb, A->nb, tempkm, A->mb,
                        1.0, A(k, m), ldak,
                             A(k, n), ldak,
                        1.0, A(m, n), ldam);
                }
            }
            for (n = 0; n < k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
                INSERT_TASK_ztrmm(
                    &options,
                    ChamLeft, uplo, ChamConjTrans, ChamNonUnit,
                    tempkm, A->nb, A->mb,
                    1.0, A(k, k), ldak,
                         A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            INSERT_TASK_zlauum(
                &options,
                uplo, tempkm, A->mb,
                A(k, k), ldak);

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZSYMM Right / Lower
         */
        for (k = 0; k < C->nt; k++) {
            RUNTIME_iteration_push(chamctxt, 3*A->nt + k);

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

                    if (k < n) {
                       INSERT_TASK_zgemm(
                           &options,
                           ChamNoTrans, ChamTrans,
                           tempmm, tempnn, tempkn, A->mb,
                           alpha, B(m, k), ldbm,  /* ldbm * K */
                                  A(n, k), ldan,  /* ldan * K */
                           zbeta, C(m, n), ldcm); /* ldcm * Y */
                    }
                    else {
                        if (k == n) {
                           INSERT_TASK_zsymm(
                               &options,
                               ChamRight, uplo,
                               tempmm, tempnn, A->mb,
                               alpha, A(k, k), ldak,  /* ldak * Y */
                                      B(m, k), ldbm,  /* ldbm * Y */
                               zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                        else {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldbm * K */
                                       A(k, n), ldak,  /* ldak * Y */
                                zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                    }
                }
                RUNTIME_data_flush( sequence, B(m, k) );
            }
            for (n = 0; n <= k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /*
     *  ChamUpper
     */
    else {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            INSERT_TASK_zpotrf(
                &options,
                ChamUpper,
                tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n - n*A->nb : A->nb;
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, ChamUpper, ChamConjTrans, ChamNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                ldam = BLKLDD(A, m);

                INSERT_TASK_zherk(
                    &options,
                    ChamUpper, ChamConjTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam);

                for (n = m+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m), ldak,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(k, m) );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, A->nt + k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, uplo, ChamNoTrans, ChamNonUnit,
                    tempkm, tempnn, A->mb,
                    mzone, A(k, k), ldak,
                           A(k, n), ldak);
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                for (m = 0; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamNoTrans,
                        A->mb, tempnn, tempkm, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                RUNTIME_data_flush( sequence, A(k, n) );
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                RUNTIME_data_flush( sequence, A(m, k) );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, uplo, ChamNoTrans, ChamNonUnit,
                    A->mb, tempkm, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            INSERT_TASK_ztrtri(
                &options,
                uplo, ChamNonUnit,
                tempkm, A->mb,
                A(k, k), ldak, A->mb*k);

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, 2*A->nt + k);

            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);

            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                INSERT_TASK_zherk(
                    &options,
                    uplo, ChamNoTrans,
                    A->mb, tempkn, A->mb,
                    1.0, A(m, k), ldam,
                    1.0, A(m, m), ldam);

                for (n = m+1; n < k; n++){
                    ldan = BLKLDD(A, n);
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        A->mb, A->nb, tempkn, A->mb,
                        1.0, A(m, k), ldam,
                             A(n, k), ldan,
                        1.0, A(m, n), ldam);
                }
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                RUNTIME_data_flush( sequence, A(m, k) );
                INSERT_TASK_ztrmm(
                    &options,
                    ChamRight, uplo, ChamConjTrans, ChamNonUnit,
                    A->mb, tempkn, A->mb,
                    1.0, A(k, k), ldak,
                         A(m, k), ldam);
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            INSERT_TASK_zlauum(
                &options,
                uplo, tempkn, A->mb,
                A(k, k), ldak);

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZSYMM Right / Upper
         */
        for (k = 0; k < C->nt; k++) {
            RUNTIME_iteration_push(chamctxt, 3*A->nt + k);

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

                    if (k < n) {
                        INSERT_TASK_zgemm(
                            &options,
                            ChamNoTrans, ChamNoTrans,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, B(m, k), ldbm,  /* ldbm * K */
                                   A(k, n), ldak,  /* ldak * Y */
                            zbeta, C(m, n), ldcm); /* ldcm * Y */
                    }
                    else {
                        if (k == n) {
                            INSERT_TASK_zsymm(
                                &options,
                                ChamRight, uplo,
                                tempmm, tempnn, A->mb,
                                alpha, A(k, k), ldak,  /* ldak * Y */
                                       B(m, k), ldbm,  /* ldbm * Y */
                                zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                        else {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, ChamTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldbm * K */
                                       A(n, k), ldan,  /* ldan * K */
                                zbeta, C(m, n), ldcm); /* ldcm * Y */
                        }
                    }
                }
                RUNTIME_data_flush( sequence, B(m, k) );
            }
            for (m = 0; m <= k; m++) {
                RUNTIME_data_flush( sequence, A(m, k) );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }

    RUNTIME_options_finalize(&options, chamctxt);
}
