/**
 *
 * @file pzunmqrrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqrrh parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A, (m), (n)
#define B(m,n)  B, (m), (n)
#define T(m,n)  T, (m), (n)
#define T2(m,n) T, (m), ((n)+A->nt)
#define D(m,n)  D, (m), (n)

/**
 *  Parallel application of Q using tile V - QR factorization (reduction
 *  Householder) - dynamic scheduling
 */
void chameleon_pzunmqrrh( int genD, int BS, cham_side_t side, cham_trans_t trans,
                          CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int K, M, RD, lastRD;
    int ldaM, ldam, ldan, ldaMRD, lddM;
    int ldbM, ldbm, ldbMRD;
    int tempMm, tempkn, tempnn, tempmm, tempMRDm, tempkmin;
    int ib;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmqr = A->nb * ib
     * ztsmqr = A->nb * ib
     * zttmqr = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);
    if (side == ChamLeft ) {
        if (trans == ChamConjTrans) {
            /*
             *  ChamLeft / ChamConjTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = chameleon_min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    lddM = BLKLDD(D, M);
                    ldbM = BLKLDD(B, M);
                    if ( genD ) {
                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempMm, tempkmin, A->nb,
                            A(M, k), ldaM,
                            D(M, k), lddM );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempMm, tempkmin,
                            0., 1.,
                            D(M, k), lddM );
#endif
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempMm, tempnn, tempkmin, ib, T->nb,
                            D(M, k), lddM,
                            T(M, k), T->mb,
                            B(M, n), ldbM);
                    }
                    RUNTIME_data_flush( sequence, D(M, k) );
                    RUNTIME_data_flush( sequence, T(M, k) );

                    for (m = M+1; m < chameleon_min(M+BS, A->mt); m++) {
                        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                        ldbm = BLKLDD(B, m);
                        ldam = BLKLDD(A, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                            RUNTIME_data_migrate( sequence, B(M, n),
                                                  B->get_rankof( B, m, n ) );
                            RUNTIME_data_migrate( sequence, B(m, n),
                                                  B->get_rankof( B, m, n ) );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(m, k), ldam,
                                T(m, k), T->mb,
                                B(M, n), ldbM,
                                B(m, n), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(m, k) );
                        RUNTIME_data_flush( sequence, T(m, k) );
                    }
                }
                for (RD = BS; RD < A->mt-k; RD *= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldbM   = BLKLDD(B, M   );
                        ldbMRD = BLKLDD(B, M+RD);
                        ldaMRD = BLKLDD(A, M+RD);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                            RUNTIME_data_migrate( sequence, B(M, n),
                                                  B->get_rankof( B, M+RD, n ) );
                            RUNTIME_data_migrate( sequence, B(M+RD, n),
                                                  B->get_rankof( B, M+RD, n ) );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempMRDm, tempnn, tempkn, tempMRDm, ib, T->nb,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb,
                                B (M,    n), ldbM,
                                B (M+RD, n), ldbMRD);
                        }
                        RUNTIME_data_flush( sequence, A (M+RD, k) );
                        RUNTIME_data_flush( sequence, T2(M+RD, k) );
                    }
                }

                /* Restore the original location of the tiles */
                for (n = 0; n < B->nt; n++) {
                    RUNTIME_data_migrate( sequence, B(k, n),
                                          B->get_rankof( B, k, n ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        } else {
            /*
             *  ChamLeft / ChamNoTrans
             */
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                lastRD = 0;
                for (RD = BS; RD < A->mt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldbM   = BLKLDD(B, M   );
                        ldbMRD = BLKLDD(B, M+RD);
                        ldaMRD = BLKLDD(A, M+RD);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                            RUNTIME_data_migrate( sequence, B(M, n),
                                                  B->get_rankof( B, M+RD, n ) );
                            RUNTIME_data_migrate( sequence, B(M+RD, n),
                                                  B->get_rankof( B, M+RD, n ) );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempMRDm, tempnn, tempkn, tempMRDm, ib, T->nb,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb,
                                B (M,    n), ldbM,
                                B (M+RD, n), ldbMRD);
                        }
                        RUNTIME_data_flush( sequence, A (M+RD, k) );
                        RUNTIME_data_flush( sequence, T2(M+RD, k) );
                    }
                }
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = chameleon_min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    lddM = BLKLDD(D, M);
                    ldbM = BLKLDD(B, M);
                    for (m = chameleon_min(M+BS, A->mt)-1; m > M; m--) {
                        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                        ldbm = BLKLDD(B, m);
                        ldam = BLKLDD(A, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                            RUNTIME_data_migrate( sequence, B(M, n),
                                                  B->get_rankof( B, m, n ) );
                            RUNTIME_data_migrate( sequence, B(m, n),
                                                  B->get_rankof( B, m, n ) );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(m, k), ldam,
                                T(m, k), T->mb,
                                B(M, n), ldbM,
                                B(m, n), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(m, k) );
                        RUNTIME_data_flush( sequence, T(m, k) );
                    }
                    if ( genD ) {
                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempMm, tempkmin, A->nb,
                            A(M, k), ldaM,
                            D(M, k), lddM );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempMm, tempkmin,
                            0., 1.,
                            D(M, k), lddM );
#endif
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(M, n),
                                              B->get_rankof( B, M, n ) );

                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            tempMm, tempnn, tempkmin, ib, T->nb,
                            D(M, k), lddM,
                            T(M, k), T->mb,
                            B(M, n), ldbM);
                    }
                    RUNTIME_data_flush( sequence, D(M, k) );
                    RUNTIME_data_flush( sequence, T(M, k) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }
    else {
        if (trans == ChamConjTrans) {
            /*
             *  ChamRight / ChamConjTrans
             */
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                lastRD = 0;
                for (RD = BS; RD < A->mt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldaMRD = BLKLDD(A, M+RD);
                        for (m = 0; m < B->mt; m++) {
                            ldbm   = BLKLDD(B, m);
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;

                            RUNTIME_data_migrate( sequence, B(m, M),
                                                  B->get_rankof( B, m, M+RD ) );
                            RUNTIME_data_migrate( sequence, B(m, M+RD),
                                                  B->get_rankof( B, m, M+RD ) );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempMRDm, tempkn, tempmm, ib, T->nb,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb,
                                B (m, M), ldbm,
                                B (m, M+RD), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A (M+RD, k) );
                        RUNTIME_data_flush( sequence, T2(M+RD, k) );
                    }
                }
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = chameleon_min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    lddM = BLKLDD(D, M);
                    for (n = chameleon_min(M+BS, A->mt)-1; n > M; n--) {
                        ldan = BLKLDD(A, n);
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        for (m = 0; m < B->mt; m++) {
                            ldbm = BLKLDD(B, m);
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;

                            RUNTIME_data_migrate( sequence, B(m, M),
                                                  B->get_rankof( B, m, n ) );
                            RUNTIME_data_migrate( sequence, B(m, m),
                                                  B->get_rankof( B, m, n ) );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(n, k), ldan,
                                T(n, k), T->mb,
                                B(m, M), ldbm,
                                B(m, n), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(n, k) );
                        RUNTIME_data_flush( sequence, T(n, k) );
                    }
                    if ( genD ) {
                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempMm, tempkmin, A->nb,
                            A(M, k), ldaM,
                            D(M, k), lddM );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempMm, tempkmin,
                            0., 1.,
                            D(M, k), lddM );
#endif
                    }
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;

                        RUNTIME_data_migrate( sequence, B(m, M),
                                              B->get_rankof( B, m, M ) );

                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempMm, tempkmin, ib, T->nb,
                            D(M, k), lddM,
                            T(M, k), T->mb,
                            B(m, M), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(M, k) );
                    RUNTIME_data_flush( sequence, T(M, k) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        } else {
            /*
             *  ChamRight / ChamNoTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = chameleon_min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    lddM = BLKLDD(D, M);
                    if ( genD ) {
                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempMm, tempkmin, A->nb,
                            A(M, k), ldaM,
                            D(M, k), lddM );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempMm, tempkmin,
                            0., 1.,
                            D(M, k), lddM );
#endif
                    }
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempMm, tempkmin, ib, T->nb,
                            D(M, k), lddM,
                            T(M, k), T->mb,
                            B(m, M), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(M, k) );
                    RUNTIME_data_flush( sequence, T(M, k) );

                    for (n = M+1; n < chameleon_min(M+BS,  A->mt); n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        ldan = BLKLDD(A, n);
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);

                            RUNTIME_data_migrate( sequence, B(m, M),
                                                  B->get_rankof( B, m, n ) );
                            RUNTIME_data_migrate( sequence, B(m, n),
                                                  B->get_rankof( B, m, n ) );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(n, k), ldan,
                                T(n, k), T->mb,
                                B(m, M), ldbm,
                                B(m, n), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A(n, k) );
                        RUNTIME_data_flush( sequence, T(n, k) );
                    }
                }
                for (RD = BS; RD < A->mt-k; RD *= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldaMRD = BLKLDD(A, M+RD);
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm   = BLKLDD(B, m);

                            RUNTIME_data_migrate( sequence, B(m, M),
                                                  B->get_rankof( B, m, M+RD ) );
                            RUNTIME_data_migrate( sequence, B(m, M+RD),
                                                  B->get_rankof( B, m, M+RD ) );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempMRDm, tempkn, tempmm, ib, T->nb,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb,
                                B (m, M   ), ldbm,
                                B (m, M+RD), ldbm);
                        }
                        RUNTIME_data_flush( sequence, A (M+RD, k) );
                        RUNTIME_data_flush( sequence, T2(M+RD, k) );
                    }
                }

                /* Restore the original location of the tiles */
                for (m = 0; m < B->mt; m++) {
                    RUNTIME_data_migrate( sequence, B(m, k),
                                          B->get_rankof( B, m, k ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
