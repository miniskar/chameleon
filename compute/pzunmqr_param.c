/**
 *
 * @file pzunmqr_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr_param parallel algorithm
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2017-05-05
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define T(m,n) T,  m,  n
#define D(m,n) D,  m,  n

/**
 *  Parallel application of Q using tile V - QR factorization - dynamic scheduling
 */
void chameleon_pzunmqr_param( int genD, const libhqr_tree_t *qrtree,
                              cham_side_t side, cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *B,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int ldam, ldan, ldbm, ldbp, lddn, lddm;
    int tempnn, tempkmin, tempmm, tempkn;
    int ib, K, L;
    int node, nbtiles, *tiles;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    K = chameleon_min(A->mt, A->nt);

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmqr  = A->nb * ib
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr  =      A->nb * ib
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of tiles */
    tiles = (int*)calloc( qrtree->mt, sizeof(int) );

    if (side == ChamLeft ) {
        if (trans == ChamConjTrans) {
            /*
             *  ChamLeft / ChamConjTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkn);
                    ldam = BLKLDD(A, m);
                    lddm = BLKLDD(D, m);
                    ldbm = BLKLDD(B, m);

                    if ( genD ) {
                        int tempDmm = m == D->mt-1 ? D->m-m*D->mb : D->mb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDmm, tempkmin, A->nb,
                            A(m, k), ldam,
                            D(m, k), lddm );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDmm, tempkmin,
                            0., 1.,
                            D(m, k), lddm );
#endif
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(m, k), lddm,
                            T(m, k), T->mb,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = 0; i < nbtiles; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);

                    if( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempmm;
                        T = TT;
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        node = B->get_rankof( B, m, n );
                        RUNTIME_data_migrate( sequence, B(p, n), node );
                        RUNTIME_data_migrate( sequence, B(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, L, ib, T->nb,
                            A(m, k), ldam,
                            T(m, k), T->mb,
                            B(p, n), ldbp,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                /* Restore the original location of the tiles */
                for (n = 0; n < B->nt; n++) {
                    RUNTIME_data_migrate( sequence, B(k, n),
                                          B->get_rankof( B, k, n ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamLeft / ChamNoTrans
         */
        else {
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = nbtiles-1; i >=0; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);

                    if( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempmm;
                        T = TT;
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        node = B->get_rankof( B, m, n );
                        RUNTIME_data_migrate( sequence, B(p, n), node );
                        RUNTIME_data_migrate( sequence, B(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, L, ib, T->nb,
                            A(m, k), ldam,
                            T(m, k), T->mb,
                            B(p, n), ldbp,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkn);
                    ldam = BLKLDD(A, m);
                    lddm = BLKLDD(D, m);
                    ldbm = BLKLDD(B, m);

                    if ( genD ) {
                        int tempDmm = m == D->mt-1 ? D->m-m*D->mb : D->mb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDmm, tempkmin, A->nb,
                            A(m, k), ldam,
                            D(m, k), lddm );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDmm, tempkmin,
                            0., 1.,
                            D(m, k), lddm );
#endif
                    }
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(m, k), lddm,
                            T(m, k), T->mb,
                            B(m, n), ldbm);
                    }

                    RUNTIME_data_flush( sequence, D(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }
    /*
     *  ChamRight / ChamConjTrans
     */
    else {
        if (trans == ChamConjTrans) {
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n - k*A->nb : A->nb;

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = nbtiles-1; i >= 0; i--) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);

                    if( qrtree->gettype(qrtree, k, n) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->mb;
                        T = TT;
                    }

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        node = B->get_rankof( B, m, n );
                        RUNTIME_data_migrate( sequence, B(m, p), node );
                        RUNTIME_data_migrate( sequence, B(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, chameleon_min( L, tempmm ), ib, T->nb,
                            A(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, p), ldbm,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    n = qrtree->getm(qrtree, k, i);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    tempkmin = chameleon_min(tempnn, tempkn);
                    ldan = BLKLDD(A, n);
                    lddn = BLKLDD(D, n);

                    if ( genD ) {
                        int tempDnn = n == D->nt-1 ? D->n-n*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDnn, tempkmin, A->nb,
                            A(n, k), ldan,
                            D(n, k), lddn );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDnn, tempkmin,
                            0., 1.,
                            D(n, k), lddn );
#endif
                    }
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;

                        RUNTIME_data_migrate( sequence, B(m, n),
                                              B->get_rankof( B, m, n ) );

                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(n, k), lddn,
                            T(n, k), T->mb,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamRight / ChamNoTrans
         */
        else {
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == B->nt-1 ? B->n-k*B->nb : B->nb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    n = qrtree->getm(qrtree, k, i);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    tempkmin = chameleon_min(tempnn, tempkn);
                    ldan = BLKLDD(A, n);
                    lddn = BLKLDD(D, n);

                    if ( genD ) {
                        int tempDnn = n == D->nt-1 ? D->n-n*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDnn, tempkmin, A->nb,
                            A(n, k), ldan,
                            D(n, k), lddn );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDnn, tempkmin,
                            0., 1.,
                            D(n, k), lddn );
#endif
                    }
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, T->nb,
                            D(n, k), lddn,
                            T(n, k), T->mb,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, D(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                /* Setting the order of tiles */
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = 0; i < nbtiles; i++) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);

                    if( qrtree->gettype(qrtree, k, n) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->mb;
                        T = TT;
                    }

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);

                        node = B->get_rankof( B, m, n );
                        RUNTIME_data_migrate( sequence, B(m, p), node );
                        RUNTIME_data_migrate( sequence, B(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkn, chameleon_min( L, tempmm ), ib, T->nb,
                            A(n, k), ldan,
                            T(n, k), T->mb,
                            B(m, p), ldbm,
                            B(m, n), ldbm);
                    }
                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
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

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
