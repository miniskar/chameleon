/**
 *
 * @file pzgeqrf_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf_param parallel algorithm
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2017-05-17
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>
#include "libhqr.h"

#define A(m,n)  A, (m), (n)
#define T(m,n)  T, (m), (n)
#define D(m,n)  D, (m), (n)

/**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 */
void chameleon_pzgeqrf_param( int genD, const libhqr_tree_t *qrtree, CHAM_desc_t *A,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int K, L, nbgeqrt;
    int ldap, ldam, lddm;
    int tempkmin, tempkn, tempnn, tempmm;
    int ib;
    int *tiles;

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
     * zgeqrt  = A->nb * (ib+1)
     * zunmqr  = A->nb * ib
     * ztpqrt  = A->nb * (ib+1)
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /*
     * zunmqr  = A->nb * ib
     * ztpmqrt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of temporary tiles array */
    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);

    /* The number of the factorization */
    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(chamctxt, k);
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

        /* The number of geqrt to apply */
        nbgeqrt = qrtree->getnbgeqrf(qrtree, k);
        for (i = 0; i < nbgeqrt; i++) {
            m = qrtree->getm(qrtree, k, i);
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            tempkmin = chameleon_min(tempmm, tempkn);
            ldam = BLKLDD(A, m);
            lddm = BLKLDD(D, m);

            T = TS;

            INSERT_TASK_zgeqrt(
                &options,
                tempmm, tempkn, ib, T->nb,
                A(m, k), ldam,
                T(m, k), T->mb);
            if ( genD ) {
                INSERT_TASK_zlacpy(
                    &options,
                    ChamLower, tempmm, tempkn, A->nb,
                    A(m, k), ldam,
                    D(m, k), lddm );
#if defined(CHAMELEON_USE_CUDA)
                INSERT_TASK_zlaset(
                    &options,
                    ChamUpper, tempmm, tempkn,
                    0., 1.,
                    D(m, k), lddm );
#endif
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                INSERT_TASK_zunmqr(
                    &options,
                    ChamLeft, ChamConjTrans,
                    tempmm, tempnn, tempkmin, ib, T->nb,
                    D(m, k), lddm,
                    T(m, k), T->mb,
                    A(m, n), ldam);
            }
            RUNTIME_data_flush( sequence, D(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        /* Setting the order of the tiles */
        libhqr_walk_stepk( qrtree, k, tiles + (k+1) );

        for (i = k+1; i < A->mt; i++) {
            m = tiles[i];
            p = qrtree->currpiv(qrtree, k, m);

            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldap = BLKLDD(A, p);
            ldam = BLKLDD(A, m);

            if (qrtree->gettype(qrtree, k, m) == 0) {
                /* TS kernel */
                T = TS;
                L = 0;
            }
            else {
                /* TT kernel */
                T = TT;
                L = tempmm;
            }

            RUNTIME_data_migrate( sequence, A(p, k),
                                  A->get_rankof( A, m, k ) );
            RUNTIME_data_migrate( sequence, A(m, k),
                                  A->get_rankof( A, m, k ) );

            INSERT_TASK_ztpqrt(
                &options,
                tempmm, tempkn, chameleon_min(L, tempkn), ib, T->nb,
                A(p, k), ldap,
                A(m, k), ldam,
                T(m, k), T->mb);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                RUNTIME_data_migrate( sequence, A(p, n),
                                      A->get_rankof( A, m, n ) );
                RUNTIME_data_migrate( sequence, A(m, n),
                                      A->get_rankof( A, m, n ) );

                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamConjTrans,
                    tempmm, tempnn, A->nb, L, ib, T->nb,
                    A(m, k), ldam,
                    T(m, k), T->mb,
                    A(p, n), ldap,
                    A(m, n), ldam);
            }
            RUNTIME_data_flush( sequence, A(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        /* Restore the original location of the tiles */
        for (n = k; n < A->nt; n++) {
            RUNTIME_data_migrate( sequence, A(k, n),
                                  A->get_rankof( A, k, n ) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
