/**
 *
 * @file pzunglq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq parallel algorithm
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2014-11-16
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define Q(m,n) Q,  m,  n
#define T(m,n) T,  m,  n
#define D(k)   D,  k,  k

/**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 */
void chameleon_pzunglq( int genD, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldak, ldqm, lddk;
    int tempnn, tempmm, tempkmin, tempkn;
    int tempAkm, tempAkn;
    int ib, minMT;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if (A->m > A->n) {
        minMT = A->nt;
    } else {
        minMT = A->mt;
    }

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmlq  = A->nb * ib
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq  =     A->nb * ib
     * ztpmlqt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = minMT-1; k >= 0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempAkm  = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempAkn  = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        tempkmin = chameleon_min( tempAkn, tempAkm );
        tempkn   = k == Q->nt-1 ? Q->n-k*Q->nb : Q->nb;
        ldak = BLKLDD(A, k);
        lddk = BLKLDD(D, k);

        for (n = Q->nt-1; n > k; n--) {
            tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;
            for (m = k; m < Q->mt; m++) {
                tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
                ldqm = BLKLDD(Q, m);

                RUNTIME_data_migrate( sequence, Q(m, k),
                                      Q->get_rankof( Q, m, n ) );

                /* TS kernel */
                INSERT_TASK_ztpmlqt(
                    &options,
                    ChamRight, ChamNoTrans,
                    tempmm, tempnn, tempAkm, 0, ib, T->nb,
                    A(k, n), ldak,
                    T(k, n), T->mb,
                    Q(m, k), ldqm,
                    Q(m, n), ldqm);
            }
            RUNTIME_data_flush( sequence, A(k, n) );
            RUNTIME_data_flush( sequence, T(k, n) );
        }

        if ( genD ) {
            int tempDkn = k == D->nt-1 ? D->n-k*D->nb : D->nb;
            INSERT_TASK_zlacpy(
                &options,
                ChamUpper, tempkmin, tempDkn, A->nb,
                A(k, k), ldak,
                D(k),    lddk );
#if defined(CHAMELEON_USE_CUDA)
            INSERT_TASK_zlaset(
                &options,
                ChamLower, tempkmin, tempDkn,
                0., 1.,
                D(k), lddk );
#endif
        }
        for (m = k; m < Q->mt; m++) {
            tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
            ldqm = BLKLDD(Q, m);

            /* Restore the original location of the tiles */
            RUNTIME_data_migrate( sequence, Q(m, k),
                                  Q->get_rankof( Q, m, k ) );

            INSERT_TASK_zunmlq(
                &options,
                ChamRight, ChamNoTrans,
                tempmm, tempkn, tempkmin, ib, T->nb,
                D(k),    lddk,
                T(k, k), T->mb,
                Q(m, k), ldqm);
        }
        RUNTIME_data_flush( sequence, D(k)    );
        RUNTIME_data_flush( sequence, T(k, k) );

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
