/**
 *
 * @file pzsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm parallel algorithm
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2014-11-16
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
#define B(m, n) B,  m,  n
#define C(m, n) C,  m,  n
#define WA(m, n) &WA,  m,  n
#define WB(m, n) &WB,  m,  n

/**
 *  Parallel tile symmetric matrix-matrix multiplication.
 *  SUMMA algorithm for 2D block-cyclic data distribution.
 */
static inline void
chameleon_pzsymm_summa( CHAM_context_t *chamctxt, cham_side_t side, cham_uplo_t uplo,
                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                        CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                        RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;
    cham_trans_t transA, transB;
    int Am, An, m, n, k, p, q, KT, K, lp, lq;
    int ldam, ldbk, ldbm, ldcm;
    int tempmm, tempnn, tempkk;
    int lookahead, myp, myq;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;
    CHAM_desc_t WA, WB;

    lookahead = chamctxt->lookahead;
    chameleon_desc_init( &WA, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                         C->mt * C->mb, C->nb * C->q * lookahead, 0, 0,
                         C->mt * C->mb, C->nb * C->q * lookahead, C->p, C->q,
                         NULL, NULL, NULL );
    chameleon_desc_init( &WB, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                         C->mb * C->p * lookahead, C->nt * C->nb, 0, 0,
                         C->mb * C->p * lookahead, C->nt * C->nb, C->p, C->q,
                         NULL, NULL, NULL );

    KT  = side == ChamLeft ? A->nt : A->mt;
    K   = side == ChamLeft ? A->n  : A->m;
    myp = C->myrank / C->q;
    myq = C->myrank % C->q;

    for (k = 0; k < KT; k++ ) {
        lp = (k % lookahead) * C->p;
        lq = (k % lookahead) * C->q;
        tempkk = k == KT - 1 ? K - k * A->nb : A->nb;
        zbeta = k == 0 ? beta : zone;
        ldbk = BLKLDD(B, k);

        /* Transfert ownership of the k column of A or B */
        for (m = 0; m < C->mt; m ++ ) {
            tempmm = m == C->mt-1 ? C->m - m * C->mb : C->mb;

            if ( side == ChamLeft ) {
                if ( (( uplo == ChamUpper ) && ( m > k )) ||
                     (( uplo == ChamLower ) && ( m < k )) ) {
                    Am = k;
                    An = m;
                }
                else {
                    Am = m;
                    An = k;
                }
                ldam = BLKLDD(A, Am);

                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempmm, tempkk, C->mb,
                    A(  Am, An ),             ldam,
                    WA( m, (k % C->q) + lq ), WA.mb );

                RUNTIME_data_flush( sequence, A( Am, An ) );

                for ( q=1; q < C->q; q++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempmm, tempkk, C->mb,
                        WA( m, ((k+q-1) % C->q) + lq ), WA.mb,
                        WA( m, ((k+q)   % C->q) + lq ), WA.mb );
                }
            }
            else {
                ldbm = BLKLDD(B, m);

                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempmm, tempkk, C->mb,
                    B(  m,  k ),              ldbm,
                    WA( m, (k % C->q) + lq ), WA.mb );

                RUNTIME_data_flush( sequence, B( m, k ) );

                for ( q=1; q < C->q; q++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempmm, tempkk, C->mb,
                        WA( m, ((k+q-1) % C->q) + lq ), WA.mb,
                        WA( m, ((k+q)   % C->q) + lq ), WA.mb );
                }
            }
        }

        /* Transfert ownership of the k row of B, or A */
        for (n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

            if ( side == ChamRight ) {
                if ( (( uplo == ChamUpper ) && ( n < k )) ||
                     (( uplo == ChamLower ) && ( n > k )) ) {
                    Am = n;
                    An = k;
                }
                else {
                    Am = k;
                    An = n;
                }
                ldam = BLKLDD(A, Am);

                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempkk, tempnn, C->mb,
                    A(   Am,            An ), ldam,
                    WB( (k % C->p) + lp, n ), WB.mb );

                RUNTIME_data_flush( sequence, A( Am, An ) );

                for ( p=1; p < C->p; p++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempkk, tempnn, C->mb,
                        WB( ((k+p-1) % C->p) + lp, n ), WB.mb,
                        WB( ((k+p)   % C->p) + lp, n ), WB.mb );
                }
            }
            else {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempkk, tempnn, C->mb,
                    B(   k,              n ), ldbk,
                    WB( (k % C->p) + lp, n ), WB.mb );

                RUNTIME_data_flush( sequence, B( k, n ) );

                for ( p=1; p < C->p; p++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempkk, tempnn, C->mb,
                        WB( ((k+p-1) % C->p) + lp, n ), WB.mb,
                        WB( ((k+p)   % C->p) + lp, n ), WB.mb );
                }
            }
        }

        /*
         *  ChamLeft / ChamLower
         */
        for (m = myp; m < C->mt; m+=C->p) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
            ldcm = BLKLDD(C, m);

            for (n = myq; n < C->nt; n+=C->q) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                if (side == ChamLeft) {
                    transB = ChamNoTrans;
                    if ( (( uplo == ChamUpper ) && ( m > k )) ||
                         (( uplo == ChamLower ) && ( m < k )) ) {
                        transA = ChamTrans;
                    }
                    else {
                        transA = ChamNoTrans;
                    }
                }
                else {
                    transA = ChamNoTrans;
                    if ( (( uplo == ChamUpper ) && ( n < k )) ||
                         (( uplo == ChamLower ) && ( n > k )) ) {
                        transB = ChamTrans;
                    }
                    else {
                        transB = ChamNoTrans;
                    }
                }

                if ( k == m ) {
                    INSERT_TASK_zsymm(
                        options, side, uplo,
                        tempmm, tempnn, A->mb,
                        alpha, WA( m,        myq + lq ), WA.mb,  /* lda * Z */
                               WB( myp + lp, n        ), WB.mb,  /* ldb * Y */
                        zbeta, C(  m,        n        ), ldcm ); /* ldc * Y */
                }
                else {
                    INSERT_TASK_zgemm(
                        options, transA, transB,
                        tempmm, tempnn, tempkk, A->mb,
                        alpha, WA( m,        myq + lq ), WA.mb,  /* lda * Z */
                               WB( myp + lp, n        ), WB.mb,  /* ldb * Y */
                        zbeta, C(  m,        n        ), ldcm ); /* ldc * Y */
                }
            }
        }
    }

    RUNTIME_desc_flush( &WA, sequence );
    RUNTIME_desc_flush( &WB, sequence );
    RUNTIME_desc_flush(  C,  sequence );
    chameleon_sequence_wait( chamctxt, sequence );
    chameleon_desc_destroy( &WA );
    chameleon_desc_destroy( &WB );
}

/**
 *  Parallel tile symmetric matrix-matrix multiplication.
 *  Generic algorithm for any data distribution.
 */
static inline void
chameleon_pzsymm_generic( CHAM_context_t *chamctxt, cham_side_t side, cham_uplo_t uplo,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                          CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                          RUNTIME_option_t *options )
{
    int k, m, n;
    int ldam, ldan, ldak, ldbk, ldbm, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    for(m = 0; m < C->mt; m++) {
        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
        ldcm = BLKLDD(C, m);
        for(n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
            /*
             *  ChamLeft / ChamLower
             */
            if (side == ChamLeft) {
                ldam = BLKLDD(A, m);
                if (uplo == ChamLower) {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = k == C->mt-1 ? C->m-k*C->mb : C->mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(m, k), ldam,  /* lda * K */
                                       B(k, n), ldbk,  /* ldb * Y */
                                zbeta, C(m, n), ldcm); /* ldc * Y */
                        }
                        else {
                            if (k == m) {
                                INSERT_TASK_zsymm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * X */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(k, m), ldak,  /* ldak * X */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  ChamLeft / ChamUpper
                 */
                else {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = k == C->mt-1 ? C->m-k*C->mb : C->mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m), ldak,  /* ldak * X */
                                       B(k, n), ldbk,  /* ldb  * Y */
                                zbeta, C(m, n), ldcm); /* ldc  * Y */
                        }
                        else {
                            if (k == m) {
                                INSERT_TASK_zsymm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * K */
                                           B(k, n), ldbk,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamNoTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(m, k), ldam,  /* lda * K */
                                           B(k, n), ldbk,  /* ldb * Y */
                                    zbeta, C(m, n), ldcm); /* ldc * Y */
                            }
                        }
                    }
                }
            }
            /*
             *  ChamRight / ChamLower
             */
            else {
                ldan = BLKLDD(A, n);
                ldbm = BLKLDD(B, m);
                if (uplo == ChamLower) {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamNoTrans, ChamTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldb * K */
                                       A(n, k), ldan,  /* lda * K */
                                zbeta, C(m, n), ldcm); /* ldc * Y */
                        }
                        else {
                            if (k == n) {
                                INSERT_TASK_zsymm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * Y */
                                           B(m, k), ldbm,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamNoTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k), ldbm,  /* ldb  * K */
                                           A(k, n), ldak,  /* ldak * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  ChamRight / ChamUpper
                 */
                else {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldbm,  /* ldb  * K */
                                       A(k, n), ldak,  /* ldak * Y */
                                zbeta, C(m, n), ldcm); /* ldc  * Y */
                        }
                        else {
                            if (k == n) {
                                INSERT_TASK_zsymm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k), ldak,  /* ldak * Y */
                                           B(m, k), ldbm,  /* ldb  * Y */
                                    zbeta, C(m, n), ldcm); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamNoTrans, ChamTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k), ldbm,  /* ldb * K */
                                           A(n, k), ldan,  /* lda * K */
                                    zbeta, C(m, n), ldcm); /* ldc * Y */
                            }
                        }
                    }
                }
            }
        }
    }
    (void)chamctxt;
}

/**
 *  Parallel tile symmetric matrix-matrix multiplication. wrapper.
 */
void
chameleon_pzsymm( cham_side_t side, cham_uplo_t uplo,
                  CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                  CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    if ( ((C->p > 1) || (C->q > 1)) &&
         (C->get_rankof == chameleon_getrankof_2d) &&
         (chamctxt->generic_enabled != CHAMELEON_TRUE) )
    {
        chameleon_pzsymm_summa(   chamctxt, side, uplo, alpha, A, B, beta, C, &options );
    }
    else {
        chameleon_pzsymm_generic( chamctxt, side, uplo, alpha, A, B, beta, C, &options );
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
