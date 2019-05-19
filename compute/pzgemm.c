/**
 *
 * @file pzgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm parallel algorithm
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
 *  Parallel tile matrix-matrix multiplication
 *  SUMMA algorithm for 2D block-cyclic data distribution.
 */
static inline void
chameleon_pzgemm_summa( CHAM_context_t *chamctxt, cham_trans_t transA, cham_trans_t transB,
                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                        CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                        RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;
    int m, n, k, p, q, KT, K, lp, lq;
    int ldam, ldak, ldbn, ldbk, ldcm;
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

    KT  = transA == ChamNoTrans ? A->nt : A->mt;
    K   = transA == ChamNoTrans ? A->n  : A->m;
    myp = A->myrank / A->q;
    myq = A->myrank % A->q;

    /*
     *  A: ChamNoTrans / B: ChamNoTrans
     */
    for (k = 0; k < KT; k++ ) {
        lp = (k % lookahead) * C->p;
        lq = (k % lookahead) * C->q;
        tempkk = k == KT - 1 ? K - k * A->nb : A->nb;
        zbeta = k == 0 ? beta : zone;
        ldak = BLKLDD(A, k);
        ldbk = BLKLDD(B, k);

        /* Transfert ownership of the k column of A */
        for (m = 0; m < C->mt; m ++ ) {
            tempmm = m == C->mt-1 ? C->m - m * C->mb : C->mb;
            ldam = BLKLDD(A, m);

            if ( transA == ChamNoTrans ) {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempmm, tempkk, C->mb,
                    A(  m,  k ),             ldam,
                    WA( m, (k % C->q) + lq ), WA.mb );

                RUNTIME_data_flush( sequence, A( m, k ) );

                for ( q=1; q < C->q; q++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempmm, tempkk, C->mb,
                        WA( m, ((k+q-1) % C->q) + lq ), WA.mb,
                        WA( m, ((k+q)   % C->q) + lq ), WA.mb );
                }
            }
            else {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempkk, tempmm, C->mb,
                    A(  k,  m ),              ldak,
                    WA( m, (k % C->q) + lq ), WA.mb );

                RUNTIME_data_flush( sequence, A( k, m ) );

                for ( q=1; q < C->q; q++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempkk, tempmm, C->mb,
                        WA( m, ((k+q-1) % C->q) + lq ), WA.mb,
                        WA( m, ((k+q)   % C->q) + lq ), WA.mb );
                }
            }
        }

        /* Transfert ownership of the k row of B */
        for (n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
            ldbn = BLKLDD(B, n);

            if ( transB == ChamNoTrans ) {
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
            else {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempnn, tempkk, C->mb,
                    B(   n,              k ), ldbn,
                    WB( (k % C->p) + lp, n ), WB.mb );

                RUNTIME_data_flush( sequence, B( n, k ) );

                for ( p=1; p < C->p; p++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempnn, tempkk, C->mb,
                        WB( ((k+p-1) % C->p) + lp, n ), WB.mb,
                        WB( ((k+p)   % C->p) + lp, n ), WB.mb );
                }
            }
        }

        for (m = myp; m < C->mt; m+=C->p) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
            ldcm = BLKLDD(C, m);

            for (n = myq; n < C->nt; n+=C->q) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                INSERT_TASK_zgemm(
                    options,
                    transA, transB,
                    tempmm, tempnn, tempkk, A->mb,
                    alpha, WA( m,        myq + lq ), WA.mb,  /* lda * Z */
                           WB( myp + lp, n        ), WB.mb,  /* ldb * Y */
                    zbeta, C(  m,        n        ), ldcm ); /* ldc * Y */
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
 *  Parallel tile matrix-matrix multiplication.
 *  Generic algorithm for any data distribution.
 */
static inline void
chameleon_pzgemm_generic( CHAM_context_t *chamctxt, cham_trans_t transA, cham_trans_t transB,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                          CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                          RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;

    int m, n, k;
    int ldam, ldak, ldbn, ldbk, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

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
                            options,
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
                            options,
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
                            options,
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
                            options,
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

    (void)chamctxt;
}

/**
 *  Parallel tile matrix-matrix multiplication. wrapper.
 */
void
chameleon_pzgemm( cham_trans_t transA, cham_trans_t transB,
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
        chameleon_pzgemm_summa(   chamctxt, transA, transB, alpha, A, B, beta, C, &options );
    }
    else {
        chameleon_pzgemm_generic( chamctxt, transA, transB, alpha, A, B, beta, C, &options );
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
