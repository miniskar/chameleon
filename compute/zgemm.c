/**
 *
 * @file zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm wrappers
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

/**
 *
 * @defgroup CHAMELEON_Complex64_t
 * @brief Linear algebra routines exposed to users. LAPACK matrix data storage
 *
 */

/**
 *
 * @defgroup CHAMELEON_Complex64_t_Tile
 * @brief Linear algebra routines exposed to users. Tile matrix data storage
 *
 */

/**
 *
 * @defgroup CHAMELEON_Complex64_t_Tile_Async
 * @brief Linear algebra routines exposed to users. Tile matrix data storage,
 *  asynchronous interface.
 *
 */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgemm - Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans:     A is transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   B is not transposed;
 *          = ChamTrans:     B is transposed;
 *          = ChamConjTrans: B is conjugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix C. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( B ) and of the matrix C. N >= 0.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ) and the number of rows of
 *          the matrix op( B ). K >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = ChamNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = ChamNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgemm_Tile
 * @sa CHAMELEON_cgemm
 * @sa CHAMELEON_dgemm
 * @sa CHAMELEON_sgemm
 *
 */
int CHAMELEON_zgemm( cham_trans_t transA, cham_trans_t transB, int M, int N, int K,
                 CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB,
                 CHAMELEON_Complex64_t beta,  CHAMELEON_Complex64_t *C, int LDC )
{
    int NB;
    int Am, An, Bm, Bn;
    int status;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;
    CHAM_desc_t descCl, descCt;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((transA < ChamNoTrans) || (transA > ChamConjTrans)) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of transA");
        return -1;
    }
    if ((transB < ChamNoTrans) || (transB > ChamConjTrans)) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of transB");
        return -2;
    }
    if ( transA == ChamNoTrans ) {
        Am = M; An = K;
    } else {
        Am = K; An = M;
    }
    if ( transB == ChamNoTrans ) {
        Bm = K; Bn = N;
    } else {
        Bm = N; Bn = K;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of N");
        return -4;
    }
    if (K < 0) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of N");
        return -5;
    }
    if (LDA < chameleon_max(1, Am)) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of LDA");
        return -8;
    }
    if (LDB < chameleon_max(1, Bm)) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of LDB");
        return -10;
    }
    if (LDC < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgemm", "illegal value of LDC");
        return -13;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (CHAMELEON_Complex64_t)0.0 || K == 0) && beta == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgemm", "chameleon_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, &descAl, &descAt, ChamDescInput, ChamUpperLower,
                     A, NB, NB, LDA, An, Am, An, sequence, &request );
    chameleon_zlap2tile( chamctxt, &descBl, &descBt, ChamDescInput, ChamUpperLower,
                     B, NB, NB, LDB, Bn, Bm, Bn, sequence, &request );
    chameleon_zlap2tile( chamctxt, &descCl, &descCt, ChamDescInout, ChamUpperLower,
                     C, NB, NB, LDC, N, M,  N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgemm_Tile_Async( transA, transB, alpha, &descAt, &descBt, beta, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descCl, &descCt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descBl, &descBt );
    chameleon_ztile2lap_cleanup( chamctxt, &descCl, &descCt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgemm_Tile - Performs matrix multiplication.
 *  Tile equivalent of CHAMELEON_zgemm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans:     A is transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   B is not transposed;
 *          = ChamTrans:     B is transposed;
 *          = ChamConjTrans: B is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = ChamNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = ChamNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgemm
 * @sa CHAMELEON_zgemm_Tile_Async
 * @sa CHAMELEON_cgemm_Tile
 * @sa CHAMELEON_dgemm_Tile
 * @sa CHAMELEON_sgemm_Tile
 *
 */
int CHAMELEON_zgemm_Tile( cham_trans_t transA, cham_trans_t transB,
                      CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                      CHAMELEON_Complex64_t beta,  CHAM_desc_t *C )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgemm_Tile_Async( transA, transB, alpha, A, B, beta, C, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zgemm_Tile_Async - Performs matrix multiplication.
 *  Non-blocking equivalent of CHAMELEON_zgemm_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgemm
 * @sa CHAMELEON_zgemm_Tile
 * @sa CHAMELEON_cgemm_Tile_Async
 * @sa CHAMELEON_dgemm_Tile_Async
 * @sa CHAMELEON_sgemm_Tile_Async
 *
 */
int CHAMELEON_zgemm_Tile_Async( cham_trans_t transA, cham_trans_t transB,
                            CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                            CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    int M, N, K;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(C) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((transA < ChamNoTrans) || (transA > ChamConjTrans)) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "illegal value of transA");
        return chameleon_request_fail(sequence, request, -1);
    }
    if ((transB < ChamNoTrans) || (transB > ChamConjTrans)) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "illegal value of transB");
        return chameleon_request_fail(sequence, request, -2);
    }

    if ( transA == ChamNoTrans ) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
        Anb = A->nb;
        Ai  = A->i;
        Aj  = A->j;
    } else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
        Anb = A->mb;
        Ai  = A->j;
        Aj  = A->i;
    }

    if ( transB == ChamNoTrans ) {
        Bm  = B->m;
        Bn  = B->n;
        Bmb = B->mb;
        Bnb = B->nb;
        Bi  = B->i;
        Bj  = B->j;
    } else {
        Bm  = B->n;
        Bn  = B->m;
        Bmb = B->nb;
        Bnb = B->mb;
        Bi  = B->j;
        Bj  = B->i;
    }

    if ( (Amb != C->mb) || (Anb != Bmb) || (Bnb != C->nb) ) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "tile sizes have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (Am != C->m) || (An != Bm) || (Bn != C->n) ) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "sizes of matrices have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (Ai != C->i) || (Aj != Bi) || (Bj != C->j) ) {
        chameleon_error("CHAMELEON_zgemm_Tile_Async", "start indexes have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    M = C->m;
    N = C->n;
    K = An;

    /* Quick return */
    if ( (M == 0) || (N == 0) ||
         (((alpha == (CHAMELEON_Complex64_t)0.0) || (K == 0)) && (beta == (CHAMELEON_Complex64_t)1.0)) )
    {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzgemm( transA, transB, alpha, A, B, beta, C, sequence, request );

    return CHAMELEON_SUCCESS;
}
