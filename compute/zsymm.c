/**
 *
 * @file zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm wrappers
 *
 * @version 1.1.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2020-11-18
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zsymm - Performs one of the matrix-matrix operations
 *
 *     \f[ C = \alpha \times A \times B + \beta \times C \f]
 *
 *  or
 *
 *     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 *  where alpha and beta are scalars, A is a symmetric matrix and  B and
 *  C are m by n matrices.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the symmetric matrix A appears on the
 *          left or right in the operation as follows:
 *          = ChamLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = ChamRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = ChamLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = ChamUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 * @param[in] M
 *          Specifies the number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          Specifies the number of columns of the matrix C. N >= 0.
 *
 * @param[in] alpha
 *          Specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is M when side = ChamLeft,
 *          and is N otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,ka).
 *
 * @param[in] B
 *          B is a LDB-by-N matrix, where the leading M-by-N part of
 *          the array B must contain the matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 * @param[in] beta
 *          Specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zsymm_Tile
 * @sa CHAMELEON_csymm
 * @sa CHAMELEON_dsymm
 * @sa CHAMELEON_ssymm
 *
 */
int CHAMELEON_zsymm( cham_side_t side, cham_uplo_t uplo, int M, int N,
                 CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB,
                 CHAMELEON_Complex64_t beta,  CHAMELEON_Complex64_t *C, int LDC )
{
    int NB;
    int Am;
    int status;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;
    CHAM_desc_t descCl, descCt;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsymm", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( (side != ChamLeft) && (side != ChamRight) ){
        chameleon_error("CHAMELEON_zsymm", "illegal value of side");
        return -1;
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        chameleon_error("CHAMELEON_zsymm", "illegal value of uplo");
        return -2;
    }
    Am = ( side == ChamLeft ) ? M : N;
    if (M < 0) {
        chameleon_error("CHAMELEON_zsymm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zsymm", "illegal value of N");
        return -4;
    }
    if (LDA < chameleon_max(1, Am)) {
        chameleon_error("CHAMELEON_zsymm", "illegal value of LDA");
        return -7;
    }
    if (LDB < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zsymm", "illegal value of LDB");
        return -9;
    }
    if (LDC < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zsymm", "illegal value of LDC");
        return -12;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (CHAMELEON_Complex64_t)0.0) && beta == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZSYMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zsymm", "chameleon_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, &descAl, &descAt, ChamDescInput, uplo,
                     A, NB, NB, LDA, Am, Am, Am, sequence, &request );
    chameleon_zlap2tile( chamctxt, &descBl, &descBt, ChamDescInput, ChamUpperLower,
                     B, NB, NB, LDB, N, M,  N, sequence, &request );
    chameleon_zlap2tile( chamctxt, &descCl, &descCt, ChamDescInout, ChamUpperLower,
                     C, NB, NB, LDC, N, M,  N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zsymm_Tile_Async(  side, uplo, alpha, &descAt, &descBt, beta, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, uplo, sequence, &request );
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
 *  CHAMELEON_zsymm_Tile - Performs symmetric matrix multiplication.
 *  Tile equivalent of CHAMELEON_zsymm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the symmetric matrix A appears on the
 *          left or right in the operation as follows:
 *          = ChamLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = ChamRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = ChamLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = ChamUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 * @param[in] alpha
 *          Specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is M when side = ChamLeft,
 *          and is N otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] B
 *          B is a LDB-by-N matrix, where the leading M-by-N part of
 *          the array B must contain the matrix B.
 *
 * @param[in] beta
 *          Specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N updated matrix.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zsymm
 * @sa CHAMELEON_zsymm_Tile_Async
 * @sa CHAMELEON_csymm_Tile
 * @sa CHAMELEON_dsymm_Tile
 * @sa CHAMELEON_ssymm_Tile
 *
 */
int CHAMELEON_zsymm_Tile( cham_side_t side, cham_uplo_t uplo,
                      CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                      CHAMELEON_Complex64_t beta,  CHAM_desc_t *C )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsymm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zsymm_Tile_Async(side, uplo, alpha, A, B, beta, C, sequence, &request );

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
 *  CHAMELEON_zsymm_Tile_Async - Performs symmetric matrix multiplication.
 *  Non-blocking equivalent of CHAMELEON_zsymm_Tile().
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
 * @sa CHAMELEON_zsymm
 * @sa CHAMELEON_zsymm_Tile
 * @sa CHAMELEON_csymm_Tile_Async
 * @sa CHAMELEON_dsymm_Tile_Async
 * @sa CHAMELEON_ssymm_Tile_Async
 *
 */
int CHAMELEON_zsymm_Tile_Async( cham_side_t side, cham_uplo_t uplo,
                            CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                            CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsymm_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zsymm_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zsymm_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(C) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (side != ChamLeft) && (side != ChamRight) ){
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "illegal value of side");
        return chameleon_request_fail(sequence, request, -1);
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -2);
    }

    /* Check matrices sizes */
    if ( (B->m != C->m) || (B->n != C->n) ) {
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "B and C must have the same size");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (A->m != A->n) ||
         ( (side == ChamLeft)  && (A->m != B->m ) ) ||
         ( (side == ChamRight) && (A->m != B->n ) ) ) {
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "Matrix A must be square of size M or N regarding side.");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check tiles sizes */
    if ( (B->mb != C->mb) || (B->nb != C->nb) ) {
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "B and C must have the same tile sizes");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (A->mb != A->nb) ||
         ( (side == ChamLeft)  && (A->mb != B->mb ) ) ||
         ( (side == ChamRight) && (A->mb != B->nb ) ) ) {
        chameleon_error("CHAMELEON_zsymm_Tile_Async", "Matrix A must be square with square tiles wich fits the reagding tile size of B and C");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check submatrix starting point */
    /* if ( (B->i != C->i) || (B->j != C->j) ) { */
    /*     chameleon_error("CHAMELEON_zsymm_Tile_Async", "B and C submatrices doesn't match"); */
    /*     return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE); */
    /* } */
    /* if ( (A->i != A->j) ||  */
    /*          ( (side == ChamLeft)  && (A->i != B->i ) ) ||  */
    /*          ( (side == ChamRight) && (A->i != B->j ) ) ) { */
    /*     chameleon_error("CHAMELEON_zsymm_Tile_Async", "Submatrix A must start on diagnonal and match submatrices B and C."); */
    /*     return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE); */
    /* } */
    if( (A->i != 0) || (A->j != 0) ||
        (B->i != 0) || (B->j != 0) ||
        (C->i != 0) || (C->j != 0) ) {
        chameleon_error("CHAMELEON_zhemm_Tile_Async", "Submatrices are not supported for now");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( (C->m == 0) || (C->n == 0) ||
         ( (alpha == (CHAMELEON_Complex64_t)0.0) && (beta == (CHAMELEON_Complex64_t)1.0) ) )
    {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzsymm( side, uplo, alpha, A, B, beta, C, sequence, request );

    return CHAMELEON_SUCCESS;
}
