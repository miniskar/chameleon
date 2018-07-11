/**
 *
 * @file zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zhemm - Performs one of the matrix-matrix operations
 *
 *     \f[ C = \alpha \times A \times B + \beta \times C \f]
 *
 *  or
 *
 *     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 *  where alpha and beta are scalars, A is an hermitian matrix and  B and
 *  C are m by n matrices.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the hermitian matrix A appears on the
 *          left or right in the operation as follows:
 *          = ChamLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = ChamRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the hermitian matrix A is to be referenced as follows:
 *          = ChamLower:     Only the lower triangular part of the
 *                             hermitian matrix A is to be referenced.
 *          = ChamUpper:     Only the upper triangular part of the
 *                             hermitian matrix A is to be referenced.
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
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhemm_Tile
 * @sa CHAMELEON_chemm
 * @sa CHAMELEON_dhemm
 * @sa CHAMELEON_shemm
 *
 */
int CHAMELEON_zhemm( cham_side_t side, cham_uplo_t uplo, int M, int N,
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
    CHAM_context_t *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zhemm", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( (side != ChamLeft) && (side != ChamRight) ){
        morse_error("CHAMELEON_zhemm", "illegal value of side");
        return -1;
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        morse_error("CHAMELEON_zhemm", "illegal value of uplo");
        return -2;
    }
    Am = ( side == ChamLeft ) ? M : N;
    if (M < 0) {
        morse_error("CHAMELEON_zhemm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        morse_error("CHAMELEON_zhemm", "illegal value of N");
        return -4;
    }
    if (LDA < chameleon_max(1, Am)) {
        morse_error("CHAMELEON_zhemm", "illegal value of LDA");
        return -7;
    }
    if (LDB < chameleon_max(1, M)) {
        morse_error("CHAMELEON_zhemm", "illegal value of LDB");
        return -9;
    }
    if (LDC < chameleon_max(1, M)) {
        morse_error("CHAMELEON_zhemm", "illegal value of LDC");
        return -12;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (CHAMELEON_Complex64_t)0.0) && beta == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(CHAMELEON_FUNC_ZHEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zhemm", "morse_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = CHAMELEON_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, ChamDescInput, uplo,
                     A, NB, NB, LDA, Am, Am, Am, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, ChamDescInput, ChamUpperLower,
                     B, NB, NB, LDB, N, M,  N, sequence, &request );
    morse_zlap2tile( morse, &descCl, &descCt, ChamDescInout, ChamUpperLower,
                     C, NB, NB, LDC, N, M,  N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zhemm_Tile_Async(  side, uplo, alpha, &descAt, &descBt, beta, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     ChamDescInput, uplo, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
    morse_ztile2lap( morse, &descCl, &descCt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descBl, &descBt );
    morse_ztile2lap_cleanup( morse, &descCl, &descCt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zhemm_Tile - Performs Hermitian matrix multiplication.
 *  Tile equivalent of CHAMELEON_zhemm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the hermitian matrix A appears on the
 *          left or right in the operation as follows:
 *          = ChamLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = ChamRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the hermitian matrix A is to be referenced as follows:
 *          = ChamLower:     Only the lower triangular part of the
 *                             hermitian matrix A is to be referenced.
 *          = ChamUpper:     Only the upper triangular part of the
 *                             hermitian matrix A is to be referenced.
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
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhemm
 * @sa CHAMELEON_zhemm_Tile_Async
 * @sa CHAMELEON_chemm_Tile
 * @sa CHAMELEON_dhemm_Tile
 * @sa CHAMELEON_shemm_Tile
 *
 */
int CHAMELEON_zhemm_Tile( cham_side_t side, cham_uplo_t uplo,
                      CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                      CHAMELEON_Complex64_t beta,  CHAM_desc_t *C )
{
    CHAM_context_t *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zhemm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    CHAMELEON_zhemm_Tile_Async(side, uplo, alpha, A, B, beta, C, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

    morse_sequence_wait( morse, sequence );
    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zhemm_Tile_Async - Performs Hermitian matrix multiplication.
 *  Non-blocking equivalent of CHAMELEON_zhemm_Tile().
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
 * @sa CHAMELEON_zhemm
 * @sa CHAMELEON_zhemm_Tile
 * @sa CHAMELEON_chemm_Tile_Async
 * @sa CHAMELEON_dhemm_Tile_Async
 * @sa CHAMELEON_shemm_Tile_Async
 *
 */
int CHAMELEON_zhemm_Tile_Async( cham_side_t side, cham_uplo_t uplo,
                            CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                            CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zhemm_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("CHAMELEON_zhemm_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("CHAMELEON_zhemm_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return morse_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "invalid first descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "invalid second descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(C) != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "invalid third descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (side != ChamLeft) && (side != ChamRight) ){
        morse_error("CHAMELEON_zhemm_Tile_Async", "illegal value of side");
        return morse_request_fail(sequence, request, -1);
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "illegal value of uplo");
        return morse_request_fail(sequence, request, -2);
    }

    /* Check matrices sizes */
    if ( (B->m != C->m) || (B->n != C->n) ) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "B and C must have the same size");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (A->m != A->n) ||
         ( (side == ChamLeft)  && (A->m != B->m ) ) ||
         ( (side == ChamRight) && (A->m != B->n ) ) ) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "Matrix A must be square of size M or N regarding side.");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check tiles sizes */
    if ( (B->mb != C->mb) || (B->nb != C->nb) ) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "B and C must have the same tile sizes");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (A->mb != A->nb) ||
         ( (side == ChamLeft)  && (A->mb != B->mb ) ) ||
         ( (side == ChamRight) && (A->mb != B->nb ) ) ) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "Matrix A must be square with square tiles wich fits the reagding tile size of B and C");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check submatrix starting point */
    /* if ( (B->i != C->i) || (B->j != C->j) ) { */
    /*     morse_error("CHAMELEON_zhemm_Tile_Async", "B and C submatrices doesn't match"); */
    /*     return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE); */
    /* } */
    /* if ( (A->i != A->j) ||  */
    /*          ( (side == ChamLeft)  && (A->i != B->i ) ) ||  */
    /*          ( (side == ChamRight) && (A->i != B->j ) ) ) { */
    /*     morse_error("CHAMELEON_zhemm_Tile_Async", "Submatrix A must start on diagnonal and match submatrices B and C."); */
    /*     return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE); */
    /* } */
    if( (A->i != 0) || (A->j != 0) ||
        (B->i != 0) || (B->j != 0) ||
        (C->i != 0) || (C->j != 0) ) {
        morse_error("CHAMELEON_zhemm_Tile_Async", "Submatrices are not supported for now");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( (C->m == 0) || (C->n == 0) ||
         ( (alpha == (CHAMELEON_Complex64_t)0.0) && (beta == (CHAMELEON_Complex64_t)1.0) ) )
    {
        return CHAMELEON_SUCCESS;
    }

    morse_pzhemm( side, uplo, alpha, A, B, beta, C, sequence, request );

    return CHAMELEON_SUCCESS;
}
