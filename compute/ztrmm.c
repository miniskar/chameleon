/**
 *
 * @file ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm wrappers
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
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_ztrmm - Computes B = alpha*op( A )*B or B = alpha*B*op( A ).
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = ChamLeft:  A*X = B
 *          = ChamRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   A is transposed;
 *          = ChamTrans:     A is not transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = ChamNonUnit: A is non unit;
 *          = ChamUnit:    A us unit.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = ChamUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = ChamLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = ChamUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztrmm_Tile
 * @sa CHAMELEON_ztrmm_Tile_Async
 * @sa CHAMELEON_ctrmm
 * @sa CHAMELEON_dtrmm
 * @sa CHAMELEON_strmm
 *
 */
int CHAMELEON_ztrmm( cham_side_t side, cham_uplo_t uplo,
                 cham_trans_t transA, cham_diag_t diag,
                 int N, int NRHS, CHAMELEON_Complex64_t alpha,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB )
{
    int NB, NA;
    int status;
    CHAM_context_t *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_ztrmm", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if (side == ChamLeft) {
        NA = N;
    } else {
        NA = NRHS;
    }

    /* Check input arguments */
    if (side != ChamLeft && side != ChamRight) {
        morse_error("CHAMELEON_ztrmm", "illegal value of side");
        return -1;
    }
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        morse_error("CHAMELEON_ztrmm", "illegal value of uplo");
        return -2;
    }
    if (transA != ChamConjTrans &&
        transA != ChamNoTrans   &&
        transA != ChamTrans )
    {
        morse_error("CHAMELEON_ztrmm", "illegal value of transA");
        return -3;
    }
    if ((diag != ChamUnit) && (diag != ChamNonUnit)) {
        morse_error("CHAMELEON_ztrmm", "illegal value of diag");
        return -4;
    }
    if (N < 0) {
        morse_error("CHAMELEON_ztrmm", "illegal value of N");
        return -5;
    }
    if (NRHS < 0) {
        morse_error("CHAMELEON_ztrmm", "illegal value of NRHS");
        return -6;
    }
    if (LDA < chameleon_max(1, NA)) {
        morse_error("CHAMELEON_ztrmm", "illegal value of LDA");
        return -8;
    }
    if (LDB < chameleon_max(1, N)) {
        morse_error("CHAMELEON_ztrmm", "illegal value of LDB");
        return -10;
    }
    /* Quick return */
    if (chameleon_min(N, NRHS) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(CHAMELEON_FUNC_ZPOSV, N, N, NRHS);
    if (status != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_ztrmm", "morse_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = CHAMELEON_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, ChamDescInput, uplo,
                     A, NB, NB, LDA, NA, NA, NA, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_ztrmm_Tile_Async(  side, uplo, transA, diag, alpha, &descAt, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     ChamDescInput, uplo, sequence, &request );
    morse_ztile2lap( morse, &descBl, &descBt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descBl, &descBt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}
/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_ztrmm_Tile - Computes triangular solve.
 *  Tile equivalent of CHAMELEON_ztrmm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = ChamLeft:  A*X = B
 *          = ChamRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   A is transposed;
 *          = ChamTrans:     A is not transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = ChamNonUnit: A is non unit;
 *          = ChamUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = ChamUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = ChamLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = ChamUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztrmm
 * @sa CHAMELEON_ztrmm_Tile_Async
 * @sa CHAMELEON_ctrmm_Tile
 * @sa CHAMELEON_dtrmm_Tile
 * @sa CHAMELEON_strmm_Tile
 *
 */
int CHAMELEON_ztrmm_Tile( cham_side_t side, cham_uplo_t uplo,
                      cham_trans_t transA, cham_diag_t diag,
                      CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_ztrmm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    CHAMELEON_ztrmm_Tile_Async(side, uplo, transA, diag, alpha, A, B, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

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
 *  CHAMELEON_ztrmm_Tile_Async - Performs triangular matrix multiplication.
 *  Non-blocking equivalent of CHAMELEON_ztrmm_Tile().
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
 * @sa CHAMELEON_ztrmm
 * @sa CHAMELEON_ztrmm_Tile
 * @sa CHAMELEON_ctrmm_Tile_Async
 * @sa CHAMELEON_dtrmm_Tile_Async
 * @sa CHAMELEON_strmm_Tile_Async
 *
 */
int CHAMELEON_ztrmm_Tile_Async( cham_side_t side, cham_uplo_t uplo,
                            cham_trans_t transA, cham_diag_t diag,
                            CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_ztrmm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("CHAMELEON_ztrmm_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("CHAMELEON_ztrmm_Tile", "NULL request");
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
        morse_error("CHAMELEON_ztrmm_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_ztrmm_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("CHAMELEON_ztrmm_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (side != ChamLeft && side != ChamRight) {
        morse_error("CHAMELEON_ztrmm_Tile", "illegal value of side");
        return morse_request_fail(sequence, request, -1);
    }
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        morse_error("CHAMELEON_ztrmm_Tile", "illegal value of uplo");
        return morse_request_fail(sequence, request, -2);
    }
    if ((transA < ChamNoTrans) || (transA > ChamConjTrans)) {
        morse_error("CHAMELEON_ztrmm_Tile", "illegal value of transA");
        return morse_request_fail(sequence, request, -3);
    }
    if ((diag != ChamUnit) && (diag != ChamNonUnit)) {
        morse_error("CHAMELEON_ztrmm_Tile", "illegal value of diag");
        return morse_request_fail(sequence, request, -4);
    }

    /* Quick return */
    morse_pztrmm( side, uplo, transA, diag, alpha, A, B, sequence, request );

    return CHAMELEON_SUCCESS;
}
