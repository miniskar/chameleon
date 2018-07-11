/**
 *
 * @file zgetrs_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrs_nopiv wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2014-11-08
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgetrs_nopiv - Solves a system of linear equations A * X = B, with a general N-by-N matrix A
 *  using the tile LU factorization computed by CHAMELEON_zgetrf_nopiv.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = ChamNoTrans:   A * X = B     (No transpose)
 *          = ChamTrans:     A**T * X = B  (Transpose)
 *          = ChamConjTrans: A**H * X = B  (Conjugate transpose)
 *          Currently only ChamNoTrans is supported.
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by CHAMELEON_zgetrf_nopiv.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \return <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrs_nopiv_Tile
 * @sa CHAMELEON_zgetrs_nopiv_Tile_Async
 * @sa CHAMELEON_cgetrs_nopiv
 * @sa CHAMELEON_dgetrs_nopiv
 * @sa CHAMELEON_sgetrs_nopiv
 * @sa CHAMELEON_zgetrf_nopiv
 *
 */
int CHAMELEON_zgetrs_nopiv( cham_trans_t trans, int N, int NRHS,
                        CHAMELEON_Complex64_t *A, int LDA,
                        CHAMELEON_Complex64_t *B, int LDB )
{
    int NB;
    int status;
    CHAM_context_t *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zgetrs_nopiv", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (trans != ChamNoTrans) {
        morse_error("CHAMELEON_zgetrs_nopiv", "only ChamNoTrans supported");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    if (N < 0) {
        morse_error("CHAMELEON_zgetrs_nopiv", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        morse_error("CHAMELEON_zgetrs_nopiv", "illegal value of NRHS");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("CHAMELEON_zgetrs_nopiv", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, N)) {
        morse_error("CHAMELEON_zgetrs_nopiv", "illegal value of LDB");
        return -9;
    }
    /* Quick return */
    if (chameleon_min(N, NRHS) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNBSIZE */
    status = morse_tune(CHAMELEON_FUNC_ZGESV, N, N, NRHS);
    if (status != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zgetrs_nopiv", "morse_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = CHAMELEON_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, ChamDescInput, ChamUpperLower,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    morse_zlap2tile( morse, &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgetrs_nopiv_Tile_Async( &descAt, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
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
 *  CHAMELEON_zgetrs_nopiv_Tile - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Tile equivalent of CHAMELEON_zgetrs_nopiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by CHAMELEON_zgetrf_nopiv.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrs_nopiv
 * @sa CHAMELEON_zgetrs_nopiv_Tile_Async
 * @sa CHAMELEON_cgetrs_nopiv_Tile
 * @sa CHAMELEON_dgetrs_nopiv_Tile
 * @sa CHAMELEON_sgetrs_nopiv_Tile
 * @sa CHAMELEON_zgetrf_nopiv_Tile
 *
 */
int CHAMELEON_zgetrs_nopiv_Tile( CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zgetrs_nopiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    CHAMELEON_zgetrs_nopiv_Tile_Async( A, B, sequence, &request );

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
 *  CHAMELEON_zgetrs_nopiv_Tile_Async - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Non-blocking equivalent of CHAMELEON_zgetrs_nopiv_Tile().
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
 * @sa CHAMELEON_zgetrs_nopiv
 * @sa CHAMELEON_zgetrs_nopiv_Tile
 * @sa CHAMELEON_cgetrs_nopiv_Tile_Async
 * @sa CHAMELEON_dgetrs_nopiv_Tile_Async
 * @sa CHAMELEON_sgetrs_nopiv_Tile_Async
 * @sa CHAMELEON_zgetrf_nopiv_Tile_Async
 *
 */
int CHAMELEON_zgetrs_nopiv_Tile_Async( CHAM_desc_t *A, CHAM_desc_t *B,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zgetrs_nopiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("CHAMELEON_zgetrs_nopiv_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("CHAMELEON_zgetrs_nopiv_Tile", "NULL request");
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
        morse_error("CHAMELEON_zgetrs_nopiv_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zgetrs_nopiv_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("CHAMELEON_zgetrs_nopiv_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(N, NRHS) == 0)
     return CHAMELEON_SUCCESS;
     */
    morse_pztrsm( ChamLeft, ChamLower, ChamNoTrans, ChamUnit, (CHAMELEON_Complex64_t)1.0, A, B, sequence, request );

    morse_pztrsm( ChamLeft, ChamUpper, ChamNoTrans, ChamNonUnit, (CHAMELEON_Complex64_t)1.0, A, B, sequence, request );

    return CHAMELEON_SUCCESS;
}
