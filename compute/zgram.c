/**
 *
 * @file zgram.c
 *
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram wrappers
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> s d c z
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgram replace a general matrix by the Gram matrix inplace
 *
 *    \f[
 *    d_{i.}^2 = (1/n) \sum_j d_{ij}^2` and :math:`d_{..}^2 = (1/n^2) \sum_{i,j} d_{ij}^2 \\
 *    A_{i,j} = -(1/2) (d_{ij}^2 - d_{i.}^2 - d_{.j}^2 + d_{..}^2)
 *    \f]
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of lines and columns of the matrix A. N >= 0. When N = 0,
 *          the returned value is set to zero.
 *
 * @param[in] A
 *          The N-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
* @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgram_Tile
 * @sa CHAMELEON_zgram_Tile_Async
 * @sa CHAMELEON_sgram
 *
 */
int CHAMELEON_zgram( cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgram", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        chameleon_error("CHAMELEON_zgram", "illegal value of N");
        return -1;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zgram", "illegal value of LDA");
        return -3;
    }

    /* Quick return */
    if (N == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_DGEMM, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgram", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, &descAl, &descAt, ChamDescInout, uplo,
                         A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgram_Tile_Async( uplo, &descAt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInout, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgram_Tile - Tile equivalent of CHAMELEON_zgram().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The N-by-N matrix A.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgram
 * @sa CHAMELEON_zgram_Tile_Async
 * @sa CHAMELEON_sgram_Tile
 *
 */
int CHAMELEON_zgram_Tile( cham_uplo_t uplo, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgram_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgram_Tile_Async( uplo, A, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

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
 *  CHAMELEON_zgram_Tile_Async - Non-blocking equivalent of CHAMELEON_zgram_Tile().
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
 * @sa CHAMELEON_zgram
 * @sa CHAMELEON_zgram_Tile
 * @sa CHAMELEON_sgram_Tile_Async
 *
 */
int CHAMELEON_zgram_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *A,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgram_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgram_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgram_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zgram_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zgram_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzgram( uplo, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
