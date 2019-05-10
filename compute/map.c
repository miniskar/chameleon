/**
 *
 * @file map.c
 *
 * @copyright 2018-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map wrappers
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @date 2018-09-24
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Tile
 *
 *  Apply a given operator on each tile of the given matrix. Operates on
 *  matrices stored by tiles.  All matrices are passed through descriptors.  All
 *  dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] uplo
 *          - ChamUpper: Only the upper triangular part of the matrix is touched
 *          - ChamLower: Only the lower triangular part of the matrix is touched
 *          - ChamUpperLower: The entire the matrix is touched
 *
 * @param[in,out] A
 *          On exit, the operator has been applied on each tile of the matrix A.
 *
 * @param[in] op_fct
 *          The operator function to apply on each tile of the matrix.
 *
 * @param[in,out] op_args
 *          The arguments structure passed to the operator function when applied
 *          on each tile. May be updated by the operator function.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_map_Tile_Async
 *
 */
int CHAMELEON_map_Tile( cham_uplo_t           uplo,
                        CHAM_desc_t          *A,
                        cham_unary_operator_t op_fct,
                        void                 *op_args )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_map_Tile_Async( uplo, A, op_fct, op_args, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Tile_Async
 *
 *  Apply a given operator on each tile of the given matrix. Non-blocking equivalent of
 *  CHAMELEON_map_Tile().  May return before the computation is finished.
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
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_map_Tile
 *
 */
int CHAMELEON_map_Tile_Async( cham_uplo_t           uplo,
                              CHAM_desc_t          *A,
                              cham_unary_operator_t op_fct,
                              void                 *op_args,
                              RUNTIME_sequence_t   *sequence,
                              RUNTIME_request_t    *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_map_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pmap( uplo, A, op_fct, op_args, sequence, request );

    return CHAMELEON_SUCCESS;
}
