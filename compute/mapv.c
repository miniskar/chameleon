/**
 *
 * @file mapv.c
 *
 * @copyright 2018-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon mapv wrappers
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2024-03-11
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
 * @param[in] access
 *          - ChamR: A is accessed in read-only mode.
 *          - ChamW: A is accessed in write-only mode.
 *           WARNING: if the descriptor is set for allocation on the fly, the
 *           flush call included in this synchronous API will free all allocated
 *           data, prefer asynchronous call if you want to initialiaze data
 *           before submitting another algorithm.
 *          - ChamRW: A is accessed in read-write mode.
 *
 * @param[in] uplo
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
 * @sa CHAMELEON_mapv_Tile_Async
 *
 */
int CHAMELEON_mapv_Tile( cham_uplo_t          uplo,
                         int                  ndata,
                         cham_map_data_t     *data,
                         cham_map_operator_t *op_fct,
                         void                *op_args )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
    int                 status, i;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_mapv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_mapv_Tile_Async( uplo, ndata, data, op_fct, op_args, sequence, &request );

    for( i=0; i<ndata; i++ ) {
        CHAMELEON_Desc_Flush( data[i].desc, sequence );
    }

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
 *  CHAMELEON_mapv_Tile().  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] access
 *          - ChamR: A is accessed in read-only mode.
 *          - ChamW: A is accessed in write-only mode.
 *          INFO: tile of A can be unallocated before the call if the
 *          descriptor is set for allocation on the fly.
 *          - ChamRW: A is accessed in read-write mode.
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
 * @sa CHAMELEON_mapv_Tile
 *
 */
int CHAMELEON_mapv_Tile_Async( cham_uplo_t          uplo,
                               int                  ndata,
                               cham_map_data_t     *data,
                               cham_map_operator_t *op_fct,
                               void                *op_args,
                               RUNTIME_sequence_t  *sequence,
                               RUNTIME_request_t   *request )
{
    CHAM_context_t *chamctxt;
    int             i;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_mapv_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_mapv_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_mapv_Tile_Async", "NULL request");
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
    for( i=0; i<ndata; i++ ) {
        if (chameleon_desc_check(data[i].desc) != CHAMELEON_SUCCESS) {
            chameleon_error("CHAMELEON_mapv_Tile_Async", "invalid descriptor");
            return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
        }
    }

    chameleon_pmap( uplo, ndata, data, op_fct, op_args, sequence, request );

    return CHAMELEON_SUCCESS;
}
