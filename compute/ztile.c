/**
 *
 * @file ztile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary routines
 *
 * @version 0.9.2
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2014-11-16
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zLapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the CHAMELEON matrix in tile layout.
 *          If CHAMELEON_TRANSLATION_MODE is set to ChamInPlace,
 *          A->mat is not used and set to Af77 when returns, else if
 *          CHAMELEON_TRANSLATION_MODE is set to ChamOutOfPlace,
 *          A->mat has to be allocated before.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zTile_to_Lapack
 * @sa CHAMELEON_cLapack_to_Tile
 * @sa CHAMELEON_dLapack_to_Tile
 * @sa CHAMELEON_sLapack_to_Tile
 *
 */
int CHAMELEON_zLapack_to_Tile( CHAMELEON_Complex64_t *Af77, int LDA, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request;
    CHAM_desc_t *B;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zLapack_to_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (chameleon_desc_check( A ) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zLapack_to_Tile", "invalid descriptor");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    CHAMELEON_Desc_Create_User( &B, Af77, ChamComplexDouble, A->mb, A->nb, A->bsiz,
                                LDA, A->n, 0, 0, A->m, A->n, 1, 1,
                                chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL );

    /* Start the computation */
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pzlacpy( ChamUpperLower, B, A, sequence, &request );

    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Destroy temporary B descriptor */
    CHAMELEON_Desc_Destroy( &B );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the CHAMELEON matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If CHAMELEON_TRANSLATION_MODE is set to ChamInPlace,
 *          Af77 has to be A->mat, else if
 *          CHAMELEON_TRANSLATION_MODE is set to ChamOutOfPlace,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zLapack_to_Tile
 * @sa CHAMELEON_cTile_to_Lapack
 * @sa CHAMELEON_dTile_to_Lapack
 * @sa CHAMELEON_sTile_to_Lapack
 *
 */
int CHAMELEON_zTile_to_Lapack( CHAM_desc_t *A, CHAMELEON_Complex64_t *Af77, int LDA )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request;
    CHAM_desc_t *B;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zTile_to_Lapack", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (chameleon_desc_check( A ) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zTile_to_Lapack", "invalid descriptor");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    CHAMELEON_Desc_Create_User( &B, Af77, ChamComplexDouble, A->mb, A->nb, A->bsiz,
                                LDA, A->n, 0, 0, A->m, A->n, 1, 1,
                                chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL );

    /* Start the computation */
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pzlacpy( ChamUpperLower, A, B, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    CHAMELEON_Desc_Destroy( &B );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
