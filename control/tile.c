/**
 *
 * @file tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon layout conversion wrappers
 *
 * @version 0.9.2
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @date 2014-11-16
 *
 ***
 *
 * @defgroup Tile
 * @brief Group routines exposed to users for matrices conversion LAPACK-Tile
 *
 */
#include "control/common.h"
#include "control/auxiliary.h"

/**
 *
 * @ingroup Tile
 *
 *  CHAMELEON_Lapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 ******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[out] A
 *          Descriptor of the CHAMELEON matrix in tile layout.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Lapack_to_Tile(void *Af77, int LDA, CHAM_desc_t *A)
{
    switch( A->dtyp ) {
    case ChamComplexDouble:
        return CHAMELEON_zLapack_to_Tile( (CHAMELEON_Complex64_t *)Af77, LDA, A );
        break;
    case ChamComplexFloat:
        return CHAMELEON_cLapack_to_Tile( (CHAMELEON_Complex32_t *)Af77, LDA, A );
        break;
    case ChamRealFloat:
        return CHAMELEON_sLapack_to_Tile( (float *)Af77, LDA, A );
        break;
    case ChamRealDouble:
    default:
        return CHAMELEON_dLapack_to_Tile( (double *)Af77, LDA, A );
    }
    return CHAMELEON_ERR_ILLEGAL_VALUE;
}

/**
 *
 * @ingroup Tile
 *
 *  CHAMELEON_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 ******************************************************************************
 *
 * @param[out] A
 *          Descriptor of the CHAMELEON matrix in tile layout.
 *
 * @param[in] Af77
 *          LAPACK matrix (only needed on proc 0).
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Tile_to_Lapack(CHAM_desc_t *A, void *Af77, int LDA)
{
    switch( A->dtyp ) {
    case ChamComplexDouble:
        return CHAMELEON_zTile_to_Lapack( A, (CHAMELEON_Complex64_t *)Af77, LDA );
        break;
    case ChamComplexFloat:
        return CHAMELEON_cTile_to_Lapack( A, (CHAMELEON_Complex32_t *)Af77, LDA );
        break;
    case ChamRealFloat:
        return CHAMELEON_sTile_to_Lapack( A, (float *)Af77, LDA );
        break;
    case ChamRealDouble:
    default:
        return CHAMELEON_dTile_to_Lapack( A, (double *)Af77, LDA );
    }
    return CHAMELEON_ERR_ILLEGAL_VALUE;
}
