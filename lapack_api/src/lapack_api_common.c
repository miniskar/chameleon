/**
 *
 * @file lapack_api_common.c
 *
 * @copyright 2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                 Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas/lapack and cblas/lapack api common functions
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2022-04-22
 *
 */
#include "lapack_api_common.h"

/**
 * @brief Convert the input char BLAS trans parameter to a compatible parameter
 * for the Cblas API.
 * @param[in] trans The input char BLAS trans parameter
 * @return The CBLAS equivalent parameter (CblasNoTrans, CblasTrans or
 * CblasConjTrans).
 */
int chameleon_blastocblas_trans(const char* trans)
{
    if ( (*trans == 'N') || (*trans == 'n') ) {
        return CblasNoTrans;
    } else if ( (*trans == 'T') || (*trans == 't') ) {
        return CblasTrans;
    } else if ( (*trans == 'C') || (*trans == 'c') ) {
        return CblasConjTrans;
    } else {
        fprintf(stderr, "CHAMELEON ERROR: %s(): %s\n", "chameleon_blastocblas_trans", "illegal value of BLAS transpose parameter");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}
