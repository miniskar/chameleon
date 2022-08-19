/**
 *
 * @file lapack_api_common.h
 *
 * @copyright 2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux.
 *                 All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas/lapack and cblas/lapack api common internal functions
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2022-04-22
 *
 */
#ifndef _lapack_api_common_h_
#define _lapack_api_common_h_

#include <stdio.h>
#include <chameleon.h>
#include <chameleon/mangling.h>
#include <coreblas/cblas_wrapper.h>

/**
 *
 * @defgroup CHAMELEON_LAPACK_API
 * @brief Linear algebra routines exposed to users. LAPACK matrix data storage
 *
 */

int chameleon_blastocblas_trans(const char* value);
int chameleon_blastocblas_side(const char* side);
int chameleon_blastocblas_uplo(const char* uplo);

#endif /* _lapack_api_common_h_ */
