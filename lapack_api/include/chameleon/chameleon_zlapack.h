/**
 *
 * @file chameleon_zlapack.h
 *
 * @copyright 2022-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas/lapack and cblas/lapack api functions
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2022-04-26
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_zlapack_h_
#define _chameleon_zlapack_h_

#include <chameleon.h>
#include <coreblas/cblas_wrapper.h>

BEGIN_C_DECLS

/**
 *  Declarations of math functions (LAPACK layout, Cblas/Lapacke interface) - alphabetical order
 */
void CHAMELEON_cblas_zgemm( const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
                            const int M, const int N, const int K,
                            const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                               const CHAMELEON_Complex64_t *B, const int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, const int ldc );

END_C_DECLS

#endif /* _chameleon_zlapack_h_ */
