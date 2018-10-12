/**
 *
 * @file time_zposv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#define _TYPE  CHAMELEON_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "CHAMELEON_zposv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRS( N, NRHS ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRS( N, NRHS ))

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    cham_uplo_t uplo = ChamUpper;

    LDA = chameleon_max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,  1,     CHAMELEON_Complex64_t, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( A2, check, CHAMELEON_Complex64_t, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( X,  1,     CHAMELEON_Complex64_t, LDB, NRHS );

    /* Initialize data and save A if check */
    if ( check ) {
        CHAMELEON_zplghe( (double)N, ChamUpperLower, N, A2, LDA, 51 );
        CHAMELEON_zlacpy( uplo, N, N, A2, LDA, A, LDA );
    }
    else {
        CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, 51 );
    }
    CHAMELEON_zplrnt( N, NRHS, X, LDB, 5673 );

    /* Save b  */
    PASTE_CODE_ALLOCATE_COPY( B, check, CHAMELEON_Complex64_t, X, LDB, NRHS );

    /* CHAMELEON ZPOSV */
    START_TIMING();
    CHAMELEON_zposv( uplo, N, NRHS, A, LDA, X, LDB );
    STOP_TIMING();

    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = z_check_solution( N, N, NRHS, A2, LDA, B, X, LDB,
                                               &(dparam[IPARAM_ANORM]),
                                               &(dparam[IPARAM_BNORM]),
                                               &(dparam[IPARAM_XNORM]) );
        free(A2); free(B);
    }

    free(A); free(X);
    return 0;
}
