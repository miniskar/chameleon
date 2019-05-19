/**
 *
 * @file time_zhemm_tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @date 2014-11-16
 * @precisions normal z -> c
 *
 */
#define _TYPE  CHAMELEON_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "CHAMELEON_zhemm_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_HEMM( ChamLeft, M, N )
#define _FADDS FADDS_HEMM( ChamLeft, M, N )

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    CHAMELEON_Complex64_t alpha, beta;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    LDA = chameleon_max(M, iparam[IPARAM_LDA]);
    LDB = chameleon_max(M, iparam[IPARAM_LDB]);
    LDC = chameleon_max(M, iparam[IPARAM_LDC]);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, M );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, CHAMELEON_Complex64_t, ChamComplexDouble, LDB, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descC, 1, CHAMELEON_Complex64_t, ChamComplexDouble, LDC, M, N );

    /* Initialize Data */
    CHAMELEON_zplghe_Tile( 0, ChamUpper, descA, 5373 );
    CHAMELEON_zplrnt_Tile( descB, 7672 );
    CHAMELEON_zplrnt_Tile( descC, 6387 );

#if !defined(CHAMELEON_SIMULATION)
    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_zlarnv_work(1, ISEED, 1, &beta);
#else
    alpha = 1.5;
    beta = -2.3;
#endif

    /* Save C for check */
    PASTE_TILE_TO_LAPACK( descC, C2, check, CHAMELEON_Complex64_t, LDC, N );

    START_TIMING();
    CHAMELEON_zhemm_Tile( ChamLeft, ChamUpper, alpha, descA, descB, beta, descC );
    STOP_TIMING();

#if !defined(CHAMELEON_SIMULATION)
    /* Check the solution */
    if (check)
    {
        PASTE_TILE_TO_LAPACK( descA, A, check, CHAMELEON_Complex64_t, LDA, M );
        PASTE_TILE_TO_LAPACK( descB, B, check, CHAMELEON_Complex64_t, LDB, N );
        PASTE_TILE_TO_LAPACK( descC, C, check, CHAMELEON_Complex64_t, LDC, N );

        dparam[IPARAM_RES] = z_check_hemm( ChamLeft, ChamUpper, M, N,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]),
                                           &(dparam[IPARAM_BNORM]),
                                           &(dparam[IPARAM_XNORM]) );

        free(A); free(B); free(C); free(C2);
    }
#endif

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    PASTE_CODE_FREE_MATRIX( descC );
    return 0;
}
