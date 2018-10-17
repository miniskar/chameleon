/**
 *
 * @file time_zpotrf_tile.c
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

#define _NAME  "CHAMELEON_zpotrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( N )
#define _FADDS FADDS_POTRF( N )

#include "./timing.c"

static int
RunTest( int *iparam, double *dparam, chameleon_time_t *t_ )
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    cham_uplo_t uplo = ChamUpper;

    LDA = chameleon_max( LDA, N );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,  1,     CHAMELEON_Complex64_t, ChamComplexDouble, LDA, N, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  check, CHAMELEON_Complex64_t, ChamComplexDouble, LDB, N, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAC, check, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, N, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  check, CHAMELEON_Complex64_t, ChamComplexDouble, LDB, N, NRHS );

    /* Initialize data and save A if check */
    if ( check ) {
        CHAMELEON_zplghe_Tile( (double)N, ChamUpperLower, descAC, 51 );
        CHAMELEON_zlacpy_Tile( uplo, descAC, descA );
    }
    else {
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA, 51 );
    }
    //RUNTIME_zlocality_allrestrict( STARPU_CUDA );

    /* CHAMELEON ZPOTRF */
    START_TIMING();
    CHAMELEON_zpotrf_Tile( uplo, descA );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        /* Initialize and save B */
        CHAMELEON_zplrnt_Tile( descB, 7672 );
        CHAMELEON_zlacpy_Tile( ChamUpperLower, descB, descX );

        /* Compute the solution */
        CHAMELEON_zpotrs_Tile( uplo, descA, descX );

        /* Check solution */
        dparam[IPARAM_ANORM] = CHAMELEON_zlange_Tile( ChamInfNorm, descAC );
        dparam[IPARAM_BNORM] = CHAMELEON_zlange_Tile( ChamInfNorm, descB  );
        dparam[IPARAM_XNORM] = CHAMELEON_zlange_Tile( ChamInfNorm, descX  );
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = CHAMELEON_zlange_Tile( ChamInfNorm, descB );

        PASTE_CODE_FREE_MATRIX( descB  );
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descX  );

    }
    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
