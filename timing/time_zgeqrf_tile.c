/**
 *
 * @file time_zgeqrf_tile.c
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

#define _NAME  "CHAMELEON_zgeqrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF( M, N )
#define _FADDS FADDS_GEQRF( M, N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    CHAM_desc_t *descT;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,  1, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  ( check && M == N ), CHAMELEON_Complex64_t, ChamComplexDouble, LDB, M, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA0, ( check && M == N ), CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  ( check && M == N ), CHAMELEON_Complex64_t, ChamComplexDouble, LDB, M, NRHS );

    CHAMELEON_zplrnt_Tile( descA, 5373 );

    /* Save A for check */
    if (check == 1 && M == N){
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descA, descA0);
    }

    /* Allocate Workspace */
    CHAMELEON_Alloc_Workspace_zgels_Tile(M, N, &descT, P, Q);

    /* CHAMELEON ZGEQRF */
    START_TIMING();
    CHAMELEON_zgeqrf_Tile( descA, descT );
    STOP_TIMING();

    /* Check the solution */
    if ( check && M == N )
    {
        /* Initialize and save B */
        CHAMELEON_zplrnt_Tile( descX, 2264 );
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descX, descB);

        /* Compute the solution */
        CHAMELEON_zgeqrs_Tile( descA, descT, descX );

        /* Check solution */
        dparam[IPARAM_ANORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descA0);
        dparam[IPARAM_BNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        dparam[IPARAM_XNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descX);
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, 1.0, descA0, descX, -1.0, descB );
        dparam[IPARAM_RES] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descX  )
        PASTE_CODE_FREE_MATRIX( descA0 )
        PASTE_CODE_FREE_MATRIX( descB  )
    }

    /* Free data */
    CHAMELEON_Dealloc_Workspace(&descT);
    PASTE_CODE_FREE_MATRIX( descA )

    return 0;
}
