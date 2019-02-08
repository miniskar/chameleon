/**
 *
 * @file time_zgeqrs_tile.c
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

#define _NAME  "CHAMELEON_zgeqrs_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRS( M, N, NRHS )
#define _FADDS FADDS_GEQRS( M, N, NRHS )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    CHAM_desc_t *descT;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    check = 1;
    M = N;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,  1, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  ( check && M == N ), CHAMELEON_Complex64_t, ChamComplexDouble, LDB, M, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAC, ( check && M == N ), CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  ( check && M == N ), CHAMELEON_Complex64_t, ChamComplexDouble, LDB, M, NRHS );

    CHAMELEON_zplrnt_Tile( descA, 5373 );

    /* Save A for check */
    if (check == 1 && M == N){
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descA, descAC);
    }

    /* Allocate Workspace */
    CHAMELEON_Alloc_Workspace_zgels_Tile(M, N, &descT, P, Q);

    /* CHAMELEON ZGEQRF */
    CHAMELEON_zgeqrf_Tile( descA, descT );

    /* Check the solution */
    if ( check && M == N )
    {
        /* Initialize and save B */
        CHAMELEON_zplrnt_Tile( descX, 2264 );
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descX, descB);

        /* Compute the solution */
        START_TIMING();
        CHAMELEON_zgeqrs_Tile( descA, descT, descX );
        STOP_TIMING();

        /* Check solution */
        dparam[IPARAM_ANORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descAC);
        dparam[IPARAM_BNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        dparam[IPARAM_XNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descX);
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descX  );
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descB  )
    }

    /* Free data */
    CHAMELEON_Dealloc_Workspace(&descT);
    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
