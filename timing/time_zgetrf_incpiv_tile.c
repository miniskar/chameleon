/**
 *
 * @file time_zgetrf_incpiv_tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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

#define _NAME  "CHAMELEON_zgetrf_incpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, N)
#define _FADDS FADDS_GETRF(M, N)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    CHAM_desc_t *descL;
    int *piv;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX,  check, CHAMELEON_Complex64_t, ChamComplexDouble, LDB, M, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAC, check, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB,  check, CHAMELEON_Complex64_t, ChamComplexDouble, LDB, M, NRHS );

    CHAMELEON_zplrnt_Tile(descA, 3456);

    /* Allocate Workspace */
    CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile(chameleon_min(M,N), &descL, &piv, P, Q);

    /* Save A for check */
    if (check == 1){
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descA, descAC);
    }

    START_TIMING();
    CHAMELEON_zgetrf_incpiv_Tile( descA, descL, piv );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        CHAMELEON_zplrnt_Tile( descX, 7732 );
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descX, descB);

        CHAMELEON_zgetrs_incpiv_Tile( descA, descL, piv, descX );

        dparam[IPARAM_ANORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descAC);
        dparam[IPARAM_BNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        dparam[IPARAM_XNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descX);
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descX  );
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descB  );
    }

    /* Deallocate Workspace */
    CHAMELEON_Dealloc_Workspace(&descL);

    PASTE_CODE_FREE_MATRIX( descA );
    free( piv );

    return 0;
}
