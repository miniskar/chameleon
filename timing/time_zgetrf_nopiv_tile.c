/**
 *
 * @file time_zgetrf_nopiv_tile.c
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
 * @precisions normal z -> c d s
 *
 */
#define _TYPE  CHAMELEON_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "CHAMELEON_zgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, N)
#define _FADDS FADDS_GETRF(M, N)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_) 
{
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

    /* Save A for check */
    if (check == 1){
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descA, descAC);
    }

    /**
     * Consider this optimization on some heterogenous platforms and matrix sizes.
     * Often, TRSM kernel on GPU yields significantly less performance rate than GEMM,
     * while performances are similar on CPU. On this algorithm it is therefore
     * recommended to execute all TRSMs (~low amount) on CPU to increase GPU efficiency.
     */
    //RUNTIME_zlocality_onerestrict( CHAMELEON_TRSM, STARPU_CPU );

    START_TIMING();
    CHAMELEON_zgetrf_nopiv_Tile( descA );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        CHAMELEON_zplrnt_Tile( descX, 7732 );
        CHAMELEON_zlacpy_Tile(ChamUpperLower, descX, descB);

        CHAMELEON_zgetrs_nopiv_Tile( descA, descX );

        dparam[IPARAM_ANORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descAC);
        dparam[IPARAM_BNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        dparam[IPARAM_XNORM] = CHAMELEON_zlange_Tile(ChamInfNorm, descX);
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, 1.0, descAC, descX, -1.0, descB );
        dparam[IPARAM_RES] = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        PASTE_CODE_FREE_MATRIX( descX  );
        PASTE_CODE_FREE_MATRIX( descAC );
        PASTE_CODE_FREE_MATRIX( descB  );
    }

    PASTE_CODE_FREE_MATRIX( descA );
 
    return 0;
}
