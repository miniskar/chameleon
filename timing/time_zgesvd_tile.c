/**
 *
 * @file time_zgesvd_tile.c
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

#define _NAME  "CHAMELEON_zheev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEBRD( M, N )
#define _FADDS FADDS_GEBRD( M, N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    CHAM_desc_t *descT;
    int jobu  = ChamVec;
    int jobvt = ChamVec;
    int INFO;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( VT, (jobvt == ChamVec), CHAMELEON_Complex64_t, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( U,  (jobu  == ChamVec), CHAMELEON_Complex64_t, M, M );
    PASTE_CODE_ALLOCATE_MATRIX( S, 1, double, N, 1 );

    /* Initialiaze Data */
    CHAMELEON_zplrnt_Tile(descA, 51 );

    /* Allocate Workspace */
    CHAMELEON_Alloc_Workspace_zgesvd(N, N, &descT, 1, 1);

    if ( jobu == ChamVec ) {
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', M, M, 0., 1., U,  M);
    }
    if ( jobvt == ChamVec ) {
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', N, N, 0., 1., VT, N);
    }

    START_TIMING();
    INFO = CHAMELEON_zgesvd_Tile(jobu, jobvt, descA, S, descT, U, M, VT, N);
    STOP_TIMING();

    if( INFO != 0 ) {
        printf(" ERROR OCCURED INFO %d\n",INFO);
    }

    /* DeAllocate Workspace */
    CHAMELEON_Dealloc_Workspace(&descT);

    if ( U != NULL ) {
        free( U );
    }
    if ( VT != NULL) {
        free( VT );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    free( S );

    (void)dparam;
    return 0;
}
