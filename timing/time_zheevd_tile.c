/**
 *
 * @file time_zheevd_tile.c
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
 * @date 2016-12-09
 * @precisions normal z -> c d s
 *
 */
#define _TYPE  CHAMELEON_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "CHAMELEON_zheevd_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((double)N * (double)N * (double)N))
#define _FADDS ((2. / 3.) * ((double)N * (double)N * (double)N))

#include "./timing.c"
/* #include <mkl_service.h> */

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    CHAM_desc_t *descT;
    cham_uplo_t uplo = ChamLower;
    int vec  = ChamVec;
    int INFO;

    LDA = chameleon_max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, CHAMELEON_Complex64_t, LDA, N);
    PASTE_CODE_ALLOCATE_MATRIX( S, 1, double, N, 1 );

    /* Allocate Workspace */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, 51 );
    CHAMELEON_Alloc_Workspace_zheevd(N, N, &descT, 1, 1);

    START_TIMING();
    INFO = CHAMELEON_zheevd(vec, uplo, N, A, LDA, S, descT);
    STOP_TIMING();

    if (INFO != 0){
        printf(" ERROR OCCURED INFO %d\n", INFO);
    }

    /* DeAllocate Workspace */
    CHAMELEON_Dealloc_Workspace(&descT);

    free( A );
    (void)dparam;
    return 0;
}
