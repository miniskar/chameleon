/**
 *
 * @file time_ztrsm.c
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

#define _NAME  "CHAMELEON_ztrsm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_TRSM( ChamLeft, N, NRHS )
#define _FADDS FADDS_TRSM( ChamLeft, N, NRHS )

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    CHAMELEON_Complex64_t alpha;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    LDA = chameleon_max( LDA, N );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,      1, CHAMELEON_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( B,      1, CHAMELEON_Complex64_t, LDB, NRHS);
    PASTE_CODE_ALLOCATE_MATRIX( B2, check, CHAMELEON_Complex64_t, LDB, NRHS);

     /* Initialiaze Data */
    CHAMELEON_zplgsy( (CHAMELEON_Complex64_t)N, ChamUpperLower, N, A, LDA, 453 );
    CHAMELEON_zplrnt( N, NRHS, B, LDB, 5673 );
    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    alpha = 10.; /*alpha * N  /  2.;*/

    if (check)
    {
        memcpy(B2, B, LDB*NRHS*sizeof(CHAMELEON_Complex64_t));
    }

    START_TIMING();
    CHAMELEON_ztrsm( ChamLeft, ChamUpper, ChamNoTrans, ChamUnit,
                  N, NRHS, alpha, A, LDA, B, LDB );
    STOP_TIMING();

    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = z_check_trsm( ChamLeft, ChamUpper, ChamNoTrans, ChamUnit,
                                           N, NRHS,
                                           alpha, A, LDA, B, B2, LDB,
                                           &(dparam[IPARAM_ANORM]),
                                           &(dparam[IPARAM_BNORM]),
                                           &(dparam[IPARAM_XNORM]));
        free(B2);
    }

    free( A );
    free( B );

    return 0;
}
