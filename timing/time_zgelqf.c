/**
 *
 * @file time_zgelqf.c
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

#define _NAME  "CHAMELEON_zgelqf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GELQF(M, N)
#define _FADDS FADDS_GELQF(M, N)

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    CHAM_desc_t *T;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, CHAMELEON_Complex64_t, LDA, N );

    /* Initialize Data */
    CHAMELEON_zplrnt(M, N, A, LDA, 3456);

    /* Allocate Workspace */
    CHAMELEON_Alloc_Workspace_zgels(M, N, &T, P, Q);
    memset(T->mat, 0, (T->llm*T->lln)*sizeof(ChamComplexDouble));

    /* Save AT in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, CHAMELEON_Complex64_t, A, LDA, N );

    START_TIMING();
    CHAMELEON_zgelqf( M, N, A, LDA, T );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX( X, 1, CHAMELEON_Complex64_t, LDB, NRHS );
        CHAMELEON_zplrnt( N, NRHS, X, LDB, 5673 );
        PASTE_CODE_ALLOCATE_COPY( B, 1, CHAMELEON_Complex64_t, X, LDB, NRHS );

        CHAMELEON_zgelqs(M, N, NRHS, A, LDA, T, X, LDB);

        dparam[IPARAM_RES] = z_check_solution(M, N, NRHS, Acpy, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]),
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));

        free( Acpy );
        free( B );
        free( X );
      }

    /* Free Workspace */
    CHAMELEON_Dealloc_Workspace( &T );
    free( A );

    return 0;
}
