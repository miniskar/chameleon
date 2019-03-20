/**
 *
 * @file time_zgels.c
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

#define _NAME  "CHAMELEON_zgels"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GEQRF( M, N ) + FMULS_GEQRS( M, N, NRHS ))
#define _FADDS (FADDS_GEQRF( M, N ) + FADDS_GEQRS( M, N, NRHS ))

#include "./timing.c"
#include "timing_zauxiliary.h"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    CHAM_desc_t *T;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,    1,     CHAMELEON_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( x,    1,     CHAMELEON_Complex64_t, LDB, NRHS);
    PASTE_CODE_ALLOCATE_MATRIX( Acpy, check, CHAMELEON_Complex64_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( b,    check, CHAMELEON_Complex64_t, LDB, NRHS);

     /* Initialize Data */
    CHAMELEON_zplrnt( M, N,    A, LDA,  453 );
    CHAMELEON_zplrnt( M, NRHS, x, LDB, 5673 );

    CHAMELEON_Alloc_Workspace_zgels(M, N, &T, P, Q);

    /* Save A and b  */
    if (check) {
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N,    A, LDA, Acpy, LDA);
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, x, LDB, b,    LDB);
    }

    START_TIMING();
    CHAMELEON_zgels( ChamNoTrans, M, N, NRHS, A, LDA, T, x, LDB );
    STOP_TIMING();

    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = z_check_solution(M, N, NRHS, Acpy, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]),
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(b);
    }

    CHAMELEON_Dealloc_Workspace( &T );
    free( A );
    free( x );

    return 0;
}
