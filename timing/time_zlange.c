/**
 *
 * @file time_zlange.c
 *
 * @Copyright 2009-2014 The University of Tennessee and The University
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

#define _NAME  "CHAMELEON_zlange"
#define _FMULS FMULS_LANGE(M, N)
#define _FADDS FADDS_LANGE(M, N)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, chameleon_time_t *t_)
{
    int    hres = 0;
    double normcham, normlapack, result;
    int    norm = ChamInfNorm;

    PASTE_CODE_IPARAM_LOCALS( iparam );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, CHAMELEON_Complex64_t, M, N );

    CHAMELEON_zplrnt( M, N, A, LDA, 3436 );

    /* CHAMELEON ZLANGE */
    START_TIMING();
    normcham = CHAMELEON_zlange(norm, M, N, A, LDA);
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        double *work = (double*) malloc(chameleon_max(M,N)*sizeof(double));
        normlapack = LAPACKE_zlange_work(LAPACK_COL_MAJOR, chameleon_lapack_const(norm), M, N, A, LDA, work);
        result = fabs(normcham - normlapack);
        switch(norm) {
        case ChamMaxNorm:
            /* result should be perfectly equal */
            break;
        case ChamInfNorm:
            /* Sum order on the line can differ */
            result = result / (double)N;
            break;
        case ChamOneNorm:
            /* Sum order on the column can differ */
            result = result / (double)M;
            break;
        case ChamFrobeniusNorm:
            /* Sum oreder on every element can differ */
            result = result / ((double)M * (double)N);
            break;
        }
        if ( CHAMELEON_Comm_rank() == 0 ) {
            dparam[IPARAM_ANORM] = normlapack;
            dparam[IPARAM_BNORM] = 0.;
            dparam[IPARAM_XNORM] = 1.;
            dparam[IPARAM_RES] = result;
        }
        free( work );
    }

    free( A );

    return hres;
}
