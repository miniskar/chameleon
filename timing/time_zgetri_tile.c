/**
 *
 * @file time_zgetri_tile.c
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

#define _NAME  "CHAMELEON_zgetri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF(M, N) + FMULS_GETRI( N ))
#define _FADDS (FADDS_GETRF(M, N) + FADDS_GETRI( N ))

//#define GETRI_SYNC

#include "./timing.c"

/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
#if 0
static int check_getri_factorization(CHAM_desc_t *descA1, CHAM_desc_t *descA2, int *IPIV)
{
    int info_factorization;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double *work = (double *)malloc((descA1->m)*sizeof(double));
    double eps = LAPACKE_dlamch_work('e');
    CHAM_desc_t        *descB, *descX;
    CHAMELEON_Complex64_t *b = (CHAMELEON_Complex64_t *)malloc((descA1->m)*sizeof(CHAMELEON_Complex64_t));
    CHAMELEON_Complex64_t *x = (CHAMELEON_Complex64_t *)malloc((descA1->m)*sizeof(CHAMELEON_Complex64_t));

    CHAMELEON_Desc_Create(&descB, b, ChamComplexDouble, descA1->mb, descA1->nb, descA1->bsiz,
                      descA1->m, 1, 0, 0, descA1->m, 1, 1, 1);
    CHAMELEON_Desc_Create(&descX, x, ChamComplexDouble, descA1->mb, descA1->nb, descA1->bsiz,
                      descA1->m, 1, 0, 0, descA1->m, 1, 1, 1);

    CHAMELEON_zplrnt_Tile( descX, 537 );
    CHAMELEON_zlacpy_Tile( ChamUpperLower, descX, descB);

    CHAMELEON_zgetrs_Tile( ChamNoTrans, descA2, IPIV, descX );

    Xnorm = CHAMELEON_zlange_Tile(ChamInfNorm, descX,  work);
    Anorm = CHAMELEON_zlange_Tile(ChamInfNorm, descA1, work);
    Bnorm = CHAMELEON_zlange_Tile(ChamInfNorm, descB,  work);

    CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans,
                       (CHAMELEON_Complex64_t)1.,  descA1, descX,
                       (CHAMELEON_Complex64_t)-1., descB);

    Rnorm = CHAMELEON_zlange_Tile(ChamInfNorm, descB, work);

    if (getenv("CHAMELEON_TESTING_VERBOSE"))
      printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*(descA1->m)*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The factorization is suspicious ! \n");
        info_factorization = 1;
     }
    else{
        printf("-- The factorization is CORRECT ! \n");
        info_factorization = 0;
    }
    free(x); free(b); free(work);
    CHAMELEON_Desc_Destroy(&descB);
    CHAMELEON_Desc_Destroy(&descX);

    return info_factorization;
}
#endif

/*------------------------------------------------------------------------
 *  Check the accuracy of the computed inverse
 */

static int check_getri_inverse(CHAM_desc_t *descA1, CHAM_desc_t *descA2, int *IPIV, double *dparam )
{
    double Rnorm, Anorm, Ainvnorm, result;
    double *W = (double *)malloc(descA1->n*sizeof(double));
    CHAMELEON_Complex64_t *work = (CHAMELEON_Complex64_t *)malloc(descA1->n*descA1->n*sizeof(CHAMELEON_Complex64_t));
    double eps = LAPACKE_dlamch_work('e');
    CHAM_desc_t        *descW;

    CHAMELEON_Desc_Create(&descW, work, ChamComplexDouble,  descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, descA1->n, 0, 0, descA1->m, descA1->n);

    CHAMELEON_zlaset_Tile( ChamUpperLower, (CHAMELEON_Complex64_t)0., (CHAMELEON_Complex64_t)1., descW);
    CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans,
                       (CHAMELEON_Complex64_t)-1., descA2, descA1,
                       (CHAMELEON_Complex64_t)1.,  descW);

    Anorm    = CHAMELEON_zlange_Tile(ChamInfNorm, descA1, W);
    Ainvnorm = CHAMELEON_zlange_Tile(ChamInfNorm, descA2, W);
    Rnorm    = CHAMELEON_zlange_Tile(ChamInfNorm, descW,  W);

    dparam[IPARAM_ANORM] = Anorm;
    dparam[IPARAM_BNORM] = Ainvnorm;

    result = Rnorm / ( (Anorm*Ainvnorm)*descA1->m*eps ) ;
    dparam[IPARAM_RES] = Rnorm;

    if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        dparam[IPARAM_XNORM] = -1.;
    }
    else{
        dparam[IPARAM_XNORM] = 0.;
    }

    CHAMELEON_Desc_Destroy(&descW);
    free(W);
    free(work);

    return CHAMELEON_SUCCESS;
}

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_)
{
    CHAM_desc_t descW;
    int ret = 0;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,      1, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA2, check, CHAMELEON_Complex64_t, ChamComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, N, 1 );

    CHAMELEON_Alloc_Workspace_zgetri_Tile_Async(descA, &descW);
    CHAMELEON_zplrnt_Tile( descA, 3453 );

    if ( check ) {
        CHAMELEON_zlacpy_Tile( ChamUpperLower, descA, descA2 );
    }

    /* CHAMELEON ZGETRF / ZTRTRI / ZTRSMRV  */
    {
#if defined(TRACE_BY_SEQUENCE)
        RUNTIME_sequence_t *sequence;
        RUNTIME_request_t request[4] = { RUNTIME_REQUEST_INITIALIZER,
                                         RUNTIME_REQUEST_INITIALIZER,
                                         RUNTIME_REQUEST_INITIALIZER,
                                         RUNTIME_REQUEST_INITIALIZER };

        CHAMELEON_Sequence_Create(&sequence);

        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            CHAMELEON_zgetrf_Tile_Async( descA, piv, sequence, &request[0] );
            CHAMELEON_Sequence_Wait(sequence);

            CHAMELEON_ztrtri_Tile_Async( ChamUpper, ChamNonUnit, descA, sequence, &request[1] );
            CHAMELEON_Sequence_Wait(sequence);

            CHAMELEON_ztrsmrv_Tile_Async( ChamRight, ChamLower, ChamNoTrans, ChamUnit,
                                      (CHAMELEON_Complex64_t) 1.0, descA, &descW,
                                      sequence, &request[2] );
            CHAMELEON_Sequence_Wait(sequence);

            CHAMELEON_zlaswpc_Tile_Async( descA, 1, descA->m, piv, -1,
                                      sequence, &request[3] );
            CHAMELEON_Desc_Flush( descA, sequence );
            CHAMELEON_Sequence_Wait(sequence);
            STOP_TIMING();

        } else {

            START_TIMING();
            CHAMELEON_zgetrf_Tile_Async( descA, piv, sequence, &request[0]);
            CHAMELEON_ztrtri_Tile_Async( ChamUpper, ChamNonUnit,
                                     descA, sequence, &request[1] );
            CHAMELEON_ztrsmrv_Tile_Async( ChamRight, ChamLower, ChamNoTrans, ChamUnit,
                                      (CHAMELEON_Complex64_t) 1.0,
                                      descA, &descW, sequence, &request[2] );
            CHAMELEON_zlaswpc_Tile_Async( descA, 1, descA->m, piv, -1,
                                      sequence, &request[3] );

            /* Wait for everything */
            CHAMELEON_Desc_Flush( descA, sequence );
            CHAMELEON_Sequence_Wait( sequence );
            STOP_TIMING();

        }

        CHAMELEON_Sequence_Destroy(sequence[0]);
        CHAMELEON_Sequence_Destroy(sequence[1]);
        CHAMELEON_Sequence_Destroy(sequence[2]);
        CHAMELEON_Sequence_Destroy(sequence[3]);

#else
        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            CHAMELEON_zgetrf_Tile(descA, piv);
            CHAMELEON_ztrtri_Tile(ChamUpper, ChamNonUnit, descA);
            CHAMELEON_ztrsmrv_Tile(ChamRight, ChamLower, ChamNoTrans, ChamUnit,
                                (CHAMELEON_Complex64_t) 1.0, descA, &descW);
            CHAMELEON_zlaswpc_Tile(descA, 1, descA->m, piv, -1);
            STOP_TIMING();

        } else {

            RUNTIME_sequence_t *sequence;
            RUNTIME_request_t request[2] = { RUNTIME_REQUEST_INITIALIZER,
                                             RUNTIME_REQUEST_INITIALIZER };

            CHAMELEON_Sequence_Create(&sequence);

            START_TIMING();
            CHAMELEON_zgetrf_Tile_Async(descA, piv, sequence, &request[0]);
            CHAMELEON_zgetri_Tile_Async(descA, piv, &descW, sequence, &request[1]);
            CHAMELEON_Desc_Flush( descA, sequence );
            CHAMELEON_Sequence_Wait(sequence);
            STOP_TIMING();

            CHAMELEON_Sequence_Destroy(sequence);
        }
#endif
    }

    /* Check the solution */
    if ( check )
    {
        ret = check_getri_inverse(descA2, descA, piv, dparam);

        PASTE_CODE_FREE_MATRIX( descA2 );
    }

    PASTE_CODE_FREE_MATRIX( descA );
    free(descW.mat);
    free( piv );

    return ret;
}
