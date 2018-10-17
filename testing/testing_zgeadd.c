/**
 *
 * @file testing_zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd testing
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <chameleon.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <coreblas.h>
#include "testing_zauxiliary.h"
#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

static int check_tr_solution(cham_uplo_t uplo, cham_trans_t trans, int M, int N,
                             CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                             CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *Bcham, int LDB);
static int check_ge_solution(cham_trans_t trans, int M, int N,
                             CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                             CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *Bcham, int LDB);

int testing_zgeadd(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 6 ) {
        USAGE("GEADD", "alpha beta M N LDA LDB",
              "   - alpha  : alpha coefficient\n"
              "   - beta   : beta coefficient\n"
              "   - M      : number of rows of matrices A and C\n"
              "   - N      : number of columns of matrices B and C\n"
              "   - LDA    : leading dimension of matrix A\n"
              "   - LDB    : leading dimension of matrix B\n" );
        return -1;
    }

    CHAMELEON_Complex64_t alpha = (CHAMELEON_Complex64_t) atol(argv[0]);
    CHAMELEON_Complex64_t beta  = (CHAMELEON_Complex64_t) atol(argv[1]);
    int M   = atoi(argv[2]);
    int N   = atoi(argv[3]);
    int LDA = atoi(argv[4]);
    int LDB = atoi(argv[5]);

    double eps;
    int info_solution;
    int t, u;
    int LDAxN = LDA*max(M,N);
    int LDBxN = LDB*N;

    CHAMELEON_Complex64_t *A      = (CHAMELEON_Complex64_t *)malloc(LDAxN*sizeof(CHAMELEON_Complex64_t));
    CHAMELEON_Complex64_t *B      = (CHAMELEON_Complex64_t *)malloc(LDBxN*sizeof(CHAMELEON_Complex64_t));
    CHAMELEON_Complex64_t *Binit  = (CHAMELEON_Complex64_t *)malloc(LDBxN*sizeof(CHAMELEON_Complex64_t));
    CHAMELEON_Complex64_t *Bfinal = (CHAMELEON_Complex64_t *)malloc(LDBxN*sizeof(CHAMELEON_Complex64_t));

    /* Check if unable to allocate memory */
    if ( (!A) || (!B) || (!Binit) || (!Bfinal) )
    {
        free(A); free(B);
        free(Binit); free(Bfinal);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZGEADD ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrices A and B are randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING ZGEADD
     */

    /* Initialize A, B */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A);
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxN, B);

#if defined(PRECISION_z) || defined(PRECISION_c)
    for (t=0; t<3; t++) {
#else
    for (t=0; t<2; t++) {
#endif
        memcpy( Binit,  B, LDBxN*sizeof(CHAMELEON_Complex64_t));
        memcpy( Bfinal, B, LDBxN*sizeof(CHAMELEON_Complex64_t));

        /* CHAMELEON ZGEADD */
        CHAMELEON_zgeadd(trans[t], M, N, alpha, A, LDA, beta, Bfinal, LDB);

        /* Check the solution */
        info_solution = check_ge_solution(trans[t], M, N,
                                          alpha, A, LDA,
                                          beta, Binit, Bfinal, LDB);

        if (info_solution == 0) {
            printf("***************************************************\n");
            printf(" ---- TESTING ZGEADD (%s) ............... PASSED !\n", transstr[t]);
            printf("***************************************************\n");
        }
        else {
            printf("************************************************\n");
            printf(" - TESTING ZGEADD (%s) ... FAILED !\n", transstr[t]);
            printf("************************************************\n");
        }
    }
#ifdef _UNUSED_
    }
#endif

    /*----------------------------------------------------------
     *  TESTING TRADD
     */

    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A);
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxN, B);

#if defined(PRECISION_z) || defined(PRECISION_c)
    for (t=0; t<3; t++) {
#else
    for (t=0; t<2; t++) {
#endif
        for (u=0; u<2; u++) {
            memcpy( Binit,  B, LDBxN*sizeof(CHAMELEON_Complex64_t));
            memcpy( Bfinal, B, LDBxN*sizeof(CHAMELEON_Complex64_t));

            /* CHAMELEON ZGEADD */
            CHAMELEON_ztradd(uplo[u], trans[t], M, N, alpha, A, LDA, beta, Bfinal, LDB);

            /* Check the solution */
            info_solution = check_tr_solution(uplo[u], trans[t], M, N,
                                              alpha, A, LDA,
                                              beta, Binit, Bfinal, LDB);

            if (info_solution == 0) {
                printf("***************************************************\n");
                printf(" ---- TESTING ZTRADD (%s, %s) ............... PASSED !\n", uplostr[u], transstr[t]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING ZTRADD (%s, %s) ... FAILED !\n", uplostr[u], transstr[t]);
                printf("************************************************\n");
            }
        }
    }
#ifdef _UNUSED_
    }
#endif

    free(A); free(B);
    free(Binit); free(Bfinal);

    return 0;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_tr_solution(cham_uplo_t uplo, cham_trans_t trans, int M, int N,
                             CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                             CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *Bcham, int LDB)
{
    int info_solution;
    double Anorm, Binitnorm, Bchamnorm, Rnorm, result;
    double eps;
    CHAMELEON_Complex64_t mzone;

    double *work = (double *)malloc(max(M, N)* sizeof(double));
    int Am, An;

    mzone = -1.0;

    if (trans == ChamNoTrans) {
        Am = M; An = N;
    } else {
        Am = N; An = M;
    }

    /* if ( ((trans == ChamNoTrans) && (uplo == ChamLower)) || */
    /*      ((trans != ChamNoTrans) && (uplo == ChamUpper)) ) */
    /* { */
    /*     Anorm = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', 'L', 'N', */
    /*                                 Am, An, A, LDA, work); */
    /* } */
    /* else */
    /* { */
    /*     Anorm = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', 'U', 'N', */
    /*                                 Am, An, A, LDA, work); */
    /* } */

    /* Binitnorm  = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo[u]), 'N', */
    /*                                  M, N, Bref,   LDB, work); */
    /* Bchamnorm = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo[u]), 'N', */
    /*                                  M, N, Bcham, LDB, work); */

    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                Am, An, A, LDA, work);
    Binitnorm   = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                      M, N, Bref,    LDB, work);
    Bchamnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                      M, N, Bcham, LDB, work);

    CORE_ztradd(uplo, trans, M, N,
                alpha, A,    LDA,
                beta,  Bref, LDB);
    cblas_zaxpy( LDB*N, CBLAS_SADDR(mzone), Bcham, 1, Bref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'M', M, N, Bref, LDB, work);

    eps = LAPACKE_dlamch_work('e');

    printf("Rnorm %e, Anorm %e, Bnorm %e, (alpha A + beta B) norm %e\n",
           Rnorm, Anorm, Binitnorm, Bchamnorm);

    result = Rnorm / (max(Anorm, Binitnorm) * eps);
    printf("============\n");
    printf("Checking the norm of the difference against reference ZGEADD \n");
    printf("-- || R||_max/(max(||A||_oo,||B||_oo).eps) = %e \n",
           result);

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
    }
    else {
        printf("-- The solution is CORRECT ! \n");
        info_solution= 0 ;
    }

    free(work);

    return info_solution;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_ge_solution(cham_trans_t trans, int M, int N,
                             CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                             CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *Bcham, int LDB)
{
    int info_solution;
    double Anorm, Binitnorm, Bchamnorm, Rnorm, result;
    double eps;
    CHAMELEON_Complex64_t mzone;

    double *work = (double *)malloc(max(M, N)* sizeof(double));
    int Am, An;

    mzone = -1.0;

    if (trans == ChamNoTrans) {
        Am = M; An = N;
    } else {
        Am = N; An = M;
    }

    Anorm      = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                     Am, An, A,      LDA, work);
    Binitnorm  = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                     M,  N,  Bref,   LDB, work);
    Bchamnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I',
                                     M,  N,  Bcham, LDB, work);

    CORE_zgeadd(trans, M, N,
                alpha, A,    LDA,
                beta,  Bref, LDB);
    cblas_zaxpy( LDB*N, CBLAS_SADDR(mzone), Bcham, 1, Bref, 1);

    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'M', M, N, Bref, LDB, work);

    eps = LAPACKE_dlamch_work('e');

    printf("Rnorm %e, Anorm %e, Bnorm %e, (alpha A + beta B) norm %e\n",
           Rnorm, Anorm, Binitnorm, Bchamnorm);

    result = Rnorm / (max(Anorm, Binitnorm) * eps);
    printf("============\n");
    printf("Checking the norm of the difference against reference ZGEADD \n");
    printf("-- || R||_max/(max(||A||_oo,||B||_oo).eps) = %e \n",
           result);

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
    }
    else {
        printf("-- The solution is CORRECT ! \n");
        info_solution= 0 ;
    }

    free(work);

    return info_solution;
}
