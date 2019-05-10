/**
 *
 * @file testing_dgram.c
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dgram testing
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2019-04-12
 * @precisions normal d -> d s
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

static int check_solution(cham_uplo_t uplo,
                          int N,
                          double *Aref,
                          double *Acham, int LDA);
static int compute_gram_sequential(cham_uplo_t uplo,
                                   int N,
                                   double *A,
                                   int LDA);

int testing_dgram(int argc, char **argv)
{
    int hres = 0;
    /* Check for number of arguments */
    if ( argc < 2) {
        USAGE("GRAM", "N LDA",
              "   - N      : number of rows of matrix A\n"
              "   - LDA    : leading dimension of matrix A\n");
        return -1;
    }
    int N     = atoi(argv[0]);
    int LDA   = atoi(argv[1]);

    double eps;
    int info_solution;
    int i, j, ua;
    int LDAxN = LDA*N;

    double *A     = (double *)malloc(LDAxN*sizeof(double));
    double *Aref  = (double *)malloc(LDAxN*sizeof(double));
    double *Acham = (double *)malloc(LDAxN*sizeof(double));

    /* Check if unable to allocate memory */
    if ( (!A) || (!Aref) || (!Acham) )
    {
        free(A); free(Aref); free(Acham);
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    if (CHAMELEON_My_Mpi_Rank() == 0){
        printf("\n");
        printf("------ TESTS FOR CHAMELEON GRAM ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", N, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n",eps);
        printf(" Computational tests pass if scaled residuals are less than 10.\n");
    }

    /*----------------------------------------------------------
     *  TESTING GRAM
     */

    /* Initialize A such that it is symmetric */
    CHAMELEON_dplgsy( (double)N, ChamUpperLower, N, A, LDA, 51 );
    /* Gram is meant to be used with A full of positive values only */
#if defined(PRECISION_d) || defined(PRECISION_s)
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if ( A[i+j*LDA] < 0. ) {
                A[i+j*LDA] = -A[i+j*LDA];
            }
        }
    }
#endif

    for (ua=0; ua<3; ua++) {
        CHAMELEON_dlacpy( ChamUpperLower, N, N, A, LDA, Aref, LDA );
        CHAMELEON_dlacpy( ChamUpperLower, N, N, A, LDA, Acham, LDA );

        /* CHAMELEON GRAM */
        CHAMELEON_dgram(uplo[ua], N, Acham, LDA);

        /* Check the solution */
        info_solution = check_solution(uplo[ua], N, Aref, Acham, LDA);

        if (CHAMELEON_My_Mpi_Rank() == 0){
            if (info_solution == 0) {
                printf("***************************************************\n");
                printf(" ---- TESTING GRAM (%s) ............... PASSED !\n", uplostr[ua]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING GRAM (%s) ... FAILED !\n", uplostr[ua]);    hres++;
                printf("************************************************\n");
            }
        }
    }
    free(A); free(Aref); free(Acham);

    return hres;
}

/*--------------------------------------------------------------
 * Check the solution
 */
static int check_solution(cham_uplo_t uplo,
                          int N,
                          double *Aref,
                          double *Acham, int LDA)
{
    int info_solution;
    double Arefnorm, Rnorm, result;
    double eps;
    double mdone;

    double *work = (double *)malloc(N * sizeof(double));

    mdone = -1.0;

    /*
     * Compute the Gram matrix sequentially
     * we consider the matrix on entry as symmetric
     */
    compute_gram_sequential(uplo, N, Aref, LDA);

    /* Compute norm of Aref to scale the result norm */
    if (uplo == ChamUpperLower) {
        Arefnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I',
                                       N,  N,  Aref,  LDA, work);
    } else {
        Arefnorm = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'I',
                                       chameleon_lapack_const(uplo), chameleon_lapack_const(ChamNonUnit),
                                       N, N, Aref, LDA, work);
    }
    /* compute the difference Aref = Aref - Acham */
    cblas_daxpy(LDA*N, mdone, Acham, 1, Aref, 1);

    /* compute the norm of the difference */
    if (uplo == ChamUpperLower) {
        Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I',
                                    N, N, Aref, LDA, work);
    } else {
        Rnorm = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'I',
                                    chameleon_lapack_const(uplo), chameleon_lapack_const(ChamNonUnit),
                                    N, N, Aref, LDA, work);
    }

    eps = LAPACKE_dlamch_work('e');
    if (CHAMELEON_My_Mpi_Rank() == 0)
        printf("Rnorm %e, Anorm %e\n", Rnorm, Arefnorm);

    /* scale the norm in respect to Aref */
    result = Rnorm / (Arefnorm * N * eps);

    if (CHAMELEON_My_Mpi_Rank() == 0){
        printf("============\n");
        printf("Checking the norm of the difference against reference GRAM \n");
        printf("-- ||Acham - Aref||_oo/((||Aref||_oo.N.eps) = %e \n",
               result);
    }

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
         if (CHAMELEON_My_Mpi_Rank() == 0)
             printf("-- The solution is suspicious ! \n");
         info_solution = 1;
    }
    else {
         if (CHAMELEON_My_Mpi_Rank() == 0)
             printf("-- The solution is CORRECT ! \n");
         info_solution= 0 ;
    }

    free(work);

    return info_solution;
}

/*--------------------------------------------------------------
 * Compute the Gram matrix sequentially
 * We consider the matrix on entry as symmetric
 */
static int compute_gram_sequential(cham_uplo_t uplo,
                                   int N,
                                   double *A,
                                   int LDA)
{
    int m, n;
    double eps;
    double squareij, mean_dij, mhalf;

    double *work = (double *)malloc(N * sizeof(double));

    mhalf = -0.5;

    /* initialize work */
    memset(work, 0., N*sizeof(double));

    /* first: compute the means of squares */
    for (n=0; n<N; n++) {
        int mmin = ( uplo == ChamLower ) ? n                     : 0;
        int mmax = ( uplo == ChamUpper ) ? chameleon_min(n+1, N) : N;
        for (m = mmin; m < mmax; m++) {
            squareij = A[m+n*LDA]*A[m+n*LDA];
            /* accumulate squares on columns */
            work[n] += squareij;
            if ( m != n && uplo != ChamUpperLower ) {
                /* accumulate squares on the symmetric part */
                work[m] += squareij;
            }
        }
    }
    mean_dij = 0.;
    for (n=0; n<N; n++) {
        /* accumulate the squares over the entire matrix */
        mean_dij += work[n];
        /* compute the mean on each column */
        work[n] /= N;
    }
    /* compute the global mean */
    mean_dij /= N*N;
    /* second: compute the Gram matrix factors */
    for (n=0; n<N; n++) {
        int mmin = ( uplo == ChamLower ) ? n                     : 0;
        int mmax = ( uplo == ChamUpper ) ? chameleon_min(n+1, N) : N;
        for (m = mmin; m < mmax; m++) {
            squareij = A[m+n*LDA]*A[m+n*LDA];
            A[m+n*LDA] = mhalf*( squareij - work[m] - work[n] + mean_dij );
        }
    }

    free(work);

    return 0;
}