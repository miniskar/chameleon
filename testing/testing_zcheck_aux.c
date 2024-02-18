/**
 *
 * @file testing_zcheck_aux.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Nathalie Furmento
 * @author Alycia Lisito
 * @author Philippe Swartvagher
 * @date 2023-01-05
 * @precisions normal z -> c d s
 *
 */
#include "../control/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <chameleon.h>

#if !defined(CHAMELEON_SIMULATION)

#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <coreblas.h>
#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

#ifndef max
#define max( _a_, _b_ ) ( (_a_) > (_b_) ? (_a_) : (_b_) )
#endif

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two matrices by their norms.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] M
 *          The number of rows of the matrices A and B.
 *
 * @param[in] N
 *          The number of columns of the matrices A and B.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] B
 *          The matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the matrix B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zmatrices_std( run_arg_list_t *args, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB )
{
    int info_solution        = 0;
    double Anorm, Rnorm, result;
    cham_fixdbl_t eps = testing_getaccuracy();

    double *work = (double *)malloc( LDA*N*sizeof(double) );

    /* Computes the norms */
    if ( uplo == ChamUpperLower ) {
        Anorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, A, LDA, work );
    }
    else {
        Anorm = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), 'N',
                                     M, N, A, LDA, work );
    }

    /* Computes the difference with the core function */
    CORE_zgeadd( ChamNoTrans, M, N, 1, A, LDA, -1, B, LDB );

    /* Computes the residual's norm */
    if ( uplo == ChamUpperLower ) {
        Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, B, LDB, work );
    }
    else {
        Rnorm = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), 'N',
                                     M, N, B, LDB, work );
    }
    if ( Anorm != 0. ) {
        result = Rnorm / (Anorm * eps);
    }
    else {
        result = Rnorm;
    }

    /* Verifies if the result is inside a threshold */
    if ( isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||B||", Rnorm );

    free(work);

    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two matrices by their norms.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zmatrices( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descB )
{
    int info_solution        = 0;
    int M                    = descA->m;
    int N                    = descB->n;
    int LDA                  = M;
    int LDB                  = M;
    int rank                 = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A = NULL;
    CHAMELEON_Complex64_t *B = NULL;

    if ( rank == 0 ) {
        A = (CHAMELEON_Complex64_t *)malloc((size_t)(LDA) * (size_t)(N) * sizeof(CHAMELEON_Complex64_t));
        B = (CHAMELEON_Complex64_t *)malloc((size_t)(LDB) * (size_t)(N) * sizeof(CHAMELEON_Complex64_t));
        if ( (A == NULL) || (B == NULL) ) {
            free( A );
            free( B );
            return CHAMELEON_ERR_OUT_OF_RESOURCES;
        }
    }

    /* Converts the matrices to LAPACK layout in order to compare them on the main process */
    CHAMELEON_zDesc2Lap( uplo, descA, A, LDA );
    CHAMELEON_zDesc2Lap( uplo, descB, B, LDB );

    if ( rank == 0 ) {
        info_solution = check_zmatrices_std( args, uplo, M, N, A, LDA, B, LDB );
        free( A );
        free( B );
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed norms.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] norm_type
 *          Whether the norm is a Max, One, Inf or Frobenius norm.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] diag
 *          Whether it is a unitary diagonal matrix or not.
 *
 * @param[in] norm_cham
 *          The computed norm.
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_znorm_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_normtype_t norm_type, cham_uplo_t uplo,
                     cham_diag_t diag, double norm_cham, int M, int N, CHAMELEON_Complex64_t *A, int LDA )
{
    int info_solution  = 0;
    double *work       = (double*) malloc(chameleon_max(M, N)*sizeof(double));
    double norm_lapack;
    double result;
    cham_fixdbl_t eps = testing_getaccuracy();

    /* Computes the norm with the LAPACK function */
    switch (matrix_type) {
    case ChamGeneral:
        norm_lapack = LAPACKE_zlange_work( LAPACK_COL_MAJOR, chameleon_lapack_const(norm_type), M, N, A, LDA, work );
        break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        norm_lapack = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, chameleon_lapack_const(norm_type), chameleon_lapack_const(uplo), M, A, LDA, work );
        break;
#endif
    case ChamSymmetric:
        norm_lapack = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, chameleon_lapack_const(norm_type), chameleon_lapack_const(uplo), M, A, LDA, work );
        break;
    case ChamTriangular:
        norm_lapack = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, chameleon_lapack_const(norm_type), chameleon_lapack_const(uplo), chameleon_lapack_const(diag), M, N, A, LDA, work );
        break;
    default:
        fprintf(stderr, "check_znorm: matrix_type(%d) unsupported\n", matrix_type );
        free( work );
        return 1;
    }

    /* Compares the norms */
    result = fabs( norm_cham - norm_lapack ) / ( norm_lapack * eps );

    run_arg_add_double( args, "||A||", norm_cham );
    run_arg_add_double( args, "||B||", norm_lapack );

    switch(norm_type) {
    case ChamInfNorm:
        /* Sum order on the line can differ */
        result = result / (double)N;
        break;
    case ChamOneNorm:
        /* Sum order on the column can differ */
        result = result / (double)M;
        break;
    case ChamFrobeniusNorm:
        /* Sum order on every element can differ */
        result = result / ((double)M * (double)N);
        break;
    case ChamMaxNorm:
    default:
#if defined(PRECISION_z) || defined(PRECISION_c)
        /* Add a factor two to macth different implementation of cabs */
        result = result / 2;
#else
        /* result should be perfectly equal */
#endif
        ;
    }

    run_arg_add_double( args, "||R||", result );
    info_solution = ( result < 1 ) ? 0 : 1;

    free(work);

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares the Chameleon computed norm with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] norm_type
 *          Whether the norm is a Max, One, Inf or Frobenius norm.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] diag
 *          Whether it is a unitary diagonal matrix or not.
 *
 * @param[in] norm_cham
 *          The Chameleon computed norm.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_znorm( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_normtype_t norm_type, cham_uplo_t uplo,
                 cham_diag_t diag, double norm_cham, CHAM_desc_t *descA )
{
    int info_solution        = 0;
    int M                    = descA->m;
    int N                    = descA->n;
    int rank                 = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A = NULL;
    int LDA                  = M;

    if ( rank == 0 ) {
        A = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
    }

    /* Converts the matrix to LAPACK layout in order to use the LAPACK norm function */
    CHAMELEON_zDesc2Lap( uplo, descA, A, LDA );
    if ( rank == 0 ) {
        info_solution = check_znorm_std( args, matrix_type, norm_type, uplo, diag, norm_cham, M, N, A, LDA );
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed sum.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] trans
 *          Whether the first matrix is transposed, conjugate transposed or not transposed.
 *
 * @param[in] M
 *          The number of rows of the matrices A and C.
 *
 * @param[in] N
 *          The number of columns of the matrices B and C.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] Bref
 *          The matrix Bref.
 *
 * @param[in] Bcham
 *          The matrix of the computed result A+B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsum_std( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A,
                     int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *Bcham, int LDB )
{
    int info_solution            = 0;
    int Am                       = (trans == ChamNoTrans)? M : N;
    int An                       = (trans == ChamNoTrans)? N : M;
    double Anorm, Binitnorm, Rnorm, result;
    CHAMELEON_Complex64_t  mzone = -1.0;
    cham_uplo_t uploA            = uplo;

    /* Computes the max norms of A, B and A+B */
    if ( uplo == ChamUpperLower ) {
        Anorm     = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', Am, An, A,    LDA );
        Binitnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', M,  N,  Bref, LDB );
    }
    else {
        if ( trans == ChamNoTrans ) {
            Anorm = LAPACKE_zlantr( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), 'N', Am, An, A, LDA );
        }
        else {
            uploA = (uplo == ChamUpper) ? ChamLower : ChamUpper;
            Anorm = LAPACKE_zlantr( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uploA), 'N', Am, An, A, LDA );
        }
        Binitnorm = LAPACKE_zlantr( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), 'N', M, N, Bref, LDB );
    }

    cham_fixdbl_t eps = testing_getaccuracy();
    double *work = malloc(chameleon_max(M, N)* sizeof(double));

    /* Makes the sum with the core function */
    if ( uplo == ChamUpperLower ) {
        CORE_zgeadd( trans, M, N, alpha, A, LDA, beta, Bref, LDB );
    }
    else {
        CORE_ztradd( uplo, trans, M, N, alpha, A, LDA, beta, Bref, LDB );
    }
    cblas_zaxpy( LDB*N, CBLAS_SADDR(mzone), Bcham, 1, Bref, 1 );

    /* Calculates the norm from the core function's result */
    if ( uplo == ChamUpperLower ) {
        Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Bref, LDB, work );
    }
    else {
        Rnorm = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), 'N',
                                        M, N, Bref, LDB, work );
    }
    result = Rnorm / (max(Anorm, Binitnorm) * eps);

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||B||", Binitnorm );
    run_arg_add_double( args, "||R||", Rnorm );

    /* Verifies if the result is inside a threshold */
    if ( isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    free(work);
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed sum with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] trans
 *          Whether the first matrix is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descBref
 *          The descriptor of the matrix Bref.
 *
 * @param[in] descBcham
 *          The descriptor of the matrix of the Chameleon computed result A+B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsum ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                 CHAMELEON_Complex64_t beta, CHAM_desc_t *descBref, CHAM_desc_t *descBcham )
{
    int info_solution            = 0;
    int M                        = descBref->m;
    int N                        = descBref->n;
    int Am                       = (trans == ChamNoTrans)? M : N;
    int An                       = (trans == ChamNoTrans)? N : M;
    int LDA                      = Am;
    int LDB                      = M;
    int rank                     = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A     = NULL;
    CHAMELEON_Complex64_t *Bref  = NULL;
    CHAMELEON_Complex64_t *Bcham = NULL;
    cham_uplo_t uploA            = uplo;

    if ( rank == 0 ) {
        A     = malloc( LDA*An*sizeof(CHAMELEON_Complex64_t) );
        Bref  = malloc( LDB*N* sizeof(CHAMELEON_Complex64_t) );
        Bcham = malloc( LDB*N* sizeof(CHAMELEON_Complex64_t) );
    }

    if ( uplo != ChamUpperLower && trans != ChamNoTrans ) {
        uploA = (uplo == ChamUpper) ? ChamLower : ChamUpper;
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uploA, descA,     A,     LDA );
    CHAMELEON_zDesc2Lap( uplo,  descBref,  Bref,  LDB );
    CHAMELEON_zDesc2Lap( uplo,  descBcham, Bcham, LDB );

    if ( rank == 0 ) {
        info_solution = check_zsum_std( args, uplo, trans, M, N, alpha, A, LDA, beta, Bref, Bcham, LDB );

        free(A);
        free(Bref);
        free(Bcham);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core functions computed scale.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] M
 *          The number of rows of the matrices A and Ainit.
 *
 * @param[in] N
 *          The number of columns of the matrices A and Ainit.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] Ainit
 *          The matrix Ainit.
 *
 * @param[in] A
 *          The scaled matrix A.
 *
 * @param[in] LDA
 *      The leading dimension of the matrices A and Ainit.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zscale_std( run_arg_list_t *args, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *Ainit, CHAMELEON_Complex64_t *A, int LDA )
{
    int info_solution;

    /* Scales using core function */
    CORE_zlascal( uplo, M, N, alpha, Ainit, LDA );

    /* Compares the two matrices */
    info_solution = check_zmatrices_std( args, uplo, M, N, Ainit, LDA, A, LDA );

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed scale with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descAinit
 *          The descriptor of the matrix Ainit.
 *
 * @param[in] descA
 *          The descriptor of the scaled matrix A.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zscale( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descAinit, CHAM_desc_t *descA )
{
    int info_solution;
    int M                        = descA->m;
    int N                        = descA->n;
    int LDA                      = M;
    int rank                     = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A     = NULL;
    CHAMELEON_Complex64_t *Ainit = NULL;

    if ( rank == 0 ) {
        A     = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
        Ainit = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
    }

    /* Converts the matrix to LAPACK layout in order to scale with BLAS */
    CHAMELEON_zDesc2Lap( uplo, descA,     A,     LDA );
    CHAMELEON_zDesc2Lap( uplo, descAinit, Ainit, LDA );

    /* Compares the two matrices */
    if ( rank == 0 ) {
        info_solution = check_zscale_std( args, uplo, M, N, alpha, Ainit, A, LDA );
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

    if ( rank == 0 ) {
        free( A );
        free( Ainit );
    }

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the rank of the matrix A is K.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[in] K
 *          The rank of the matrix A.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @retval 0 success, else failure
 *
 *******************************************************************************
 */
int check_zrankk_std( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA )
{
    int info_solution = 0;
    int minMN = chameleon_min(M, N);
    double Anorm, Rnorm, result;
    cham_fixdbl_t eps = testing_getaccuracy();

    Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'F', M, N, A, LDA );

    /* check rank of A using SVD, value K+1 of Sigma must be small enough */
    CHAMELEON_Complex64_t *U  = malloc( M * M * sizeof(CHAMELEON_Complex64_t) );
    CHAMELEON_Complex64_t *VT = malloc( N * N * sizeof(CHAMELEON_Complex64_t) );
    double *S    = malloc( minMN * sizeof(double) );
    double *work = malloc( minMN * sizeof(double) );

    LAPACKE_zgesvd( LAPACK_COL_MAJOR, 'A', 'A', M, N, A, LDA, S, U, M, VT, N, work );

    /* Computes the residual's norm */
    if ( K >= minMN ) {
        Rnorm = 0.;
    } else {
        Rnorm = S[K];
    }

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||R||", Rnorm );

    result = Rnorm / (Anorm * eps);

    /* Verifies if the result is inside a threshold */
    if ( isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    free(S);
    free(U);
    free(VT);
    free(work);

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the rank of the matrix A is K.
 *
 *******************************************************************************
 *
 * @param[in] K
 *          The rank of the matrix A.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @retval 0 success, else failure
 *
 *******************************************************************************
 */
int check_zrankk( run_arg_list_t *args, int K, CHAM_desc_t *descA )
{
    int info_solution = 0;
    int M = descA->m;
    int N = descA->n;
    int LDA = descA->m;
    int rank = CHAMELEON_Comm_rank();

    /* Converts the matrices to LAPACK layout in order to check values on the main process */
    CHAMELEON_Complex64_t *A = NULL;
    if ( rank == 0 ) {
        A = malloc( M*N*sizeof(CHAMELEON_Complex64_t) );
    }
    CHAMELEON_Desc2Lap( ChamUpperLower, descA, A, LDA );

    /* check rank of A using SVD, value K+1 of Sigma must be small enough */
    if ( rank == 0 ) {
        info_solution = check_zrankk_std( args, M, N, K, A, LDA );
        free( A );
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

    return info_solution;
}

#endif /* defined(CHAMELEON_SIMULATION) */
