/**
 *
 * @file testing_zcheck.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2022-02-22
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
    double eps               = LAPACKE_dlamch_work('e');

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
        A = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
        B = (CHAMELEON_Complex64_t *)malloc(LDB*N*sizeof(CHAMELEON_Complex64_t));
    }

    /* Converts the matrices to LAPACK layout in order to compare them on the main process */
    CHAMELEON_zDesc2Lap( uplo, descA, A, LDA );
    CHAMELEON_zDesc2Lap( uplo, descB, B, LDB );

    if ( rank == 0 ) {
        info_solution = check_zmatrices_std( args, uplo, M, N, A, LDA, B, LDB );
        free(A);
        free(B);
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
    double eps         = LAPACKE_dlamch_work('e');

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

    double  eps  = LAPACKE_dlamch_work('e');
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
 * @brief Compares two core function computed products.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] transB
 *          Whether the second product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] M
 *          The number of rows of the matrices A and C.
 *
 * @param[in] N
 *          The number of columns of the matrices B and C.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the  umber of rows of the matrxi B.
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
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] Cref
 *          The matrix Cref.
 *
 * @param[in] C
 *          The matrix of the computed result alpha*A*B+beta*C.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices C and Cref.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgemm_std( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA,
                     CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC )
{
    int info_solution           = 0;
    double Anorm, Bnorm, Crefnorm, Rnorm, result;
    CHAMELEON_Complex64_t mzone = -1.0;

    /* Calculates the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', M, K, A, LDA );
    } else {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', K, M, A, LDA );
    }
    if ( transB == ChamNoTrans ) {
        Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', K, N, B, LDB );
    } else {
        Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', N, K, B, LDB );
    }

    /* Computes the norms for comparing */
    Crefnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );

    double eps = LAPACKE_dlamch_work('e');

    /* Makes the multiplication with the core function */
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K,
                 CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
    cblas_zaxpy( LDC * N, CBLAS_SADDR(mzone), C, 1, Cref, 1 );

    /* Calculates the norm with the core function's result */
    Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );

    if ( ( alpha != 0. ) || (beta != 0. ) ) {
        result = Rnorm / ( ( cabs(alpha) * max(Anorm, Bnorm) + cabs(beta) * Crefnorm ) * K * eps );
    }
    else {
        result = Rnorm;
    }
    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||B||", Bnorm );
    run_arg_add_double( args, "||C||", Crefnorm );
    run_arg_add_double( args, "||R||", Rnorm );

    /* Verifies if the result is inside a threshold */
    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] transB
 *          Whether the second product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the matrix Cref.
 *
 * @param[in] descC
 *          The descriptor of the matrix of the Chameleon computed result alpha*A*B+beta*C.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgemm( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                 CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int info_solution           = 0;
    int M                       = descC->m;
    int N                       = descC->n;
    int K                       = (transA != ChamNoTrans)? descA->m : descA->n;
    int rank                    = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *C    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    int An                      = ( transA == ChamNoTrans ) ? K : M;
    int LDA                     = ( transA == ChamNoTrans ) ? M : K;
    int Bn                      = ( transB == ChamNoTrans ) ? N : K;
    int LDB                     = ( transB == ChamNoTrans ) ? K : N;
    int LDC                     = M;

    /* Creates the LAPACK version of the matrices */
    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc(LDA*An*sizeof(CHAMELEON_Complex64_t));
        B    = (CHAMELEON_Complex64_t *)malloc(LDB*Bn*sizeof(CHAMELEON_Complex64_t));
        Cref = (CHAMELEON_Complex64_t *)malloc(LDC*N *sizeof(CHAMELEON_Complex64_t));
        C    = (CHAMELEON_Complex64_t *)malloc(LDC*N *sizeof(CHAMELEON_Complex64_t));
    }

    CHAMELEON_zDesc2Lap( ChamUpperLower, descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descC,    C,    LDC );

    if ( rank == 0 ) {

        info_solution = check_zgemm_std( args, transA, transB, alpha, M, N, K, A, LDA, B, LDB, beta, Cref, C, LDC );

        free(A);
        free(B);
        free(C);
        free(Cref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed symmetric products.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] side
 *          Whether the symmetric matrix A appears on the left or right in the operation.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] M
 *          The order of the matrix A and number of rows of the matrices B and C.
 *
 * @param[in] N
 *          The number of columns of the matrices B and C.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The symmetric matrix A.
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
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] Cref
 *          The matrix Cref.
 *
 * @param[in] C
 *          The matrix of the computed result alpha*A*B+beta*C or alpha*B*A+beta*C.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices C and Cref.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsymm_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_side_t side, cham_uplo_t uplo, int M, int N,
                     CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB,
                     CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC )
{
    int info_solution = 0;
    int An;
    double Anorm, Bnorm, Crefnorm, Cchamnorm, Clapacknorm, Rnorm, result;
    CHAMELEON_Complex64_t mzone = -1.0;
    double *work;
    char normTypeA, normTypeB;

    if ( side == ChamLeft ) {
        normTypeA = 'I';
        normTypeB = 'O';
        An = M;
    }
    else {
        normTypeA = 'O';
        normTypeB = 'I';
        An = N;
    }

    work = malloc( sizeof(double) * An );
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        Anorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, normTypeA, chameleon_lapack_const(uplo), An, A, LDA, work );
    }
    else
#endif
    {
        Anorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, normTypeA, chameleon_lapack_const(uplo), An, A, LDA, work );
    }
    Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, normTypeB, M, N, B, LDB );
    free( work );

    /* Computes the norms for comparing */
    Crefnorm  = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );
    Cchamnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, C,    LDC, NULL );

    double eps = LAPACKE_dlamch_work('e');

    /* Makes the multiplication with the core function */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        cblas_zhemm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                        M, N, CBLAS_SADDR(alpha),
                        A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
    }
    else
#endif
    {
        cblas_zsymm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB,
                     CBLAS_SADDR(beta), Cref, LDC );
    }
    cblas_zaxpy(LDC * N, CBLAS_SADDR(mzone), C, 1, Cref, 1);

    /* Computes the norm with the core function's result */
    Clapacknorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );
    Rnorm       = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );

    if ( ( alpha != 0. ) || (beta != 0. ) ) {
        result = Rnorm / ((cabs(alpha) * max(Anorm, Bnorm) + cabs(beta) * Crefnorm) * An * eps);
    }
    else {
        result = Rnorm;
    }

    /* Verifies if the result is inside a threshold */
    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    (void)Clapacknorm;
    (void)Cchamnorm;
    (void)matrix_type;
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed symmetric product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] side
 *          Whether the symmetric matrix A appears on the left or right in the operation.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the symmetric matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the matrix Cref.
 *
 * @param[in] descC
 *          The descriptor of the matrix of the Chameleon computed result alpha*A*B+beta*C or alpha*B*A+beta*C.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsymm( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_side_t side, cham_uplo_t uplo,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                 CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int info_solution = 0;
    int M             = descC->m;
    int N             = descC->n;
    int rank          = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    CHAMELEON_Complex64_t *C    = NULL;
    int An                      = ( side == ChamLeft )? M : N;
    int LDA                     = An;
    int LDB                     = M;
    int LDC                     = M;

    if ( rank == 0 ) {
        A    = ( CHAMELEON_Complex64_t * ) malloc( LDA*An*sizeof( CHAMELEON_Complex64_t ) );
        B    = ( CHAMELEON_Complex64_t * ) malloc( LDB*N *sizeof( CHAMELEON_Complex64_t ) );
        Cref = ( CHAMELEON_Complex64_t * ) malloc( LDC*N *sizeof( CHAMELEON_Complex64_t ) );
        C    = ( CHAMELEON_Complex64_t * ) malloc( LDC*N *sizeof( CHAMELEON_Complex64_t ) );
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uplo,           descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descC,    C,    LDC );

    if ( rank == 0 ) {
        info_solution = check_zsymm_std( args, matrix_type, side, uplo,
                                         M, N, alpha, A, LDA, B, LDB,
                                         beta, Cref, C, LDC );
        free(A);
        free(B);
        free(C);
        free(Cref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed matrices rank k operations.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] trans
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] N
 *          The order of the matrix C and number of rows of the matrices A and B.
 *
 * @param[in] K
 *          The number of columns of the matrices A and B.
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
 * @param[in] B
 *          The matrix B - only used for her2k and syr2k.
 *
 * @param[in] LDB
 *          The leading dimension of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] Cref
 *          The symmetric matrix Cref.
 *
 * @param[in] C
 *          The symmetric matrix of the computed result.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices Cref and C.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsyrk_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha,
                     CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref,
                     CHAMELEON_Complex64_t *C, int LDC )
{
    int     info_solution = 0;
    double  Anorm, Bnorm, Crefnorm, Cchamnorm, Clapacknorm, Rnorm, result;
    double *work = malloc(sizeof(double)*N);

    Bnorm = 0.;
    if ( trans == ChamNoTrans ) {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', N, K, A, LDA );
        if ( B != NULL ) {
            Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', N, K, B, LDB );
        }
    }
    else {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', K, N, A, LDA );
        if ( B != NULL ) {
            Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', K, N, B, LDB );
        }
    }

    /* Computes the norms for comparing */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        Crefnorm  = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
        Cchamnorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, C,    LDC, work );
    }
    else
#endif
    {
        Crefnorm  = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
        Cchamnorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, C,    LDC, work );
    }

    double eps = LAPACKE_dlamch_work('e');
    double ABnorm;

    /* Makes the multiplication with the core function */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        if ( B == NULL ) {
            cblas_zherk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, creal(alpha), A, LDA, creal(beta), Cref, LDC );
            ABnorm = Anorm * Anorm;
        }
        else {
            cblas_zher2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, CBLAS_SADDR(alpha), A, LDA, B, LDB, creal(beta), Cref, LDC );
            ABnorm = 2. * Anorm * Bnorm;
        }

        Clapacknorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
    }
    else
#endif
    {
        if ( B == NULL ) {
            cblas_zsyrk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, CBLAS_SADDR(alpha), A, LDA, CBLAS_SADDR(beta), Cref, LDC );
            ABnorm = Anorm * Anorm;
        }
        else {
            cblas_zsyr2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
            ABnorm = 2. * Anorm * Bnorm;
        }

        Clapacknorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
    }

    CORE_ztradd( uplo, ChamNoTrans, N, N, -1., C, LDC, 1., Cref, LDC );

    /* Computes the norm with the core function's result */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        Rnorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), N, Cref, LDC, NULL );
    }
    else
#endif
    {
        Rnorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), N, Cref, LDC, NULL );
    }
    result = Rnorm / ((ABnorm + Crefnorm) * K * eps);
    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||B||", Bnorm );
    run_arg_add_double( args, "||C||", Crefnorm );
    run_arg_add_double( args, "||R||", Rnorm );
    /* Verifies if the result is inside a threshold */
    if ( isinf(Clapacknorm) || isinf(Cchamnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    free(work);

    (void)matrix_type;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed matrix rank k operation with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] trans
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B - only used for her2k and syr2k.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the symmetric matrix C.
 *
 * @param[in] descC
 *          The descriptor of the symmetric matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsyrk( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_trans_t trans,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                 CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int info_solution           = 0;
    int An, Bn, K, N;
    int LDA, LDB, LDC;
    int rank                    = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    CHAMELEON_Complex64_t *C    = NULL;

    if ( trans == ChamNoTrans ) {
        N   = descA->m;
        K   = descA->n;
        An  = K;
        Bn  = K;
        LDA = N;
        LDB = N;
    }
    else {
        N   = descA->n;
        K   = descA->m;
        An  = N;
        Bn  = N;
        LDA = K;
        LDB = K;
    }
    LDC = N;

    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc( LDA * An * sizeof(CHAMELEON_Complex64_t) );
        if ( descB != NULL ) {
            B = (CHAMELEON_Complex64_t *)malloc( LDB * Bn * sizeof(CHAMELEON_Complex64_t) );
        }
        Cref = (CHAMELEON_Complex64_t *)malloc( LDC * N * sizeof(CHAMELEON_Complex64_t) );
        C    = (CHAMELEON_Complex64_t *)malloc( LDC * N * sizeof(CHAMELEON_Complex64_t) );
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( ChamUpperLower, descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( uplo,           descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( uplo,           descC,    C,    LDC );
    if ( descB != NULL ) {
        CHAMELEON_zDesc2Lap( ChamUpperLower, descB, B, LDB );
    }

    if ( rank == 0 ) {

        info_solution = check_zsyrk_std( args, matrix_type, uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, Cref, C, LDC );

        free(A);
        free(C);
        free(Cref);
        if ( descB != NULL ) {
            free(B);
        }
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed matrix triangular product.
 *
 *******************************************************************************
 *
 * @param[in] check_func
 *          Whether it is a triangular product or a triangular linear solution.
 *
 * @param[in] side
 *          Whether A appears on the left or on the right of the product.
 *
 * @param[in] uplo
 *          Whether A is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] trans
 *          Whether A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] diag
 *          Whether A is a unitary diagonal matrix or not.
 *
 * @param[in] M
 *          The number of rows of the matrix B and the order of the matrix A if side is left.
 *
 * @param[in] N
 *          The number of columns of the matrix B and the order of the matrix A if side is right.
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
 * @param[in] Bref
 *          The matrix Bref.
 *
 * @param[in] B
 *          The matrix of the computed result.
 *
 * @param[in] LDB
 *          The leading dimension of the matrices Bref and B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrmm_std( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag, int M, int N,
                     CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *B, int LDB )
{
    CHAMELEON_Complex64_t mzone = -1.0;
    int     info_solution, An;
    double  Anorm, Bnorm, Rnorm, result;
    char    normTypeA, normTypeB;
    double *work;
    double  eps = LAPACKE_dlamch_work('e');

    /* Computes the norms for comparing */
    if ( side == ChamLeft ) {
        normTypeA = 'O';
        if ( trans == ChamNoTrans ) {
            normTypeA = 'I';
        }
        normTypeB = 'O';
        An = M;
    }
    else {
        normTypeA = 'O';
        if ( trans != ChamNoTrans ) {
            normTypeA = 'I';
        }
        normTypeB = 'I';
        An = N;
    }

    work  = malloc( sizeof(double) * An );
    Anorm = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, normTypeA,
                                 chameleon_lapack_const(uplo), chameleon_lapack_const(diag),
                                 An, An, A, LDA, work );
    free( work );

    Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, normTypeB, M, N, B, LDB );

    /* Makes the multiplication with the core function */
    if (check_func == CHECK_TRMM) {
        cblas_ztrmm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                    (CBLAS_DIAG)diag, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB);
    }
    else {
        cblas_ztrsm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                    (CBLAS_DIAG)diag, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB);
    }

    /* Computes the norm with the core function's result */
    cblas_zaxpy( LDB * N, CBLAS_SADDR(mzone), Bref, 1, B, 1 );
    Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, B, LDB, NULL );

    result = Rnorm / ((Anorm + Bnorm) * An * eps);

    /* Verifies if the result is inside a threshold */
    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    (void)Bnorm;
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed matrix triangular product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] check_func
 *          Whether it is a triangular product or a triangular linear solution.
 *
 * @param[in] side
 *          Whether A appears on the left or on the right of the product.
 *
 * @param[in] uplo
 *          Whether A is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] trans
 *          Whether A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] diag
 *          Whether A is a unitary diagonal matrix or not.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descBref
 *          The descriptor of the matrix Bref.
 *
 * @param[in] descB
 *          The descriptor of the matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrmm( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descBref, CHAM_desc_t *descB )
{
    int info_solution           = 0;
    int M                       = descBref->m;
    int N                       = descBref->n;
    int rank                    = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *Bref = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    int An                      = ( side == ChamLeft )? M : N;
    int LDA                     = An;
    int LDB                     = M;

    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc( LDA*An*sizeof(CHAMELEON_Complex64_t) );
        Bref = (CHAMELEON_Complex64_t *)malloc( LDB*N *sizeof(CHAMELEON_Complex64_t) );
        B    = (CHAMELEON_Complex64_t *)malloc( LDB*N *sizeof(CHAMELEON_Complex64_t) );
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uplo,           descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descBref, Bref, LDB );

    if ( rank == 0 ) {

        info_solution = check_ztrmm_std( args, check_func, side, uplo, trans, diag, M, N, alpha, A, LDA, Bref, B, LDB );

        free(A);
        free(B);
        free(Bref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed product U*U' or L'*L result with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether the upper or lower triangle of A is stored.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descAAt
 *          The descriptor of the matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zlauum( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descAAt )
{
    int info_local, info_global;
    int N       = descA->n;
    double eps  = LAPACKE_dlamch_work('e');
    double result, Anorm, AAtnorm, Rnorm;
    CHAM_desc_t *descAt;

    Anorm   = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, ChamNonUnit, descA );
    AAtnorm = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, ChamNonUnit, descAAt );

    if ( uplo == ChamUpper ) {
        descAt = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zlaset_Tile( ChamLower, 0., 0., descAt );
        CHAMELEON_zlacpy_Tile( ChamUpper, descA, descAt );

        /* Computes U * U' */
        CHAMELEON_ztrmm_Tile( ChamRight, ChamUpper, ChamConjTrans, ChamNonUnit, 1., descA, descAt );
    }
    else {
        descAt = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zlaset_Tile( ChamUpper, 0., 0., descAt );
        CHAMELEON_zlacpy_Tile( ChamLower, descA, descAt );

        /* Computes L' * L */
        CHAMELEON_ztrmm_Tile( ChamLeft, ChamLower, ChamConjTrans, ChamNonUnit, 1., descA, descAt );
    }

    /* Computes AAt - A * A' */
    CHAMELEON_ztradd_Tile( uplo, ChamNoTrans, -1., descAAt, 1., descAt );

    Rnorm = CHAMELEON_zlantr_Tile( ChamMaxNorm, uplo, ChamNonUnit, descAt );

    CHAMELEON_Desc_Destroy( &descAt );

    /* Compares the residual's norm */
    result = Rnorm / ( Anorm * Anorm * N * eps );
    if (  isnan(AAtnorm) || isinf(AAtnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two product U*U' or L'*L result.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether the upper or lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrices A and AAt.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] AAt
 *          The matrix of the computed result.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A and AAt.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zlauum_std( run_arg_list_t *args, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AAt, int LDA )
{
    int          info;
    CHAM_desc_t *descA;
    CHAM_desc_t *descAAt;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,   CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAAt, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );

    CHAMELEON_zLap2Desc( uplo, A,   LDA, descA   );
    CHAMELEON_zLap2Desc( uplo, AAt, LDA, descAAt );

    info = check_zlauum( args, uplo, descA, descAAt );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descAAt );

    (void)args;
    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a Chameleon computed factorization is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] descA
 *          The descriptor of the symmetric matrix A.
 *
 * @param[in] descLU
 *          The descriptor of the matrix of the Chameleon factorisation of the matrix A.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxtrf( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descLU )
{
    int info_local, info_global;
    int M      = descA->m;
    int N      = descA->n;
    double Anorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');

    CHAM_desc_t *descL, *descU;
    cham_trans_t transL = ChamNoTrans;
    cham_trans_t transU = ChamNoTrans;

    descL = CHAMELEON_Desc_Copy( descA, NULL );
    descU = CHAMELEON_Desc_Copy( descA, NULL );

    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descL );
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descU );

    switch ( uplo ) {
    case ChamUpper:
#if defined(PRECISION_z) || defined(PRECISION_c)
        transL = ( matrix_type == ChamHermitian ) ? ChamConjTrans : ChamTrans;
#else
        transL = ChamTrans;
#endif
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descL );
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descU );
        break;
    case ChamLower:
#if defined(PRECISION_z) || defined(PRECISION_c)
        transU = ( matrix_type == ChamHermitian ) ? ChamConjTrans : ChamTrans;
#else
        transU = ChamTrans;
#endif
        CHAMELEON_zlacpy_Tile( ChamLower, descLU, descL );
        CHAMELEON_zlacpy_Tile( ChamLower, descLU, descU );
        break;
    case ChamUpperLower:
    default:
        CHAMELEON_zlacpy_Tile( ChamLower, descLU, descL );
        CHAMELEON_zlaset_Tile( ChamUpper, 0., 1., descL );
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descU );
    }

    switch ( matrix_type ) {
    case ChamGeneral: {
        CHAM_desc_t *subL, *subU;
        subL = chameleon_desc_submatrix( descL, 0, 0, M, chameleon_min(M, N) );
        subU = chameleon_desc_submatrix( descU, 0, 0, chameleon_min(M, N), N );

        Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
        CHAMELEON_zgemm_Tile( transL, transU, -1., subL, subU, 1., descA );
        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );

        free( subL );
        free( subU );
    }
        break;

#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        Anorm = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descA );
        CHAMELEON_zgemm_Tile( transL, transU, -1., descL, descU, 1., descA );
        Rnorm = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descA );
        break;
#endif

    case ChamSymmetric:
        Anorm = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descA );
        CHAMELEON_zgemm_Tile( transL, transU, -1., descL, descU, 1., descA );
        Rnorm = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descA );
        break;

    default:
        fprintf(stderr, "check_zxxtrf: matrix_type(%d) unsupported\n", matrix_type );
        return 1;
    }

    result = Rnorm / ( Anorm * N * eps );
    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else{
        info_local = 0;
    }

    CHAMELEON_Desc_Destroy( &descL );
    CHAMELEON_Desc_Destroy( &descU );

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a core function computed factorization is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] M
 *          The number of rows of the matrices A and LU.
 *
 * @param[in] N
 *          The number of columns of the matrices A and LU.
 *
 * @param[in] A
 *          The symmetric matrix A.
 *
 * @param[in] LU
 *          The matrix with the computed factorisation of the matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A and LU.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxtrf_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, int M, int N,
                      CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *LU, int LDA )
{
    int          info;
    CHAM_desc_t *descA;
    CHAM_desc_t *descLU;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descLU, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A,  LDA, descA  );
    CHAMELEON_zLap2Desc( ChamUpperLower, LU, LDA, descLU );

    info = check_zxxtrf( args, matrix_type, uplo, descA, descLU );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descLU );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the  linear solution of op(A) * x = b is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] trans
 *          Whether the A matrix is non transposed, tranposed or conjugate transposed.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] descA
 *          The descriptor of the symmetric matrix A.
 *
 * @param[in] descX
 *          The descriptor of the matrix X.
 *
 * @param[inout] descB
 *          The descriptor of the matrix B = A*X. On exit, it contains the remainder from A*x-B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsolve( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_trans_t trans, cham_uplo_t uplo,
                  CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB )
{
    int info_local, info_global;
    int M                = descA->m;
    int N                = descA->n;
    double Anorm, Bnorm, Xnorm, Rnorm, result = 0;
    double eps           = LAPACKE_dlamch_work('e');
    cham_normtype_t norm = (trans == ChamNoTrans) ? ChamOneNorm : ChamInfNorm;

    /* Computes the norms */
    Bnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
    Xnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descX );

    switch ( matrix_type ) {
    case ChamGeneral:
        Anorm = CHAMELEON_zlange_Tile( norm, descA );
        CHAMELEON_zgemm_Tile( trans, ChamNoTrans, -1., descA, descX, 1., descB );
        break;

#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        Anorm = CHAMELEON_zlanhe_Tile( norm, uplo, descA );
        CHAMELEON_zhemm_Tile( ChamLeft, uplo, -1., descA, descX, 1., descB );
        break;
#endif

    case ChamSymmetric:
        Anorm = CHAMELEON_zlansy_Tile( norm, uplo, descA );
        CHAMELEON_zsymm_Tile( ChamLeft, uplo, -1., descA, descX, 1., descB );
        break;

    default:
        fprintf(stderr, "check_zsolve: matrix_type(%d) unsupported\n", matrix_type );
        return 1;
    }

    Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
    result = Rnorm / ( Anorm * Xnorm * chameleon_max( M, N ) * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||X||", Xnorm );
    run_arg_add_double( args, "||B||", Bnorm );
    run_arg_add_double( args, "||Ax-b||", Rnorm );
    run_arg_add_double( args, "||Ax-b||/N/eps/(||A||||x||+||b||", result );

    if ( isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the  linear solution of op(A) * x = b is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] trans
 *          Whether the A matrix is non transposed, tranposed or conjugate transposed.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] N
 *          The order of the matrix A and the number of rows of the matrices X and B.
 *
 * @param[in] NRHS
 *          The number of columns of the matrices X and B.
 *
 * @param[in] A
 *          The symmetric matrix A.
 *
 * @param[in] LDA
 *          The leading dimenson of the matrix A.
 *
 * @param[in] X
 *          The matrix X.
 *
 * @param[inout] B
 *          The matrix B = A*X. On exit, it contains the remainder from A*x-B.
 *
 * @param[in] LDB
 *          The leading dimension of the matrices X and B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsolve_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_trans_t trans, cham_uplo_t uplo, int N, int NRHS,
                      CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, CHAMELEON_Complex64_t *B, int LDB )
{
    int          info;
    CHAM_desc_t *descA, *descX, *descB;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N,    0, 0, N, N,    1, 1 );
    CHAMELEON_Desc_Create(
        &descB, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, 1, 1 );
    CHAMELEON_Desc_Create(
        &descX, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, 1, 1 );

    CHAMELEON_zLap2Desc( uplo,           A, LDA, descA );
    CHAMELEON_zLap2Desc( ChamUpperLower, B, LDB, descB );
    CHAMELEON_zLap2Desc( ChamUpperLower, X, LDB, descX );

    info = check_zsolve( args, matrix_type, trans, uplo, descA, descX, descB );

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descX );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the matrix A0 is the inverse of Ai.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether they are upper triangular matrices or lower triangular matrices.
 *
 * @param[in] diag
 *          Whether they are unitary diagonal matrices or not.
 *
 * @param[in] descA0
 *          The descriptor of the matrix A0.
 *
 * @param[in] descAi
 *          The descriptor of the matrix Ai.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrtri( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_diag_t diag,
                  CHAM_desc_t *descA0, CHAM_desc_t *descAi )
{
    int          info_local, info_global;
    cham_uplo_t  uplo_inv;
    CHAM_desc_t *descI, *descB = NULL;
    double       Rnorm, Anorm, Ainvnorm, result;
    double       eps = LAPACKE_dlamch_work('e');
    int          N   = descA0->m;

    /* Creates an identity matrix */
    descI = CHAMELEON_Desc_Copy( descA0, NULL );
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., descI );

    /* Calculates the residual I - A*(A**-1) */
    switch ( matrix_type ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        /* Ainv comes from potri and is hermitian */
        assert( uplo != ChamUpperLower );

        Anorm    = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descA0 );
        Ainvnorm = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descAi );

        /*
         * Expand Ainv into a full matrix and call ZHEMM to multiply
         * Ainv on the left by A.
         */
        uplo_inv = ( uplo == ChamUpper ) ? ChamLower : ChamUpper;
        descB = CHAMELEON_Desc_Copy( descAi, NULL );
        CHAMELEON_ztradd_Tile( uplo_inv, ChamConjTrans, 1., descAi, 0., descB );
        CHAMELEON_zlacpy_Tile( uplo, descAi, descB );

        CHAMELEON_zhemm_Tile( ChamLeft, uplo, -1., descA0, descB, 1., descI );
        break;
#endif

    case ChamSymmetric:
        /* Ainv comes from potri and is symmetric */
        assert( uplo != ChamUpperLower );

        Anorm    = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descA0 );
        Ainvnorm = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descAi );

        /*
         * Expand Ainv into a full matrix and call ZHEMM to multiply
         * Ainv on the left by A.
         */
        uplo_inv = ( uplo == ChamUpper ) ? ChamLower : ChamUpper;
        descB = CHAMELEON_Desc_Copy( descAi, NULL );
        CHAMELEON_ztradd_Tile( uplo_inv, ChamTrans, 1., descAi, 0., descB );
        CHAMELEON_zlacpy_Tile( uplo, descAi, descB );

        CHAMELEON_zsymm_Tile( ChamLeft, uplo, -1., descA0, descB, 1., descI );
        break;

    case ChamTriangular:
        /* Ainv comes from trtri */
        assert( uplo != ChamUpperLower );

        Anorm    = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, diag, descA0 );
        Ainvnorm = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, diag, descAi );

        /*
         * Expand Ainv into a full matrix and call ZHEMM to multiply
         * Ainv on the left by A.
         */
        uplo_inv = ( uplo == ChamUpper ) ? ChamLower : ChamUpper;
        descB = CHAMELEON_Desc_Copy( descAi, NULL );

        if ( diag == ChamUnit ) {
            //CHAMELEON_ztradd_Tile( uplo, ChamNoTrans, 1., descAi, 0., descB );
            CHAMELEON_zlacpy_Tile( uplo, descAi, descB );
            CHAMELEON_zlaset_Tile( uplo_inv, 0., 1., descB );
        }
        else {
            CHAMELEON_zlaset_Tile( uplo_inv, 0., 1., descB );
            CHAMELEON_zlacpy_Tile( uplo, descAi, descB );
            //CHAMELEON_ztradd_Tile( uplo, ChamNoTrans, 1., descAi, 0., descB );
        }

        /* Computes - A * A^-1 */
        CHAMELEON_ztrmm_Tile( ChamLeft, uplo, ChamNoTrans, diag, -1., descA0, descB );
        /* Computes I - A * A^-1 */
        CHAMELEON_zgeadd_Tile( ChamNoTrans, 1., descB, 1., descI );
        break;

    case ChamGeneral:
    default:
        /* Ainv comes from getri */
        assert( uplo == ChamUpperLower );

        Anorm    = CHAMELEON_zlange_Tile( ChamOneNorm, descA0 );
        Ainvnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descAi );

        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, -1., descA0, descAi, 1., descI );
        break;
    }

    Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descI );

    /* Compares the residual's norm */
    result = Rnorm / ( Anorm * Ainvnorm * N * eps );
    if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    CHAMELEON_Desc_Destroy( &descI );
    if ( descB != NULL ) {
        CHAMELEON_Desc_Destroy( &descB );
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the matrix A0 is the inverse of Ai.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether they are upper triangular matrices or lower triangular matrices.
 *
 * @param[in] diag
 *          Whether they are unitary diagonal matrices or not.
 *
 * @param[in] N
 *          the order of the matrices A0 and Ai.
 *
 * @param[in] A0
 *          The matrix A0.
 *
 * @param[in] Ai
 *          The matrix Ai.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A0 and Ai.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrtri_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_diag_t diag,
                      int N, CHAMELEON_Complex64_t *A0, CHAMELEON_Complex64_t *Ai, int LDA )
{
    int          info;
    CHAM_desc_t *descA0, *descAi;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA0, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAi, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );

    CHAMELEON_zLap2Desc( uplo, A0, LDA, descA0 );
    CHAMELEON_zLap2Desc( uplo, Ai, LDA, descAi );

    info = check_ztrtri( args, matrix_type, uplo, diag, descA0, descAi );

    CHAMELEON_Desc_Destroy( &descA0 );
    CHAMELEON_Desc_Destroy( &descAi );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a matrix is orthogonal.
 *
 *******************************************************************************
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zortho( run_arg_list_t *args, CHAM_desc_t *descQ )
{
    int info_local, info_global;
    int M      = descQ->m;
    int N      = descQ->n;
    int minMN  = chameleon_min(M, N);
    double result, normR;
    double eps = LAPACKE_dlamch_work('e');
    CHAM_desc_t *descI, *subI;

    /* Builds the identity */
    descI = CHAMELEON_Desc_Copy( descQ, NULL );
    subI = chameleon_desc_submatrix( descI, 0, 0, minMN, minMN );
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., subI );

    /* Performs Id - Q'Q */
    if ( M >= N ) {
        CHAMELEON_zherk_Tile( ChamUpper, ChamConjTrans, -1., descQ, 1., subI );
    }
    else {
        CHAMELEON_zherk_Tile( ChamUpper, ChamNoTrans, -1., descQ, 1., subI );
    }

    /* Verifies the residual's norm */
    normR = CHAMELEON_zlansy_Tile( ChamOneNorm, ChamUpper, subI );
    result = normR / ( (double)minMN * eps );

    run_arg_add_double( args, "||I-QQ'||", normR );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    free( subI );
    CHAMELEON_Desc_Destroy( &descI );

    /* Reduces the result on all processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a matrix is orthogonal.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix Q.
 *
 * @param[in] N
 *          The number of columns of the matrix Q.
 *
 * @param[in] descQ
 *          The matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zortho_std( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int          info;
    CHAM_desc_t *descQ;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descQ, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, N, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, Q, LDQ, descQ );

    info = check_zortho( args, descQ );

    CHAMELEON_Desc_Destroy( &descQ );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q generated with a call to ungqr and
 *          the computed factorisation of the matrix A (descAF).
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgelqf( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ )
{
    int info_local, info_global;
    int M = descQ->m;
    int N = descQ->n;
    int K = chameleon_min( descA->m, descA->n );
    double result, Anorm, Rnorm;
    double eps = LAPACKE_dlamch_work('e');
    CHAM_desc_t *descL;
    int full_lq = ( M == N ) ? 1 : 0;

    assert( descA->n == N );
    assert( descA->m == descAF->m );

    descL = CHAMELEON_Desc_Copy( descA, NULL );

    if ( full_lq ) {
        /*
         * Cas lapack zlqt01.f
         * A full N-by-N Q has been generated
         */
        assert( descAF->n == N );

        /* Copy L */
        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descL );
        CHAMELEON_zlacpy_Tile( ChamLower, descAF, descL );

        /* Compute L - A * Q' */
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamConjTrans, -1., descA, descQ, 1., descL );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descL );
    }
    else {
        /*
         * Cas lapack zlqt02.f
         * A partial Q has been generated (K < min(M, N))
         */
        CHAM_desc_t *subL, *subAF;

        assert( descAF->n >= M );
        assert( N >= M );

        /* Copy L(1:k,1:m) */
        subL  = chameleon_desc_submatrix( descL,  0, 0, K, M );
        subAF = chameleon_desc_submatrix( descAF, 0, 0, K, M );

        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., subL );
        CHAMELEON_zlacpy_Tile( ChamLower, subAF, subL );

        /* Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)' */
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamConjTrans, -1., descA, descQ, 1., subL );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, subL );

        free( subL );
        free( subAF );
    }

    CHAMELEON_Desc_Destroy(&descL);

    Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    result = Rnorm / ( (double)N * Anorm * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Reduces the result on all processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
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
 *          The number of elementary reflectors whose product defines the matrix Q.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] Q
 *          The matrix Q generated with a call to ungqr and the computed factorisation of the matrix A (descAF).
 *
 * @param[in] LDQ
 *          The leading dimension fo the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgelqf_std( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int          info;
    CHAM_desc_t *descA, *descAF, *descQ;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, K, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAF, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, K, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descQ,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, N, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A,  LDA, descA  );
    CHAMELEON_zLap2Desc( ChamUpperLower, AF, LDA, descAF );
    CHAMELEON_zLap2Desc( ChamUpperLower, Q,  LDQ, descQ  );

    info = check_zgelqf( args, descA, descAF, descQ );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descAF );
    CHAMELEON_Desc_Destroy( &descQ  );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q generated with a call to ungqr and
 *          the computed factorisation of the matrix A (descAF).
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgeqrf( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ )
{
    int info_local, info_global;
    int M = descQ->m;
    int N = descQ->n;
    int K = chameleon_min( descA->m, descA->n );
    double result, Anorm, Rnorm;
    double eps = LAPACKE_dlamch_work('e');
    CHAM_desc_t *descR;
    int full_qr = ( M == N ) ? 1 : 0;

    assert( descA->m == M );
    assert( descA->n == descAF->n );

    descR = CHAMELEON_Desc_Copy( descA, NULL );

    if ( full_qr ) {
        /*
         * Cas lapack zqrt01.f
         * A full M-by-M Q has been generated
         */
        assert( descAF->m == M );

        /* Copy R */
        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descR );
        CHAMELEON_zlacpy_Tile( ChamUpper, descAF, descR );

        /* Compute R - Q'*A */
        CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ, descA, 1., descR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descR );
    }
    else {
        /*
         * Cas lapack zqrt02.f
         * A partial Q has been generated (K < min(M, N))
         */
        CHAM_desc_t *subR, *subAF;

        assert( descAF->m >= N );
        assert( N <= M );

        /* Copy R(1:n,1:k) */
        subR  = chameleon_desc_submatrix( descR,  0, 0, N, K );
        subAF = chameleon_desc_submatrix( descAF, 0, 0, N, K );

        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., subR );
        CHAMELEON_zlacpy_Tile( ChamUpper, subAF, subR );

        /* Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k) */
        CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ, descA, 1., subR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, subR );

        free( subR );
        free( subAF );
    }

    CHAMELEON_Desc_Destroy(&descR);

    Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    result = Rnorm / ( (double)M * Anorm * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Reduces the result on all processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
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
 *          The number of elementary reflectors whose product defines the matrix Q.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] Q
 *          The matrix Q generated with a call to ungqr and
 *          the computed factorisation of the matrix A (descAF).
 *
 * @param[in] LDQ
 *          The leading dimension fo the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgeqrf_std( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int          info;
    CHAM_desc_t *descA, *descAF, *descQ;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, M, K, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAF, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, M, K, 1, 1 );
    CHAMELEON_Desc_Create(
        &descQ,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, M, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A,  LDA, descA  );
    CHAMELEON_zLap2Desc( ChamUpperLower, AF, LDA, descAF );
    CHAMELEON_zLap2Desc( ChamUpperLower, Q,  LDQ, descQ  );

    info = check_zgeqrf( args, descA, descAF, descQ );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descAF );
    CHAMELEON_Desc_Destroy( &descQ  );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the decomposition is correct.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Whether the matrix Q appears on the left or on the right of the product.
 *
 * @param[in] trans
 *          Whether the matrix Q is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descC
 *          The descriptor of the matrix C.
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q.
 *
 * @param[in] descCC
 *          The descriptor of the matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zqc( run_arg_list_t *args, cham_side_t side, cham_trans_t trans,
               CHAM_desc_t *descC, CHAM_desc_t *descQ, CHAM_desc_t *descCC )
{
    int info_local, info_global;
    int M = descQ->m;
    double Cnorm, Qnorm, CCnorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');

    Cnorm  = CHAMELEON_zlange_Tile( ChamOneNorm, descC );
    Qnorm  = CHAMELEON_zlange_Tile( ChamOneNorm, descQ );
    CCnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descCC );

    if ( side == ChamLeft ) {
        CHAMELEON_zgemm_Tile( trans, ChamNoTrans, -1., descQ, descC, 1., descCC );
    }
    else {
        CHAMELEON_zgemm_Tile( ChamNoTrans, trans, -1., descC, descQ, 1., descCC );
    }

    Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descCC );
    result = Rnorm / ( M * Cnorm * eps );

    run_arg_add_double( args, "||A||", Qnorm );
    run_arg_add_double( args, "||B||", CCnorm );
    run_arg_add_double( args, "||C||", Cnorm );
    run_arg_add_double( args, "||R||", Rnorm );

    if ( isnan(CCnorm) || isinf(CCnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)Qnorm;
    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the decomposition is correct.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Whether the matrix Q appears on the left or on the right of the product.
 *
 * @param[in] trans
 *          Whether the matrix Q is transposed, conjugate transposed or not transposed.
 *
 * @param[in] M
 *          The number of rows of the matrices C and CC and the order of the matrix Q if side = ChamLeft.
 *
 * @param[in] N
 *          The number of columns of the matrices C and CC and the order of the matrix Q if side = ChamRight.
 *
 * @param[in] C
 *          The matrix C.
 *
 * @param[in] CC
 *          The matrix of the computed result.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices C and CC.
 *
 * @param[in] Q
 *          The matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zqc_std( run_arg_list_t *args, cham_side_t side, cham_trans_t trans, int M, int N,
                   CHAMELEON_Complex64_t *C, CHAMELEON_Complex64_t *CC, int LDC, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int info;
    int Qm = ( side == ChamLeft ) ? M : N;
    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAM_desc_t *descC, *descCC, *descQ;

    CHAMELEON_Desc_Create(
        &descC,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N,  0, 0, M,  N,  1, 1 );
    CHAMELEON_Desc_Create(
        &descCC, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N,  0, 0, M,  N,  1, 1 );
    CHAMELEON_Desc_Create(
        &descQ,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, Qm, 0, 0, Qm, Qm, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, C,  LDC, descC );
    CHAMELEON_zLap2Desc( ChamUpperLower, CC, LDC, descCC );
    CHAMELEON_zLap2Desc( ChamUpperLower, Q,  LDQ, descQ );

    info = check_zqc( args, side, trans, descC, descQ, descCC );

    CHAMELEON_Desc_Destroy( &descC );
    CHAMELEON_Desc_Destroy( &descCC );
    CHAMELEON_Desc_Destroy( &descQ );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] Bnorm
 *          TODO
 *
 * @param[in] descR
 *          The descriptor of the matrix R.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgeqrs( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR )
{
    int info_local, info_global, nb;
    int M = descA->m;
    int N = descA->n;
    int NRHS = descR->n;
    int maxMNK = chameleon_max( M, chameleon_max( N, NRHS ) );
    double Rnorm, result;
    double Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    double eps = LAPACKE_dlamch_work('e');

    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    if ( trans == ChamNoTrans ) {
        CHAM_desc_t *descRR;
        /*
         * Corresponds to lapack/testings/lin/[sdcz]qrt17.f
         *
         * ZQRT17 computes the ratio
         *
         *    || R'*op(A) ||/(||A||*alpha*max(M,N,NRHS)*eps)
         *
         * where R = op(A)*X - B, op(A) is A or A', and alpha = ||B||
         *
         */
        CHAMELEON_Desc_Create( &descRR, NULL, ChamComplexDouble, nb, nb, nb*nb,
                               NRHS, N, 0, 0, NRHS, N, descA->p, descA->q );

        CHAMELEON_zgemm_Tile( ChamConjTrans, trans, 1., descR, descA, 0., descRR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descRR );
        result = Rnorm / (Anorm * Bnorm * eps * maxMNK);
        CHAMELEON_Desc_Destroy( &descRR );
    }
    else {
        /*
         * To implement this test, we need to look at LAPACK working note 41, page 29
         * and more especially to lapack/testings/lin/[sdcz]qrt14.f
         */
        fprintf(stderr, "GEQRS testing not implemented with M >= N when transA = ChamConjTrans\n");
        return 0;
    }

    if ( isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] Bnorm
 *          TODO
 *
 * @param[in] descR
 *          The descriptor of the matrix R.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgelqs( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR )
{
    int info_local, info_global, nb;
    int M = descA->m;
    int N = descA->n;
    int NRHS = descR->n;
    int maxMNK = chameleon_max( M, chameleon_max( N, NRHS ) );
    double Rnorm, result;
    double Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    double eps = LAPACKE_dlamch_work('e');

    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    if ( trans == ChamNoTrans ) {
        /*
         * To implement this test, we need to look at LAPACK working note 41, page 29
         * and more especially to lapack/testings/lin/[sdcz]lqt14.f
         */
        fprintf(stderr, "GELQS testing not implemented with N > M when transA = ChamNoTrans\n");
        return 0;
    }
    else {
        CHAM_desc_t *descRR;
        /*
         * Corresponds to lapack/testings/lin/[sdcz]qrt17.f
         *
         * ZQRT17 computes the ratio
         *
         *    || R'*op(A) ||/(||A||*alpha*max(M,N,NRHS)*eps)
         *
         * where R = op(A)*X - B, op(A) is A or A', and alpha = ||B||
         *
         */
        CHAMELEON_Desc_Create( &descRR, NULL, ChamComplexDouble, nb, nb, nb*nb, NRHS, M, 0, 0, NRHS, M, descA->p, descA->q );

        CHAMELEON_zgemm_Tile( ChamConjTrans, trans, 1., descR, descA, 0., descRR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descRR );
        result = Rnorm / (Anorm * Bnorm * eps * maxMNK);
        CHAMELEON_Desc_Destroy( &descRR );
    }

    if ( isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descX
 *          The descriptor of the matrix X.
 *
 * @param[in] descB
 *          The descriptor of the matrix B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgels( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB )
{
    int info_solution;
    int M = descA->m;
    int N = descA->n;
    double Bnorm = CHAMELEON_zlange_Tile( ChamInfNorm, descB );

    info_solution = check_zsolve( args, ChamGeneral, trans, ChamUpperLower,
                                  descA, descX, descB );

    if ( M >= N ) {
        info_solution = check_zgeqrs( args, trans, descA, Bnorm, descB );
    }
    else {
        info_solution = check_zgelqs( args, trans, descA, Bnorm, descB );
    }

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
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] M
 *        The number of rows of the matrix A.
 *
 * @param[in] N
 *        The number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] X
 *          The matrix X.
 * 
 * @param[in] LDX
 *          The leading dimension of the matrix X.
 *
 * @param[in] B
 *          The matrix B.
 * 
 * @param[in] LDB
 *          The leading dimension fo the matrix B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgels_std( run_arg_list_t *args, cham_trans_t trans, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, int LDX, CHAMELEON_Complex64_t *B, int LDB )
{
    int info_solution;
    CHAM_desc_t *descA, *descX, *descB;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );
    int Xm = (trans == ChamNoTrans) ? N : M;
    int Bm = (trans == ChamNoTrans) ? M : N;

    CHAMELEON_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N,    0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descX, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDX, NRHS, 0, 0, Xm, NRHS, 1, 1 );
    CHAMELEON_Desc_Create(
        &descB, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, Bm, NRHS, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A, LDA, descA );
    CHAMELEON_zLap2Desc( ChamUpperLower, X, LDX, descX );
    CHAMELEON_zLap2Desc( ChamUpperLower, B, LDB, descB );

    info_solution = check_zgels( args, trans, descA, descX, descB );

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );
    CHAMELEON_Desc_Destroy( &descB );

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
    double eps = LAPACKE_dlamch_work('e');

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

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Extends the check_zgeqrf and check_zortho to check the specific QR
 * factorization and Q generation used in Polar decompositions.
 *
 * [ A1 ] = [ Q1 ] [ AF1 ]
 * [ A2 ]   [ Q2 ]
 *
 *******************************************************************************
 *
 * @param[in] descA1
 *          The descriptor of the top of the matrix A.
 *
 * @param[in] descA2
 *          The descriptor of the bottom of the matrix A.
 *
 * @param[in] descQ1
 *          The descriptor of the top of the matrix Q generated from the
 *          QR factorization of the matrix A.
 *
 * @param[in] descQ2
 *          The descriptor of the bottom of the matrix Q generated from the
 *          QR factorization of the matrix A.
 *
 * @param[in] descAF1
 *          The descriptor of the top of the QR factorization of the matrix A that holds R.
 *
 * @retval 0  on success
 * @retval >0 on failure
 *
 *******************************************************************************
 */
int check_zgepdf_qr( run_arg_list_t *args, CHAM_desc_t *descA1, CHAM_desc_t *descA2, CHAM_desc_t *descQ1,
                     CHAM_desc_t *descQ2, CHAM_desc_t *descAF1 )
{
    int info_local, info_global;
    int M = (descQ1->m + descQ2->m);
    int N = descQ1->n;
    int K = descAF1->n;
    double result, Anorm, A1norm, A2norm, Rnorm;
    double eps = LAPACKE_dlamch_work('e');
    CHAM_desc_t *descR, *subR, *subAF;

    /*
     * We exploit the fact that the lower matrix is supposed to be smaller than
     * the upper one, and the fact that K is always equal to N in this specific
     * problem.
     */
    assert( descAF1->m >= N );
    assert( K == N );

    /*
     * First check: || R - Q' A ||
     */
    descR = CHAMELEON_Desc_Copy( descAF1, NULL );

    /* Copy R(1:n,1:k) */
    subR  = chameleon_desc_submatrix( descR,   0, 0, N, K );
    subAF = chameleon_desc_submatrix( descAF1, 0, 0, N, K );

    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., subR );
    CHAMELEON_zlacpy_Tile( ChamUpper, subAF, subR );
    free( subAF );

    /*
     * Compute R(1:n,1:k) - Q(:,1:n)' * A(:,1:k)
     *       = R(1:n,1:k) - Q1(:,1:n)' * A1(:,1:k) - Q2(,1:n)' * A2(,1:k)
     */
    CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ1, descA1, 1., subR );
    CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ2, descA2, 1., subR );

    Rnorm  = CHAMELEON_zlange_Tile( ChamOneNorm, subR );
    A1norm = CHAMELEON_zlange_Tile( ChamOneNorm, descA1 );
    A2norm = CHAMELEON_zlange_Tile( ChamOneNorm, descA2 );
    Anorm  = A1norm + A2norm; /* This is an upper bound exact when A2 is diagonal due to OneNorm */
    result = Rnorm / ( (double)M * Anorm * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /*
     * Second check: || I - Q' Q ||
     */
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., subR );

    /* Performs I - Q'Q = I - Q1'Q1 - Q2'Q2 */
    CHAMELEON_zherk_Tile( ChamUpper, ChamConjTrans, -1., descQ1, 1., subR );
    CHAMELEON_zherk_Tile( ChamUpper, ChamConjTrans, -1., descQ2, 1., subR );

    /* Verifies the residual's norm */
    Rnorm = CHAMELEON_zlansy_Tile( ChamOneNorm, ChamUpper, subR );
    result = Rnorm / ( (double)K * eps );

    run_arg_add_double( args, "||I-QQ'||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local++;
    }

    free( subR );
    CHAMELEON_Desc_Destroy( &descR );

    /* Reduces the result on all processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a Chameleon Polar Decomposition is correct.
 *
 * @Warning Check only the general case for now.
 *
 *******************************************************************************
 *
 * @param[in,out] descA
 *          The descriptor of the matrix A, on exit the matrix is modified.
 *
 * @param[in] descU
 *          The descriptor of the orthogonal polar factor of the decomposition.
 *
 * @param[in] descH
 *          The descriptor of the symmetric/hermitian polar factor of the decomposition.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxpd( run_arg_list_t *args,
                 CHAM_desc_t *descA, CHAM_desc_t *descU, CHAM_desc_t *descH )
{
    int info_local, info_global;
    double Anorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');

    /* Compute ||A|| */
    Anorm = CHAMELEON_zlange_Tile( ChamFrobeniusNorm, descA );

    /* R = A - U * H */
    CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, 1., descU, descH, -1., descA );

    /* Compute ||R|| */
    Rnorm = CHAMELEON_zlange_Tile( ChamFrobeniusNorm, descA );

    result = Rnorm / (Anorm * eps);
    run_arg_add_double( args, "||A||",         Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else{
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a Polar Decomposition is correct.
 *
 * @Warning Check only the general case for now.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrices A and U.
 *
 * @param[in] N
 *          The number of columns of the matrices A and U and the order of the matrix H.
 *
 * @param[in,out] A
 *          The matrix A, on exit the matrix is modified.
 *
 * @param[in] U
 *          The matrix with the orthogonal polar factor of the decomposition.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A and U.
 *
 * @param[in] H
 *          The matrix of the symmetric/hermitian polar factor of the decomposition.
 *
 * @param[in] LDH
 *          The leading dimension of the matrix H.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxpd_std( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *A,
                     CHAMELEON_Complex64_t *U, int LDA, CHAMELEON_Complex64_t *H, int LDH )
{
    int          info;
    CHAM_desc_t *descA, *descU, *descH;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descU, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descH, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDH, N, 0, 0, N, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A, LDA, descA );
    CHAMELEON_zLap2Desc( ChamUpperLower, U, LDA, descU );
    CHAMELEON_zLap2Desc( ChamUpperLower, H, LDH, descH );

    info = check_zxxpd( args, descA, descU, descH );

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descU );
    CHAMELEON_Desc_Destroy( &descH );

    return info;
}

#endif /* defined(CHAMELEON_SIMULATION) */
