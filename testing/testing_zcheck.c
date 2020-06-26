/**
 *
 * @file testing_zcheck.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
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
#include "../control/common.h"
#include "testings.h"
#include "testing_zcheck.h"
#include "flops.h"

#ifndef max
#define max( _a_, _b_ ) ( (_a_) > (_b_) ? (_a_) : (_b_) )
#endif

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Compares two matrices by their norms.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Wether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] descA
 *          The first matrix descriptor.
 *
 * @param[in] descA2
 *          The second matrix descriptor.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zmatrices( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descB )
{
    int info_solution = 0;
    int M = descA->m;
    int N = descB->n;
    int LDA = descA->m;
    int rank = CHAMELEON_Comm_rank();
    double Anorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    CHAMELEON_Complex64_t *A = NULL;
    CHAMELEON_Complex64_t *B = NULL;

    if ( rank == 0 ) {
        A = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
        B = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
    }

    /* Converts the matrices to LAPACK layout in order to compare them on the main process */
    CHAMELEON_zDesc2Lap( uplo, descA, A, LDA );
    CHAMELEON_zDesc2Lap( uplo, descB, B, LDA );

    if ( rank == 0 ) {
        double *work = (double *)malloc(LDA*N*sizeof(double));

        /* Computes the norms */
        if ( uplo == ChamUpperLower ) {
            Anorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, A, LDA, work );
        }
        else {
            Anorm = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), 'N',
                                         M, N, A, LDA, work );
        }

        /* Computes the difference with the core function */
        CORE_zgeadd( ChamNoTrans, M, N, 1, A, LDA, -1, B, LDA );

        /* Computes the residual's norm */
        if ( uplo == ChamUpperLower ) {
            Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, B, LDA, work );
        }
        else {
            Rnorm = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), 'N',
                                         M, N, B, LDA, work );
        }
        if ( Anorm != 0. ) {
            result = Rnorm / (Anorm * eps);
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

        free(work);
        free(A);
        free(B);
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
 * @ingroup new_testing
 *
 * @brief Compares the Chameleon computed norm with a Lapack computed one.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Wether it is a general, triangular or symmetric matrix.
 *
 * @param[in] norm_type
 *          Wether it should compare a Max, One, Inf or Frobenius norm.
 *
 * @param[in] uplo
 *          Wether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] diag
 *          Wether it is a unitary diagonal matrix or not.
 *
 * @param[in] norm_cham
 *          The Chameleon computed norm.
 *
 * @param[in] descA
 *          The matrix descriptor.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_znorm( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_normtype_t norm_type, cham_uplo_t uplo,
                 cham_diag_t diag, double norm_cham, CHAM_desc_t *descA )
{
    int info_solution = 0;
    int M = descA->m;
    int N = descA->n;
    int LDA = descA->m;
    int rank = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A = NULL;

    if ( rank == 0 ) {
        A = (CHAMELEON_Complex64_t *)malloc(N*LDA*sizeof(CHAMELEON_Complex64_t));
    }

    /* Converts the matrix to LAPACK layout in order to use the LAPACK norm function */
    CHAMELEON_zDesc2Lap( uplo, descA, A, LDA );

    if ( rank == 0 ) {
        double *work = (double*) malloc(chameleon_max(M, N)*sizeof(double));
        double norm_lapack;
        double result;
        double eps = LAPACKE_dlamch_work('e');

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
            fprintf(stderr, "check_znorm: mtxtype(%d) unsupported\n", matrix_type );
            free( work );
            return 1;
        }

        /* Compares the norms */
        result = fabs( norm_cham - norm_lapack ) / ( norm_lapack * eps );

        switch(norm_type) {
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
            /* Sum order on every element can differ */
            result = result / ((double)M * (double)N);
            break;
        }

        info_solution = ( result < 1 ) ? 0 : 1;

        free(work);
        free(A);
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
 * @ingroup new_testing
 *
 * @brief Compares a Chameleon computed sum with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Wether it is a upper triangular matrix, a lower triangular matrix or
 *          a general matrix.
 *
 * @param[in] trans
 *          Wether the first matrix is transposed, conjugate transposed or not
 *          transposed.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descBref
 *          The descriptor of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descBcham
 *          The matrix descriptor of the Chameleon computed result A+B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsum ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                 CHAMELEON_Complex64_t beta, CHAM_desc_t *descBref, CHAM_desc_t *descBcham )
{
    int info_solution = 0;
    int M = descBref->m;
    int N = descBref->n;
    int Am = (trans == ChamNoTrans)? M : N;
    int An = (trans == ChamNoTrans)? N : M;
    int LDA = Am;
    int LDB = M;
    int rank = CHAMELEON_Comm_rank();
    double Anorm, Binitnorm, Rnorm, result;
    CHAMELEON_Complex64_t *A     = NULL;
    CHAMELEON_Complex64_t *Bref  = NULL;
    CHAMELEON_Complex64_t *Bcham = NULL;
    CHAMELEON_Complex64_t  mzone = -1.0;
    cham_uplo_t uploA = uplo;

    if ( rank == 0 ) {
        A     = malloc(An*LDA*sizeof(CHAMELEON_Complex64_t));
        Bref  = malloc(N*LDB*sizeof(CHAMELEON_Complex64_t));
        Bcham = malloc(N*LDB*sizeof(CHAMELEON_Complex64_t));
    }

    /* Computes the max norms of A, B and A+B */
    if ( uplo == ChamUpperLower ) {
        Anorm     = CHAMELEON_zlange_Tile( ChamMaxNorm, descA     );
        Binitnorm = CHAMELEON_zlange_Tile( ChamMaxNorm, descBref  );
    }
    else {
        if ( trans == ChamNoTrans ) {
            Anorm = CHAMELEON_zlantr_Tile( ChamMaxNorm, uplo, ChamNonUnit, descA );
        }
        else {
            uploA = (uplo == ChamUpper) ? ChamLower : ChamUpper;
            Anorm = CHAMELEON_zlantr_Tile( ChamMaxNorm, uploA, ChamNonUnit, descA );
        }
        Binitnorm = CHAMELEON_zlantr_Tile( ChamMaxNorm, uplo, ChamNonUnit, descBref );
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uploA, descA,     A,     LDA );
    CHAMELEON_zDesc2Lap( uplo,  descBref,  Bref,  LDB );
    CHAMELEON_zDesc2Lap( uplo,  descBcham, Bcham, LDB );

    if ( rank == 0 ) {
        double  eps  = LAPACKE_dlamch_work('e');
        double *work = malloc(chameleon_max(M, N)* sizeof(double));

        /* Makes the sum with the core function */
        if ( uplo == ChamUpperLower ) {
            CORE_zgeadd( trans, M, N,
                         alpha, A,    LDA,
                         beta,  Bref, LDB );
        }
        else {
            CORE_ztradd( uplo, trans, M, N,
                         alpha, A,    LDA,
                         beta,  Bref, LDB );
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

        /* Verifies if the result is inside a threshold */
        if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
            info_solution = 1;
        }
        else {
            info_solution = 0;
        }

        free(work);
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
 * @ingroup new_testing
 *
 * @brief Compares a Chameleon computed scale with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Wether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA1
 *          The original matrix descriptor.
 *
 * @param[in] descA2
 *          The scaled matrix descriptor.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zscale( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA1, CHAM_desc_t *descA2 )
{
    int info_solution;
    int M = descA1->m;
    int N = descA1->n;
    int rank = CHAMELEON_Comm_rank();
    CHAM_desc_t *descBlas;
    CHAMELEON_Complex64_t *Ainit = NULL;

    if ( rank == 0 ) {
        Ainit = (CHAMELEON_Complex64_t *)malloc(M*N*sizeof(CHAMELEON_Complex64_t));
    }

    /* Converts the matrix to LAPACK layout in order to scale with BLAS */
    CHAMELEON_zDesc2Lap( uplo, descA1, Ainit, M );

    if ( rank == 0 ) {
        /* Scales using core function */
        CORE_zlascal( uplo, M, N, alpha, Ainit, M );
    }

    /* Converts back into Chameleon to compare with check_zmatrices */
    descBlas = CHAMELEON_Desc_CopyOnZero( descA1, NULL );
    CHAMELEON_zLap2Desc( uplo, Ainit, M, descBlas );

    /* Compares the two matrices */
    info_solution = check_zmatrices( args, uplo, descA2, descBlas );

    if ( rank == 0 ) {
        free( Ainit );
    }

    CHAMELEON_Desc_Destroy( &descBlas );

    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Compares a Chameleon computed product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Wether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] transB
 *          Wether the second product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descBref
 *          The descriptor of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the matrix C.
 *
 * @param[in] descC
 *          The matrix descriptor of the Chameleon computed result alpha*A*B+beta*C.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgemm( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                 CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int An, LDA, Bn, LDB;
    int info_solution = 0;
    int M = descC->m;
    int N = descC->n;
    int K = (transA != ChamNoTrans)? descA->m : descA->n;
    int LDC = descC->m;
    int rank = CHAMELEON_Comm_rank();
    double Anorm, Bnorm, Crefnorm, Rnorm, result;
    CHAMELEON_Complex64_t *A = NULL;
    CHAMELEON_Complex64_t *B = NULL;
    CHAMELEON_Complex64_t *C = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    CHAMELEON_Complex64_t mzone = -1.0;

    /* Calculates the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Anorm = CHAMELEON_zlange_Tile(ChamInfNorm, descA);
        LDA = M;
        An  = K;
    } else {
        Anorm = CHAMELEON_zlange_Tile(ChamOneNorm, descA);
        LDA = K;
        An  = M;
    }
    if ( transB == ChamNoTrans ) {
        Bnorm = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        LDB = K;
        Bn  = N;
    } else {
        Bnorm = CHAMELEON_zlange_Tile(ChamOneNorm, descB);
        LDB = N;
        Bn  = K;
    }

    /* Computes the norms for comparing */
    Crefnorm = CHAMELEON_zlange_Tile(ChamMaxNorm, descCref);

    /* Creates the LAPACK version of the matrices */
    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc(An*LDA*sizeof(CHAMELEON_Complex64_t));
        B    = (CHAMELEON_Complex64_t *)malloc(Bn*LDB*sizeof(CHAMELEON_Complex64_t));
        Cref = (CHAMELEON_Complex64_t *)malloc(N *LDC*sizeof(CHAMELEON_Complex64_t));
        C    = (CHAMELEON_Complex64_t *)malloc(N *LDC*sizeof(CHAMELEON_Complex64_t));
    }

    CHAMELEON_zDesc2Lap( ChamUpperLower, descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descC,    C,    LDC );

    if ( rank == 0 ) {
        double eps = LAPACKE_dlamch_work('e');

        /* Makes the multiplication with the core function */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K,
                     CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
        cblas_zaxpy(LDC * N, CBLAS_SADDR(mzone), C, 1, Cref, 1);

        /* Calculates the norm with the core function's result */
        Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );

        if ( ( alpha != 0. ) || (beta != 0. ) ) {
            result = Rnorm / ((cabs(alpha) * max(Anorm, Bnorm) + cabs(beta) * Crefnorm) * K * eps);
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

        free(A);
        free(B);
        free(C);
        free(Cref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Compares a Chameleon computed hermitian product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Wether the hermitian matrix A appears on the left or right in the operation.
 *
 * @param[in] uplo
 *          Wether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the hermitian matrix A.
 *
 * @param[in] descB
 *          The descriptor of the hermitian matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the hermitian matrix C.
 *
 * @param[in] descC
 *          The matrix descriptor of the Chameleon computed result alpha*A*B+beta*C or alpha*B*A+beta*C.
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
    int An, LDA;
    int M = descC->m;
    int N = descC->n;
    int LDB = M;
    int LDC = M;
    int rank = CHAMELEON_Comm_rank();
    double Anorm, Bnorm, Crefnorm, Cchamnorm, Clapacknorm, Rnorm, result;
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    CHAMELEON_Complex64_t *C    = NULL;
    CHAMELEON_Complex64_t mzone = -1.0;

    if ( side == ChamLeft ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
        if ( matrix_type == ChamHermitian ) {
            Anorm = CHAMELEON_zlanhe_Tile(ChamInfNorm, uplo, descA);
        }
        else
#endif
        {
            Anorm = CHAMELEON_zlansy_Tile(ChamInfNorm, uplo, descA);
        }
        Bnorm = CHAMELEON_zlange_Tile(ChamOneNorm, descB);
        LDA = M;
        An  = M;
    } else {
#if defined(PRECISION_z) || defined(PRECISION_c)
        if ( matrix_type == ChamHermitian ) {
            Anorm = CHAMELEON_zlanhe_Tile(ChamOneNorm, uplo, descA);
        }
        else
#endif
        {
            Anorm = CHAMELEON_zlansy_Tile(ChamOneNorm, uplo, descA);
        }
        Bnorm = CHAMELEON_zlange_Tile(ChamInfNorm, descB);
        LDA = N;
        An  = N;
    }

    /* Computes the norms for comparing */
    Crefnorm  = CHAMELEON_zlange_Tile( ChamMaxNorm, descCref );
    Cchamnorm = CHAMELEON_zlange_Tile( ChamMaxNorm, descC    );

    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc(LDA * An * sizeof(CHAMELEON_Complex64_t));
        B    = (CHAMELEON_Complex64_t *)malloc(LDB * N  * sizeof(CHAMELEON_Complex64_t));
        Cref = (CHAMELEON_Complex64_t *)malloc(LDC * N  * sizeof(CHAMELEON_Complex64_t));
        C    = (CHAMELEON_Complex64_t *)malloc(LDC * N  * sizeof(CHAMELEON_Complex64_t));
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uplo,           descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descC,    C,    LDC );

    if ( rank == 0 ) {
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
            cblas_zsymm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                         M, N, CBLAS_SADDR(alpha),
                         A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
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
            info_solution= 0 ;
        }

        free(A);
        free(B);
        free(C);
        free(Cref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    (void)Clapacknorm;
    (void)Cchamnorm;
    (void)matrix_type;
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Compares a Chameleon computed matrix rank k operation with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Wether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] trans
 *          Wether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B - only used for her2k and sy2k.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the matrix C.
 *
 * @param[in] descC
 *          The matrix descriptor of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsyrk( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_trans_t trans,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                 CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int LDA, info_solution = 0;
    int An, K, N = descC->n;
    int LDC = N;
    int rank = CHAMELEON_Comm_rank();
    double Anorm, Bnorm, Crefnorm, Cchamnorm, Clapacknorm, Rnorm, result;
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    CHAMELEON_Complex64_t *C    = NULL;

    Bnorm = 0.;
    if ( trans == ChamNoTrans ) {
        Anorm = CHAMELEON_zlange_Tile( ChamInfNorm, descA );
        if ( descB != NULL ) {
            Bnorm = CHAMELEON_zlange_Tile( ChamInfNorm, descB );
        }
        K = descA->n;
        LDA = N;
        An  = K;
    }
    else {
        Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
        if ( descB != NULL ) {
            Bnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
        }
        K = descA->m;
        LDA = K;
        An  = N;
    }

    /* Computes the norms for comparing */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        Crefnorm  = CHAMELEON_zlanhe_Tile( ChamInfNorm, uplo, descCref );
        Cchamnorm = CHAMELEON_zlanhe_Tile( ChamInfNorm, uplo, descC    );
    }
    else
#endif
    {
        Crefnorm  = CHAMELEON_zlansy_Tile( ChamInfNorm, uplo, descCref );
        Cchamnorm = CHAMELEON_zlansy_Tile( ChamInfNorm, uplo, descC    );
    }

    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc(LDA * An * sizeof(CHAMELEON_Complex64_t));
        if ( descB != NULL ) {
            B = (CHAMELEON_Complex64_t *)malloc(LDA * An * sizeof(CHAMELEON_Complex64_t));
        }
        Cref = (CHAMELEON_Complex64_t *)malloc(LDC * N * sizeof(CHAMELEON_Complex64_t));
        C    = (CHAMELEON_Complex64_t *)malloc(LDC * N * sizeof(CHAMELEON_Complex64_t));
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( ChamUpperLower, descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( uplo,           descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( uplo,           descC,    C,    LDC );
    if ( descB != NULL ) {
        CHAMELEON_zDesc2Lap( ChamUpperLower, descB, B, LDA );
    }

    if ( rank == 0 ) {
        double eps = LAPACKE_dlamch_work('e');
        double ABnorm;
        double *work = malloc(sizeof(double)*N);

        /* Makes the multiplication with the core function */
#if defined(PRECISION_z) || defined(PRECISION_c)
        if ( matrix_type == ChamHermitian ) {
            if ( descB == NULL ) {
                cblas_zherk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                             N, K, creal(alpha), A, LDA, creal(beta), Cref, LDC );
                ABnorm = Anorm * Anorm;
            }
            else {
                cblas_zher2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                              N, K, CBLAS_SADDR(alpha), A, LDA, B, LDA, creal(beta), Cref, LDC );
                ABnorm = 2. * Anorm * Bnorm;
            }

            Clapacknorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
        }
        else
#endif
        {
            if ( descB == NULL ) {
                cblas_zsyrk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                             N, K, CBLAS_SADDR(alpha), A, LDA, CBLAS_SADDR(beta), Cref, LDC );
                ABnorm = Anorm * Anorm;
            }
            else {
                cblas_zsyr2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                              N, K, CBLAS_SADDR(alpha), A, LDA, B, LDA, CBLAS_SADDR(beta), Cref, LDC );
                ABnorm = 2. * Anorm * Bnorm;
            }

            Clapacknorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
        }

        CORE_ztradd( uplo, ChamNoTrans, N, N,
                     -1.,  C,    LDC,
                      1.,  Cref, LDC );

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

        /* Verifies if the result is inside a threshold */
        if ( isinf(Clapacknorm) || isinf(Cchamnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
            info_solution = 1;
        }
        else {
            info_solution = 0;
        }

        free(work);
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

    (void)matrix_type;
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Compares a Chameleon computed matrix triangular product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] check_func
 *          Wether it is a triangular product or a triangular linear solution.
 *
 * @param[in] side
 *          Whether A appears on the left or on the right of the product.
 *
 * @param[in] uplo
 *          Wether A is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] trans
 *          Wether A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] diag
 *          Wether A is a unitary diagonal matrix or not.
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
 * @param[in] descBref
 *          The descriptor of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrmm( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB, CHAM_desc_t *descBref )
{
    int info_solution = 0;
    int M = descB->m;
    int N = descB->n;
    int An, LDA, LDB = M;
    int rank = CHAMELEON_Comm_rank();
    double Anorm, Bnorm, Brefnorm, Rnorm, result;
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *Bref = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t mzone = -1.0;

    /* Computes the norms for comparing */
    if ( ((side == ChamLeft)  && (trans == ChamNoTrans)) ||
         ((side == ChamRight) && (trans != ChamNoTrans)) ) {
        Anorm = CHAMELEON_zlantr_Tile( ChamInfNorm, uplo, diag, descA );
    }
    else {
        Anorm = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, diag, descA );
    }
    if ( side == ChamLeft ) {
        An = M;
        LDA = M;
        Brefnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descBref );
        Bnorm    = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
    }
    else {
        An = N;
        LDA = N;
        Brefnorm = CHAMELEON_zlange_Tile( ChamInfNorm, descBref );
        Bnorm    = CHAMELEON_zlange_Tile( ChamInfNorm, descB );
    }

    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc(An*LDA*sizeof(CHAMELEON_Complex64_t));
        Bref = (CHAMELEON_Complex64_t *)malloc(N *LDB*sizeof(CHAMELEON_Complex64_t));
        B    = (CHAMELEON_Complex64_t *)malloc(N *LDB*sizeof(CHAMELEON_Complex64_t));
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uplo,           descA,    A,    LDA);
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB);
    CHAMELEON_zDesc2Lap( ChamUpperLower, descBref, Bref, LDB);

    if ( rank == 0 ) {
        double eps = LAPACKE_dlamch_work('e');

        /* Makes the multiplication with the core function */
        if (check_func == CHECK_TRMM) {
            cblas_ztrmm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                        (CBLAS_DIAG)diag, M, N, CBLAS_SADDR(alpha), A, LDA, Bref, LDB);
        }
        else {
            cblas_ztrsm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                        (CBLAS_DIAG)diag, M, N, CBLAS_SADDR(alpha), A, LDA, Bref, LDB);
        }

        /* Computes the norm with the core function's result */
        cblas_zaxpy( LDB * N, CBLAS_SADDR(mzone), B, 1, Bref, 1 );
        Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Bref, LDB, NULL );

        result = Rnorm / ((Anorm + Brefnorm) * An * eps);

        /* Verifies if the result is inside a threshold */
        if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
            info_solution = 1;
        }
        else {
            info_solution = 0;
        }

        free(A);
        free(B);
        free(Bref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    (void)Bnorm;
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Compares a Chameleon computed product U*U' or L'*L result with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Wether the upper or lower triangle of A is stored.
 *
 * @param[in] descA1
 *          The descriptor of the A matrix.
 *
 * @param[in] descA2
 *          The descriptor of the Chameleon computed result matrix.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zlauum( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descAAt )
{
    int info_local, info_global;
    int N = descA->n;
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
 * @ingroup new_testing
 *
 * @brief Checks if a Chameleon computed factorization is correct.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Wether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] descA
 *          The descriptor of the original matrix.
 *
 * @param[in] descB
 *          The descriptor of the Chameleon factorized matrix.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxtrf( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo,
                  CHAM_desc_t *descA, CHAM_desc_t *descLU )
{
    int info_local, info_global;
    int M = descA->m;
    int N = descA->n;
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
        transL = (mtxtype == ChamHermitian) ? ChamConjTrans : ChamTrans;
#else
        transL = ChamTrans;
#endif
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descL );
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descU );
        break;
    case ChamLower:
#if defined(PRECISION_z) || defined(PRECISION_c)
        transU = (mtxtype == ChamHermitian) ? ChamConjTrans : ChamTrans;
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

    switch ( mtxtype ) {
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
        fprintf(stderr, "check_zxxtrf: mtxtype(%d) unsupported\n", mtxtype );
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

    (void)args;
    return info_global;
}


/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Checks if the  linear solution of op(A) * x = b is correct.
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 *
 * @param[in] trans
 *          Wether the A matrix is non transposed, tranposed or conjugate transposed.
 *
 * @param[in] uplo
 *
 * @param[in] descA
 *          The descriptor of the A matrix.
 *
 * @param[in] descX
 *          The descriptor of the X matrix.
 *
 * @param[inout] descB
 *          The descriptor of the B = A*X matrix. On exit, it contains the remainder from A*x-B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsolve( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_trans_t trans, cham_uplo_t uplo,
                  CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB )
{
    int info_local, info_global;
    int M = descA->m;
    int N = descA->n;
    double Anorm, Bnorm, Xnorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    cham_normtype_t norm = (trans == ChamNoTrans) ? ChamOneNorm : ChamInfNorm;

    /* Computes the norms */
    Bnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
    Xnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descX );

    switch ( mtxtype ) {
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
        fprintf(stderr, "check_zsolve: mtxtype(%d) unsupported\n", mtxtype );
        return 1;
    }

    Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
    result = Rnorm / ( Anorm * Xnorm * chameleon_max( M, N ) * eps );

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
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

    (void)Bnorm;
    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Checks if the A1 matrix is the inverse of A2.
 *
 *******************************************************************************
 *
 * @param[in] is_herm
 *          Wether the matrices are hermitian.
 *
 * @param[in] uplo
 *          Wether they are upper triangular matrices or lower triangular matrices.
 *
 * @param[in] diag
 *          Wether they are unitary diagonal matrices or not.
 *
 * @param[in] descA1
 *          The descriptor of the A1 matrix.
 *
 * @param[in] descX
 *          The descriptor of the A2 matrix.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrtri( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_diag_t diag,
                  CHAM_desc_t *descA0, CHAM_desc_t *descAi )
{
    int info_local, info_global;
    int N = descA0->m;
    cham_uplo_t uplo_inv;
    CHAM_desc_t *descI, *descB = NULL;
    double Rnorm, Anorm, Ainvnorm, result;
    double eps = LAPACKE_dlamch_work('e');

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

int check_zortho( run_arg_list_t *args, CHAM_desc_t *descQ )
{
    int info_local, info_global;
    int M = descQ->m;
    int N = descQ->n;
    int minMN = chameleon_min(M, N);
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

    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup new_testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *          The descriptor of the initial matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the factorized matrix A.
 *
 * @param[in] descQ
 *          The descriptor of the Q matrix generated with a call to ungqr and
 *          the factorized matrix A (descAF).
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

    descL = CHAMELEON_Desc_Copy( descA, NULL );

    if ( (K < chameleon_min( M, N )) && (N >= M) ) {
        /*
         * Cas lapack zlqt02.f
         */
        CHAM_desc_t *subL, *subAF;

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
    else {
        /*
         * Cas lapack zlqt01.f
         */

        /* Copy L */
        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descL );
        CHAMELEON_zlacpy_Tile( ChamLower, descAF, descL );

        /* Compute L - A * Q' */
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamConjTrans, -1., descA, descQ, 1., descL );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descL );
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
 * @ingroup new_testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *          The descriptor of the initial matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the factorized matrix A.
 *
 * @param[in] descQ
 *          The descriptor of the Q matrix generated with a call to ungqr and
 *          the factorized matrix A (descAF).
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

    descR = CHAMELEON_Desc_Copy( descA, NULL );

    if ( (K < chameleon_min( M, N )) && (M >= N) ) {
        /*
         * Cas lapack zqrt02.f
         */
        CHAM_desc_t *subR, *subAF;

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
    else {
        /*
         * Cas lapack zqrt01.f
         */

        /* Copy R */
        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descR );
        CHAMELEON_zlacpy_Tile( ChamUpper, descAF, descR );

        /* Compute R - Q'*A */
        CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ, descA, 1., descR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descR );
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

    if (  isnan(CCnorm) || isinf(CCnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
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

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
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
        CHAMELEON_Desc_Create( &descRR, NULL, ChamComplexDouble, nb, nb, nb*nb,
                               NRHS, M, 0, 0, NRHS, M, descA->p, descA->q );

        CHAMELEON_zgemm_Tile( ChamConjTrans, trans, 1., descR, descA, 0., descRR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descRR );
        result = Rnorm / (Anorm * Bnorm * eps * maxMNK);
        CHAMELEON_Desc_Destroy( &descRR );
    }

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
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

#endif /* defined(CHAMELEON_SIMULATION) */
