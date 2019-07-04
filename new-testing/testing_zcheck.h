/**
 *
 * @file testing_zcheck.h
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings header
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-07-16
 * @precisions normal z -> c d s
 *
 */

#ifndef _testing_zcheck_h_
#define _testing_zcheck_h_

#include "testings.h"
#include <math.h>

#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif

#define CHECK_TRMM 3
#define CHECK_TRSM 4

void print_zmatrix      ( int M, int N, CHAMELEON_Complex64_t *A, int LDA );
void print_zdesc_matrix ( CHAM_desc_t *descA );
void zsabotage          ( CHAM_desc_t *descA );
void potri_product      ( cham_uplo_t uplo, CHAM_desc_t *descA1, CHAM_desc_t *descA2 );

int check_zmatrices     ( cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descB );
int check_znorm         ( cham_mtxtype_t mtxtype, cham_normtype_t norm_type, cham_uplo_t uplo,
                          cham_diag_t diag, double norm_cham, CHAM_desc_t *descA );
int check_zsum          ( cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                          CHAMELEON_Complex64_t beta, CHAM_desc_t *descBref, CHAM_desc_t *descBcham );
int check_zscale        ( cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA1, CHAM_desc_t *descA2 );
int check_zgemm         ( cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                          CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC );
int check_zsymm         ( cham_mtxtype_t mtxtype, cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                          CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC );
int check_zsyrk         ( cham_mtxtype_t mtxtype, cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                          CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC );
int check_ztrmm         ( int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB, CHAM_desc_t *descBref );
int check_zlauum        ( cham_uplo_t uplo, CHAM_desc_t *descA1, CHAM_desc_t *descA2 );
int check_zxxtrf        ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, CHAM_desc_t *descA1, CHAM_desc_t *descA2 );
int check_zsolve        ( cham_mtxtype_t mtxtype, cham_trans_t trans, cham_uplo_t uplo,
                          CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB );
int check_ztrtri        ( cham_mtxtype_t mtxtype, cham_uplo_t uplo, cham_diag_t diag,
                          CHAM_desc_t *descA, CHAM_desc_t *descAi );

/* Using QR factorization */
int check_zortho        ( CHAM_desc_t *descQ );
int check_zgeqrf        ( CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ );
int check_zgelqf        ( CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ );
int check_zgels         ( cham_trans_t trans, CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB );
int check_zgeqrs        ( cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR );
int check_zgelqs        ( cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR );
int check_zqc           ( cham_side_t side, cham_trans_t trans, CHAM_desc_t *descC, CHAM_desc_t *descQ, CHAM_desc_t *descCC );

#endif

