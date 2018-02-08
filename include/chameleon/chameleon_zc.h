/**
 *
 * @file chameleon_zc.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon mixed precision wrappers header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#ifndef _CHAMELEON_ZC_H_
#define _CHAMELEON_ZC_H_

BEGIN_C_DECLS

/**
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 */
//int CHAMELEON_zcgesv(int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, int * IPIV, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t *X, int LDX, int *ITER);
//int CHAMELEON_zcposv(CHAMELEON_enum uplo, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t *X, int LDX, int *ITER);
//int CHAMELEON_zcgels(CHAMELEON_enum trans, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t *X, int LDX, int *ITER);
//int CHAMELEON_zcungesv(CHAMELEON_enum trans, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t *X, int LDX, int *ITER);

/**
 *  Declarations of math functions (tile layout) - alphabetical order
 */
//int CHAMELEON_zcgesv_Tile(CHAMELEON_desc_t *A, int *IPIV, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER);
//int CHAMELEON_zcposv_Tile(CHAMELEON_enum uplo, CHAMELEON_desc_t *A, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER);
/* int CHAMELEON_zcgels_Tile(CHAMELEON_enum trans, CHAMELEON_desc_t *A, CHAMELEON_desc_t *T, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER); */
//int CHAMELEON_zcungesv_Tile(CHAMELEON_enum trans, CHAMELEON_desc_t *A, CHAMELEON_desc_t *T, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER);

/**
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 */
//int CHAMELEON_zcgesv_Tile_Async(CHAMELEON_desc_t *A, int *IPIV, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zcposv_Tile_Async(CHAMELEON_enum uplo, CHAMELEON_desc_t *A, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zcgels_Tile_Async(CHAMELEON_enum trans, CHAMELEON_desc_t *A, CHAMELEON_desc_t *T, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zcungesv_Tile_Async(CHAMELEON_enum trans, CHAMELEON_desc_t *A, CHAMELEON_desc_t *T, CHAMELEON_desc_t *B, CHAMELEON_desc_t *X, int *ITER, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

END_C_DECLS

#endif
