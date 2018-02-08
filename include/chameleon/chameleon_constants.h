/**
 *
 * @file chameleon_constants.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global constants
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */
#ifndef _chameleon_constants_h_
#define _chameleon_constants_h_

/**
 *
 * @brief Chameleon constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 */
#define ChameleonByte              0
#define ChameleonInteger           1
#define ChameleonRealFloat         2
#define ChameleonRealDouble        3
#define ChameleonComplexFloat      4
#define ChameleonComplexDouble     5

#define ChameleonCM              101
#define ChameleonRM              102
#define ChameleonCCRB            103
#define ChameleonCRRB            104
#define ChameleonRCRB            105
#define ChameleonRRRB            106

#define ChameleonNoTrans         111
#define ChameleonTrans           112
#define ChameleonConjTrans       113

#define ChameleonUpper           121
#define ChameleonLower           122
#define ChameleonUpperLower      123

#define ChameleonNonUnit         131
#define ChameleonUnit            132

#define ChameleonLeft            141
#define ChameleonRight           142

#define ChameleonOneNorm         171
#define ChameleonRealOneNorm     172
#define ChameleonTwoNorm         173
#define ChameleonFrobeniusNorm   174
#define ChameleonInfNorm         175
#define ChameleonRealInfNorm     176
#define ChameleonMaxNorm         177
#define ChameleonRealMaxNorm     178

#define ChameleonDistUniform     201
#define ChameleonDistSymmetric   202
#define ChameleonDistNormal      203

#define ChameleonHermGeev        241
#define ChameleonHermPoev        242
#define ChameleonNonsymPosv      243
#define ChameleonSymPosv         244

#define ChameleonNoPacking       291
#define ChameleonPackSubdiag     292
#define ChameleonPackSupdiag     293
#define ChameleonPackColumn      294
#define ChameleonPackRow         295
#define ChameleonPackLowerBand   296
#define ChameleonPackUpeprBand   297
#define ChameleonPackAll         298

#define ChameleonNoVec           301
#define ChameleonVec             302
#define ChameleonIvec            303

#define ChameleonForward         391
#define ChameleonBackward        392

#define ChameleonColumnwise      401
#define ChameleonRowwise         402
#define ChameleonTrd            1001
#define ChameleonBrd            1002

#define ChameleonW               501
#define ChameleonA2              502

/**
 *  CHAMELEON constants - boolean
 */
#define CHAMELEON_FALSE  0
#define CHAMELEON_TRUE   1

#define CHAMELEON_CPU    ((1ULL)<<1)
#define CHAMELEON_CUDA   ((1ULL)<<3)

/**
 *  State machine switches
 */
#define CHAMELEON_WARNINGS        1
#define CHAMELEON_ERRORS          2
#define CHAMELEON_AUTOTUNING      3
#define CHAMELEON_DAG             4
#define CHAMELEON_PROFILING_MODE  5
#define CHAMELEON_PARALLEL_MODE   6
#define CHAMELEON_BOUND           7
#define CHAMELEON_PROGRESS        8
#define CHAMELEON_GEMM3M          9

/**
 *  CHAMELEON constants - configuration parameters
 */
#define CHAMELEON_CONCURRENCY       1
#define CHAMELEON_TILE_SIZE         2
#define CHAMELEON_INNER_BLOCK_SIZE  3
#define CHAMELEON_HOUSEHOLDER_MODE  5
#define CHAMELEON_HOUSEHOLDER_SIZE  6
#define CHAMELEON_TRANSLATION_MODE  7

#define CHAMELEON_FLAT_HOUSEHOLDER  1
#define CHAMELEON_TREE_HOUSEHOLDER  2

#define CHAMELEON_INPLACE           1
#define CHAMELEON_OUTOFPLACE        2

/**
 *  CHAMELEON constants - success & error codes
 */
#define CHAMELEON_SUCCESS                 0
#define CHAMELEON_ERR_NOT_INITIALIZED  -101
#define CHAMELEON_ERR_REINITIALIZED    -102
#define CHAMELEON_ERR_NOT_SUPPORTED    -103
#define CHAMELEON_ERR_ILLEGAL_VALUE    -104
#define CHAMELEON_ERR_NOT_FOUND        -105
#define CHAMELEON_ERR_OUT_OF_RESOURCES -106
#define CHAMELEON_ERR_INTERNAL_LIMIT   -107
#define CHAMELEON_ERR_UNALLOCATED      -108
#define CHAMELEON_ERR_FILESYSTEM       -109
#define CHAMELEON_ERR_UNEXPECTED       -110
#define CHAMELEON_ERR_SEQUENCE_FLUSHED -111

#endif
