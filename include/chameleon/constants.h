/**
 *
 * @file chameleon_constants.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2020-03-03
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
/**
 * @brief Matrix floating point arithmetic
 */
typedef enum chameleon_flttype_e {
    ChamByte          = 0,
    ChamInteger       = 1,
    ChamRealFloat     = 2,
    ChamRealDouble    = 3,
    ChamComplexFloat  = 4,
    ChamComplexDouble = 5,
} cham_flttype_t;

/**
 * @brief Matrix tile storage
 */
typedef enum chameleon_storage_e {
    ChamCM   = 101,
    ChamRM   = 102,
    ChamCCRB = 103,
    ChamCRRB = 104,
    ChamRCRB = 105,
    ChamRRRB = 106,
} cham_storage_t;

/**
 * @brief Transpostion
 */
typedef enum chameleon_trans_e {
    ChamNoTrans   = 111, /**< Use A         */
    ChamTrans     = 112, /**< Use A^t       */
    ChamConjTrans = 113  /**< Use conj(A^t) */
} cham_trans_t;

static inline int
isValidTrans( cham_trans_t trans )
{
    if ( (trans >= ChamNoTrans) && (trans <= ChamConjTrans) ) {
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * @brief Upper/Lower part
 */
typedef enum chameleon_uplo_e {
    ChamUpper      = 121, /**< Use lower triangle of A */
    ChamLower      = 122, /**< Use upper triangle of A */
    ChamUpperLower = 123  /**< Use the full A          */
} cham_uplo_t;

/**
 * @brief Diagonal
 */
typedef enum chameleon_diag_e {
    ChamNonUnit = 131, /**< Diagonal is non unitary */
    ChamUnit    = 132  /**< Diagonal is unitary     */
} cham_diag_t;

/**
 * @brief Side of the operation
 */
typedef enum chameleon_side_e {
    ChamLeft  = 141, /**< Apply operator on the left  */
    ChamRight = 142  /**< Apply operator on the right */
} cham_side_t;

/**
 * @brief Norms
 */
typedef enum chameleon_normtype_e {
    ChamOneNorm       = 171, /**< One norm:       max_j( sum_i( |a_{ij}| ) )   */
    /* ChamRealOneNorm   = 172, */
    /* ChamTwoNorm       = 173, */
    ChamFrobeniusNorm = 174, /**< Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) ) */
    ChamInfNorm       = 175, /**< Inifinite norm: max_i( sum_j( |a_{ij}| ) )   */
    /* ChamRealInfNorm   = 176, */
    ChamMaxNorm       = 177, /**< Inifinite norm: max_{i,j}( | a_{ij} | )      */
    /* ChamRealMaxNorm   = 178 */
} cham_normtype_t;

/**
 * @brief Random distribution for matrix generator
 */
typedef enum chameleon_dist_e {
    ChamDistUniform   = 201,
    ChamDistSymmetric = 202,
    ChamDistNormal    = 203,
} cham_dist_t;

/**
 * @brief Main matrix type
 */
typedef enum chameleon_mtxtype_e {
    ChamGeneral    = 231,
    ChamSymmetric  = 232,
    ChamHermitian  = 233,
    ChamTriangular = 234,
} cham_mtxtype_t;

/**
 * @brief Eigen and singular values generator format
 */
#define ChamHermGeev        241
#define ChamHermPoev        242
#define ChamNonsymPosv      243
#define ChamSymPosv         244

#define ChamNoPacking       291
#define ChamPackSubdiag     292
#define ChamPackSupdiag     293
#define ChamPackColumn      294
#define ChamPackRow         295
#define ChamPackLowerBand   296
#define ChamPackUpeprBand   297
#define ChamPackAll         298

/**
 * @brief Singular/Eigen vector job description
 */
typedef enum chameleon_job_e {
    ChamNoVec = 301,
    ChamVec   = 302,
    ChamIvec  = 303,
} cham_job_t;

/**
 * @brief Algorithm Direction
 */
typedef enum chameleon_dir_e {
    ChamDirForward  = 391, /**< Forward direction   */
    ChamDirBackward = 392, /**< Backward direction  */
} cham_dir_t;

/**
 * @brief Direction of the main vectors as for the householder reflectors in QR/LQ factorizations.
 */
typedef enum chameleon_store_e {
    ChamColumnwise = 401, /**< Column wise storage  */
    ChamRowwise    = 402, /**< Row wise storage     */
    ChamEltwise    = 403, /**< Element by element storage */
} cham_store_t;


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
#define CHAMELEON_GENERIC        10

/**
 *  CHAMELEON constants - configuration parameters
 */
#define CHAMELEON_CONCURRENCY       1
#define CHAMELEON_TILE_SIZE         2
#define CHAMELEON_INNER_BLOCK_SIZE  3
#define CHAMELEON_HOUSEHOLDER_MODE  5
#define CHAMELEON_HOUSEHOLDER_SIZE  6
#define CHAMELEON_TRANSLATION_MODE  7
#define CHAMELEON_LOOKAHEAD         8

/**
 * @brief QR/LQ factorization trees
 */
typedef enum chameleon_householder_e {
    ChamFlatHouseholder = 1,
    ChamTreeHouseholder = 2,
} cham_householder_t;

/**
 * @brief Translation types
 */
typedef enum chameleon_translation_e {
    ChamInPlace    = 1,
    ChamOutOfPlace = 2,
} cham_translation_t;

/**
 * @brief Constant to describe how to initialize the mat pointer in descriptors
 */
#define CHAMELEON_MAT_ALLOC_GLOBAL NULL
#define CHAMELEON_MAT_ALLOC_TILE   ((void*)-1)
#define CHAMELEON_MAT_OOC          ((void*)-2)

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

#endif /* _chameleon_constants_h_ */
