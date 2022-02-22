/**
 *
 * @file chameleon_struct.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon structures
 *
 * @version 1.1.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2021-03-16
 *
 */
#ifndef _chameleon_struct_h_
#define _chameleon_struct_h_

#include "chameleon/config.h"
#include "chameleon/types.h"
#include "chameleon/constants.h"
#include "chameleon/runtime_struct.h"

BEGIN_C_DECLS

#define CHAMELEON_TILE_FULLRANK (1 << 0)
#define CHAMELEON_TILE_DESC     (1 << 1)
#define CHAMELEON_TILE_HMAT     (1 << 2)

typedef struct chameleon_tile_s {
#if defined(CHAMELEON_KERNELS_TRACE)
    char  *name;
#endif
    int8_t format;
    int    m, n, ld;
    void  *mat;
} CHAM_tile_t;

/**
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containning in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%nb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */
struct chameleon_desc_s;
typedef struct chameleon_desc_s CHAM_desc_t;

typedef void*        (*blkaddr_fct_t)  ( const CHAM_desc_t*, int, int );
typedef int          (*blkldd_fct_t)   ( const CHAM_desc_t*, int );
typedef int          (*blkrankof_fct_t)( const CHAM_desc_t*, int, int );
typedef CHAM_tile_t* (*blktile_fct_t)  ( const CHAM_desc_t*, int, int );

struct chameleon_desc_s {
    const char *name;
    // function to get chameleon tiles address
    blktile_fct_t  get_blktile;
    // function to get chameleon tiles address
    blkaddr_fct_t   get_blkaddr;
    // function to get chameleon tiles leading dimension
    blkldd_fct_t    get_blkldd;
    // function to get chameleon tiles MPI rank
    blkrankof_fct_t get_rankof;
    CHAM_tile_t *tiles; // pointer to the array of tiles descriptors
    void *mat;        // pointer to the beginning of the matrix
    size_t A21;       // pointer to the beginning of the matrix A21
    size_t A12;       // pointer to the beginning of the matrix A12
    size_t A22;       // pointer to the beginning of the matrix A22
    cham_storage_t styp;  // storage layout of the matrix
    cham_flttype_t dtyp;  // precision of the matrix
    int mb;           // number of rows in a tile
    int nb;           // number of columns in a tile
    int bsiz;         // size in elements including padding
    int lm;  	      // number of rows of the entire matrix
    int ln;           // number of columns of the entire matrix
    int lmt;          // number of tile rows of the entire matrix - derived parameter
    int lnt;          // number of tile columns of the entire matrix - derived parameter
    int i;            // row index to the beginning of the submatrix
    int j;            // column index to the beginning of the submatrix
    int m;            // number of rows of the submatrix
    int n;            // number of columns of the submatrix
    int mt;           // number of tile rows of the submatrix - derived parameter
    int nt;           // number of tile columns of the submatrix - derived parameter
                      // Data for distributed cases
    int p;            // number of rows of the 2D distribution grid
    int q;            // number of columns of the 2D distribution grid
    int llm;          // number of rows of the 2D distribution grid
    int lln;          // number of columns of the 2D distribution grid
    int llm1;         // number of tile rows of the A11 matrix - derived parameter
    int lln1;         // number of tile columns of the A11 matrix - derived parameter
    int llmt;         // number of tile rows of the local (to a node) matrix
    int llnt;         // number of tile columns of the local (to a node) matrix
    int id;           // identification number of the descriptor
    int occurences;   // identify main matrix desc (occurances=1) or
                      // submatrix desc (occurances>1) to avoid unregistering
                      // GPU data twice
    int use_mat;      // 1 if we have a pointer to the overall data mat - else 0
    int alloc_mat;    // 1 if we handle the allocation of mat - else 0
    int register_mat; // 1 if we have to register mat - else 0 (handled by the application)
    int myrank;       // MPI rank of the descriptor
    int ooc;          // 1 if the matrix is not to fit in memory
    void *schedopt;   // scheduler (QUARK|StarPU) specific structure
};


/**
 *  CHAMELEON request uniquely identifies each asynchronous function call.
 */
typedef struct chameleon_context_s {
    RUNTIME_id_t       scheduler;
    int                nworkers;
    int                ncudas;
    int                nthreads_per_worker;

    /* Boolean flags */
    cham_bool_t        warnings_enabled;
    cham_bool_t        autotuning_enabled;
    cham_bool_t        parallel_enabled;
    cham_bool_t        profiling_enabled;
    cham_bool_t        progress_enabled;
    cham_bool_t        generic_enabled;
    cham_bool_t        runtime_paused;

    cham_householder_t householder;        // "domino" (flat) or tree-based (reduction) Householder
    cham_translation_t translation;        // In place or Out of place layout conversion

    int                nb;
    int                ib;
    int                rhblock;            // block size for tree-based (reduction) Householder
    int                lookahead;          // depth of the look ahead in algorithms
    void              *schedopt;           // structure for runtimes
    int                mpi_outer_init;     // MPI has been initialized outside our functions
} CHAM_context_t;

static inline void *
CHAM_tile_get_ptr( const CHAM_tile_t *tile )
{
    if ( tile->format & CHAMELEON_TILE_DESC ) {
        return ((CHAM_desc_t*)(tile->mat))->mat;
    }
    return tile->mat;
}

static inline const char *
CHAM_tile_get_typestr( const CHAM_tile_t *tile )
{
    if ( tile->format & CHAMELEON_TILE_DESC ) {
        return "Desc";
    }

    if ( tile->format & CHAMELEON_TILE_HMAT ) {
        return "HMat";
    }

    return "Full";
}

END_C_DECLS

#endif /* _chameleon_struct_h_ */
