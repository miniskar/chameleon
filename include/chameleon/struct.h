/**
 *
 * @file chameleon_struct.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon structures
 *
 * @version 1.3.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Matthieu Kuhn
 * @author Lionel Eyraud-Dubois
 * @date 2024-03-16
 *
 */
#ifndef _chameleon_struct_h_
#define _chameleon_struct_h_

#include "chameleon/config.h"
#include "chameleon/types.h"
#include "chameleon/constants.h"
#include "chameleon/runtime_struct.h"

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#else
#ifndef MPI_Comm
typedef uintptr_t MPI_Comm;
#endif
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif
#endif

BEGIN_C_DECLS

#define CHAMELEON_TILE_FULLRANK (1 << 0)
#define CHAMELEON_TILE_DESC     (1 << 1)
#define CHAMELEON_TILE_HMAT     (1 << 2)

/**
 * @brief CHAMELEON structure to hold pivot informations for the LU factorization with partial pivoting
 */
typedef struct chameleon_pivot_s {
    int   blkm0;   /**> The row index of the first row in the tile where the pivot has been selected */
    int   blkidx;  /**> The relative row index in the tile where the pivot has been selected         */
    void *pivrow;  /**> The copy of the row with the selected pivot                                  */
    void *diagrow; /**> The copy of the diagonal row to permute                                      */
} CHAM_pivot_t;

typedef struct chameleon_tile_s {
#if defined(CHAMELEON_KERNELS_TRACE)
    char  *name;
#endif
    void  *mat;
    int    rank, m, n, ld;
    int8_t format;
    int8_t flttype;
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
    blktile_fct_t   get_blktile;     /**> function to get chameleon tiles address           */
    blkaddr_fct_t   get_blkaddr;     /**> function to get chameleon tiles address           */
    blkldd_fct_t    get_blkldd;      /**> function to get chameleon tiles leading dimension */
    blkrankof_fct_t get_rankof;      /**> function to get chameleon tiles MPI rank          */
    blkrankof_fct_t get_rankof_init; /**> function to get chameleon tiles MPI rank          */

    void* get_rankof_init_arg;
    CHAM_tile_t *tiles;  /**> pointer to the array of tiles descriptors  */
    void *mat;           /**> pointer to the beginning of the matrix     */
    size_t A21;          /**> pointer to the beginning of the matrix A21 */
    size_t A12;          /**> pointer to the beginning of the matrix A12 */
    size_t A22;          /**> pointer to the beginning of the matrix A22 */
    cham_storage_t styp; /**> storage layout of the matrix               */
    cham_flttype_t dtyp; /**> precision of the matrix                    */
    int mb;              /**> number of rows in a tile                   */
    int nb;              /**> number of columns in a tile                */
    int bsiz;            /**> size in elements including padding         */

    /* Matrix sizes in single rows/columns for the full problem */
    int i;            /**> row index to the beginning of the submatrix    */
    int j;            /**> column index to the beginning of the submatrix */
    int m;            /**> number of rows of the submatrix                */
    int n;            /**> number of columns of the submatrix             */
    int lm;  	      /**> number of rows of the entire matrix            */
    int ln;           /**> number of columns of the entire matrix         */

    /* Number of rows/columns tiles for the full problem */
    int mt;           /**> number of tile rows    of the submatrix - derived parameter     */
    int nt;           /**> number of tile columns of the submatrix - derived parameter     */
    int lmt;          /**> number of tile rows    of the entire matrix - derived parameter */
    int lnt;          /**> number of tile columns of the entire matrix - derived parameter */

    /* Distributed case */
    int p;            /**> number of rows of the 2D distribution grid                          */
    int q;            /**> number of columns of the 2D distribution grid                       */
    int llm;          /**> local number of rows         of the full matrix - derived parameter */
    int lln;          /**> local number of columns      of the full matrix - derived parameter */
    int llm1;         /**> local number of tile rows    of the A11  matrix - derived parameter */
    int lln1;         /**> local number of tile columns of the A11  matrix - derived parameter */
    int llmt;         /**> local number of tile rows    of the full matrix - derived parameter */
    int llnt;         /**> local number of tile columns of the full matrix - derived parameter */

    int id;           /**> identification number of the descriptor                            */
    int occurences;   /**> identify main matrix desc (occurances=1) or                        */
                      /**> submatrix desc (occurances>1) to avoid unregistering               */
                      /**> GPU data twice                                                     */
    int use_mat;      /**> 1 if we have a pointer to the overall data mat - else 0            */
    int alloc_mat;    /**> 1 if we handle the allocation of mat - else 0                      */
    int register_mat; /**> 1 if we have to register mat - else 0 (handled by the application) */
    int myrank;       /**> MPI rank of the descriptor                                         */
    int ooc;          /**> 1 if the matrix is not to fit in memory                            */
    int64_t mpitag;   /**> First MPI tag used by the descriptor                               */
    void *schedopt;   /**> scheduler (QUARK|StarPU) specific structure                        */
};

/**
 *  CHAMELEON structure to hold pivot informations for the LU factorization with partial pivoting
 */
typedef struct chameleon_piv_s {
    const CHAM_desc_t *desc;   /**> Reference descriptor to compute data mapping based on diagonal tiles,
                              and get floating reference type                                        */
    int    *data;    /**> Pointer to the data                                                    */
    void   *ipiv;    /**> Opaque array of pointers for the runtimes to handle the ipiv array     */
    void   *nextpiv; /**> Opaque array of pointers for the runtimes to handle the pivot computation structure */
    void   *prevpiv; /**> Opaque array of pointers for the runtimes to handle the pivot computation structure */
    void   *perm;    /**> Opaque array of pointers for the runtimes to handle the temporary permutation array */
    void   *invp;    /**> Opaque array of pointers for the runtimes to handle the temporary inverse permutation array */
    int64_t mpitag_ipiv;    /**> Initial mpi tag values for the ipiv handles    */
    int64_t mpitag_nextpiv; /**> Initial mpi tag values for the nextpiv handles */
    int64_t mpitag_prevpiv; /**> Initial mpi tag values for the prevpiv handles */
    int64_t mpitag_perm;    /**> Initial mpi tag values for the nextpiv handles */
    int64_t mpitag_invp;    /**> Initial mpi tag values for the prevpiv handles */
    int     i;              /**> row index to the beginning of the submatrix    */
    int     m;              /**> The number of row in the vector ipiv           */
    int     mb;             /**> The number of row per block                    */
    int     mt;             /**> The number of tiles                            */
    int     n;              /**> The number of column considered (must be updated for each panel) */
} CHAM_ipiv_t;

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
    cham_bool_t        statistics_enabled;
    cham_bool_t        progress_enabled;
    cham_bool_t        generic_enabled;
    cham_bool_t        autominmax_enabled;
    cham_bool_t        runtime_paused;

    cham_householder_t householder;        // "domino" (flat) or tree-based (reduction) Householder
    cham_translation_t translation;        // In place or Out of place layout conversion

    int                nb;
    int                ib;
    int                rhblock;            // block size for tree-based (reduction) Householder
    int                lookahead;          // depth of the look ahead in algorithms
    void              *schedopt;           // structure for runtimes
    int                mpi_outer_init;     // MPI has been initialized outside our functions
    MPI_Comm           comm;               // MPI communicator
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
