/**
 *
 * @file descriptor.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @date 2021-03-16
 *
 ***
 *
 * @defgroup Descriptor
 * @brief Group descriptor routines exposed to users
 *
 */
#define _GNU_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "control/common.h"
#include "control/descriptor.h"
#include "chameleon/runtime.h"

static int nbdesc = 0;

/**
 *
 */
int chameleon_desc_mat_alloc( CHAM_desc_t *desc )
{
    size_t size = (size_t)(desc->llm) * (size_t)(desc->lln)
        * (size_t)CHAMELEON_Element_Size(desc->dtyp);
    if ((desc->mat = RUNTIME_malloc(size)) == NULL) {
        chameleon_error("chameleon_desc_mat_alloc", "malloc() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    /* The matrix has already been registered by the Runtime alloc */
    desc->register_mat = 0;

    return CHAMELEON_SUCCESS;
}

/**
 *
 */
int chameleon_desc_mat_free( CHAM_desc_t *desc )
{
    if ( (desc->mat       != NULL) &&
         (desc->use_mat   == 1   ) &&
         (desc->alloc_mat == 1   ) )
    {
        size_t size = (size_t)(desc->llm) * (size_t)(desc->lln)
            * (size_t)CHAMELEON_Element_Size(desc->dtyp);

        RUNTIME_free(desc->mat, size);
        desc->mat = NULL;
    }

    if ( desc->tiles ) {
#if defined(CHAMELEON_KERNELS_TRACE)
        CHAM_tile_t *tile = desc->tiles;
        int ii, jj;
        for( jj=0; jj<desc->lnt; jj++ ) {
            for( ii=0; ii<desc->lmt; ii++, tile++ ) {
                if ( tile->name ) {
                    free( tile->name );
                }
            }
        }
#endif
        free( desc->tiles );
    }
    return CHAMELEON_SUCCESS;
}

void chameleon_desc_init_tiles( CHAM_desc_t *desc )
{
    CHAM_tile_t *tile;
    int ii, jj;

    desc->tiles = malloc( desc->lmt * desc->lnt * sizeof(CHAM_tile_t) );

    tile = desc->tiles;
    for( jj=0; jj<desc->lnt; jj++ ) {
        for( ii=0; ii<desc->lmt; ii++, tile++ ) {
            int rank = desc->get_rankof( desc, ii, jj );
            tile->format = CHAMELEON_TILE_FULLRANK;
            tile->m   = ii == desc->lmt-1 ? desc->lm - ii * desc->mb : desc->mb;
            tile->n   = jj == desc->lnt-1 ? desc->ln - jj * desc->nb : desc->nb;
            tile->mat = (rank == desc->myrank) ? desc->get_blkaddr( desc, ii, jj ) : NULL;
            tile->ld  = desc->get_blkldd( desc, ii );
#if defined(CHAMELEON_KERNELS_TRACE)
            asprintf( &(tile->name), "%s(%d,%d)", desc->name, ii, jj);
#endif
        }
    }
}

/**
 *  Internal function to return MPI rank of element A(m,n) with m,n = block indices
 */
int chameleon_getrankof_2d( const CHAM_desc_t *A, int m, int n )
{
    int mm = m + A->i / A->mb;
    int nn = n + A->j / A->nb;
    return (mm % A->p) * A->q + (nn % A->q);
}

/**
 *  Internal function to return MPI rank of element DIAG(m,0) with m,n = block indices
 */
int chameleon_getrankof_2d_diag( const CHAM_desc_t *A, int m, int n )
{
    int mm = m + A->i / A->mb;
    assert( m == n );
    return (mm % A->p) * A->q + (mm % A->q);
}


/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to create tiled matrix descriptor
 * with generic function for data distribution and storage format.
 *
 ******************************************************************************
 *
 * @param[in] name
 *          Name of the descriptor for debug purpose.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg ChamRealFloat:     single precision real (S),
 *          @arg ChamRealDouble:    double precision real (D),
 *          @arg ChamComplexFloat:  single precision complex (C),
 *          @arg ChamComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 * @param[in] get_blkaddr
 *          A function which return the address of the data corresponding to
 *          the tile A(m,n).
 *
 * @param[in] get_blkldd
 *          A function that return the leading dimension of the tile A(m,*).
 *
 * @param[in] get_rankof
 *          A function that return the MPI rank of the tile A(m,n).
 *
 ******************************************************************************
 *
 * @return  The descriptor with the matrix description parameters set.
 *
 */
int chameleon_desc_init_internal( CHAM_desc_t *desc, const char *name, void *mat,
                                  cham_flttype_t dtyp, int mb, int nb,
                                  int lm, int ln, int m, int n, int p, int q,
                                  blkaddr_fct_t   get_blkaddr,
                                  blkldd_fct_t    get_blkldd,
                                  blkrankof_fct_t get_rankof )
{
    CHAM_context_t *chamctxt;
    int rc = CHAMELEON_SUCCESS;

    memset( desc, 0, sizeof(CHAM_desc_t) );

    desc->name = name;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* If one of the function get_* is NULL, we switch back to the default */
    desc->get_blktile = chameleon_desc_gettile;
    desc->get_blkaddr = get_blkaddr ? get_blkaddr : chameleon_getaddr_ccrb;
    desc->get_blkldd  = get_blkldd  ? get_blkldd  : chameleon_getblkldd_ccrb;
    desc->get_rankof  = get_rankof  ? get_rankof  : chameleon_getrankof_2d;

    /* Matrix properties */
    desc->dtyp = dtyp;
    /* Should be given as parameter to follow get_blkaddr (unused) */
    desc->styp = (get_blkaddr == chameleon_getaddr_cm ) ? ChamCM : ChamCCRB;
    desc->mb   = mb;
    desc->nb   = nb;
    desc->bsiz = mb * nb;

    /* Matrix parameters */
    desc->i = 0;
    desc->j = 0;
    desc->m = m;
    desc->n = n;

    /* Matrix stride parameters */
    desc->lm = m;
    desc->ln = n;

    /* Matrix derived parameters */
    desc->mt  = chameleon_ceil( m, mb );
    desc->nt  = chameleon_ceil( n, nb );
    desc->lmt = desc->mt;
    desc->lnt = desc->nt;

    desc->id = nbdesc;
    nbdesc++;
    desc->occurences = 0;

    desc->myrank = RUNTIME_comm_rank( chamctxt );

    /* Grid size */
    desc->p = p;
    desc->q = q;

    /* Local dimensions in tiles */
    if ( desc->myrank < (p*q) ) {
        int gmt, gnt;

        /* Compute the fictive full number of tiles to derivate the local leading dimension */
        gmt = chameleon_ceil( lm, mb );
        gnt = chameleon_ceil( ln, nb );

        desc->llmt = chameleon_ceil( gmt, p );
        desc->llnt = chameleon_ceil( gnt, q );

        // Local dimensions
        if ( ((desc->lmt-1) % p) == (desc->myrank / q) ) {
            desc->llm  = ( desc->llmt - 1 ) * mb + ((lm%mb==0) ? mb : (lm%mb));
        } else {
            desc->llm  =  desc->llmt * mb;
        }

        if ( ((desc->lnt-1) % q) == (desc->myrank % q) ) {
            desc->lln  = ( desc->llnt - 1 ) * nb + ((ln%nb==0) ? nb : (ln%nb));
        } else {
            desc->lln  =  desc->llnt * nb;
        }

        desc->llm1 = desc->llm / mb;
        desc->lln1 = desc->lln / nb;
    } else {
        desc->llmt = 0;
        desc->llnt = 0;
        desc->llm  = 0;
        desc->lln  = 0;
        desc->llm1 = 0;
        desc->lln1 = 0;
    }

    /* memory of the matrix is handled by the user */
    desc->alloc_mat    = 0;
    /* if the user gives a pointer to the overall data (tiles) we can use it */
    desc->use_mat      = 0;
    /* users data can have multiple forms: let him register tiles */
    desc->register_mat = 0;
    /* The matrix is alocated tile by tile with out of core */
    desc->ooc = 0;

    switch ( (intptr_t)mat ) {
    case (intptr_t)CHAMELEON_MAT_ALLOC_TILE:
        if ( chamctxt->scheduler == RUNTIME_SCHED_STARPU ) {
            /* Let's use the allocation on the fly as in OOC */
            desc->get_blkaddr = chameleon_getaddr_null;
            desc->mat = NULL;
            break;
        }
        /* Otherwise we switch back to the full allocation */
        chameleon_attr_fallthrough;

    case (intptr_t)CHAMELEON_MAT_ALLOC_GLOBAL:
        rc = chameleon_desc_mat_alloc( desc );
        desc->alloc_mat = 1;
        desc->use_mat   = 1;
        break;

    case (intptr_t)CHAMELEON_MAT_OOC:
        if ( chamctxt->scheduler != RUNTIME_SCHED_STARPU ) {
            chameleon_error("CHAMELEON_Desc_Create", "CHAMELEON Out-of-Core descriptors are supported only with StarPU");
            return CHAMELEON_ERR_NOT_SUPPORTED;
        }
        desc->mat = NULL;
        desc->ooc = 1;
        break;

    default:
        /* memory of the matrix is handled by users */
        desc->mat     = mat;
        desc->use_mat = 1;
    }

    desc->A21 = (size_t)(desc->llm - desc->llm%mb)*(size_t)(desc->lln - desc->lln%nb);
    desc->A12 = (size_t)(            desc->llm%mb)*(size_t)(desc->lln - desc->lln%nb) + desc->A21;
    desc->A22 = (size_t)(desc->llm - desc->llm%mb)*(size_t)(            desc->lln%nb) + desc->A12;

    chameleon_desc_init_tiles( desc );

    /* Create runtime specific structure like registering data */
    RUNTIME_desc_create( desc );

    return rc;
}

/**
 *  Internal static descriptor initializer for submatrices
 */
CHAM_desc_t* chameleon_desc_submatrix( CHAM_desc_t *descA, int i, int j, int m, int n )
{
    CHAM_desc_t *descB = malloc(sizeof(CHAM_desc_t));
    int mb, nb;

    if ( (descA->i + i + m) > descA->m ) {
        chameleon_error("chameleon_desc_submatrix", "The number of rows (i+m) of the submatrix doesn't fit in the parent matrix");
        assert((descA->i + i + m) > descA->m);
    }
    if ( (descA->j + j + n) > descA->n ) {
        chameleon_error("chameleon_desc_submatrix", "The number of rows (j+n) of the submatrix doesn't fit in the parent matrix");
        assert((descA->j + j + n) > descA->n);
    }

    memcpy( descB, descA, sizeof(CHAM_desc_t) );
    mb = descA->mb;
    nb = descA->nb;
    // Submatrix parameters
    descB->i = descA->i + i;
    descB->j = descA->j + j;
    descB->m = m;
    descB->n = n;
    // Submatrix derived parameters
    descB->mt = (m == 0) ? 0 : (descB->i+m-1)/mb - descB->i/mb + 1;
    descB->nt = (n == 0) ? 0 : (descB->j+n-1)/nb - descB->j/nb + 1;

    // Increase the number of occurences to avoid multiple free of runtime specific data structures.
    descB->occurences++;

    return descB;
}

void chameleon_desc_destroy( CHAM_desc_t *desc )
{
    RUNTIME_desc_destroy( desc );
    chameleon_desc_mat_free( desc );
}

/**
 *  Check for descriptor correctness
 */
int chameleon_desc_check(const CHAM_desc_t *desc)
{
    if (desc == NULL) {
        chameleon_error("chameleon_desc_check", "NULL descriptor");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (desc->mat == NULL && desc->use_mat == 1) {
        chameleon_error("chameleon_desc_check", "NULL matrix pointer");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (desc->dtyp != ChamRealFloat &&
        desc->dtyp != ChamRealDouble &&
        desc->dtyp != ChamComplexFloat &&
        desc->dtyp != ChamComplexDouble  ) {
        chameleon_error("chameleon_desc_check", "invalid matrix type");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if (desc->mb <= 0 || desc->nb <= 0) {
        chameleon_error("chameleon_desc_check", "negative tile dimension");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if (desc->bsiz < desc->mb*desc->nb) {
        chameleon_error("chameleon_desc_check", "tile memory size smaller than the product of dimensions");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if (desc->lm <= 0 || desc->ln <= 0) {
        chameleon_error("chameleon_desc_check", "negative matrix dimension");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ((desc->lm < desc->m) || (desc->ln < desc->n)) {
        chameleon_error("chameleon_desc_check", "matrix dimensions larger than leading dimensions");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ((desc->i > 0 && desc->i >= desc->lm) || (desc->j > 0 && desc->j >= desc->ln)) {
        chameleon_error("chameleon_desc_check", "beginning of the matrix out of scope");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if (desc->i+desc->m > desc->lm || desc->j+desc->n > desc->ln) {
        chameleon_error("chameleon_desc_check", "submatrix out of scope");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    return CHAMELEON_SUCCESS;
}

CHAM_desc_t *
CHAMELEON_Desc_SubMatrix( CHAM_desc_t *descA, int i, int j, int m, int n )
{
    return chameleon_desc_submatrix( descA, i, j, m, n );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Create - Create tiled matrix descriptor.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg ChamRealFloat:     single precision real (S),
 *          @arg ChamRealDouble:    double precision real (D),
 *          @arg ChamComplexFloat:  single precision complex (C),
 *          @arg ChamComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in number of elements of each tile, including internal padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Create( CHAM_desc_t **descptr, void *mat, cham_flttype_t dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n, int p, int q )
{
    return CHAMELEON_Desc_Create_User( descptr, mat, dtyp, mb, nb, bsiz,
                                       lm, ln, i, j, m, n, p, q,
                                       NULL, NULL, NULL );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Create_User - Create generic tiled matrix descriptor for general
 *  applications.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg ChamRealFloat:     single precision real (S),
 *          @arg ChamRealDouble:    double precision real (D),
 *          @arg ChamComplexFloat:  single precision complex (C),
 *          @arg ChamComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in number of elements of each tile, including internal padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 * @param[in] get_blkaddr
 *          A function which return the address of the data corresponding to
 *          the tile A(m,n).
 *
 * @param[in] get_blkldd
 *          A function that return the leading dimension of the tile A(m,*).
 *
 * @param[in] get_rankof
 *          A function that return the MPI rank of the tile A(m,n).
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Create_User( CHAM_desc_t **descptr, void *mat, cham_flttype_t dtyp, int mb, int nb, int bsiz,
                                int lm, int ln, int i, int j, int m, int n, int p, int q,
                                blkaddr_fct_t   get_blkaddr,
                                blkldd_fct_t    get_blkldd,
                                blkrankof_fct_t get_rankof )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t *desc;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Create_User", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Allocate memory and initialize the descriptor */
    desc = (CHAM_desc_t*)malloc(sizeof(CHAM_desc_t));
    if (desc == NULL) {
        chameleon_error("CHAMELEON_Desc_Create_User", "malloc() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    chameleon_desc_init( desc, mat, dtyp, mb, nb, bsiz,
                         lm, ln, i, j, m, n, p, q,
                         get_blkaddr, get_blkldd, get_rankof );

    status = chameleon_desc_check( desc );
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_Desc_Create_User", "invalid descriptor");
        CHAMELEON_Desc_Destroy( &desc );
        return status;
    }

    *descptr = desc;
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Create_OOC_User - Create matrix descriptor for tiled matrix which
 *  may not fit memory.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg ChamRealFloat:     single precision real (S),
 *          @arg ChamRealDouble:    double precision real (D),
 *          @arg ChamComplexFloat:  single precision complex (C),
 *          @arg ChamComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in number of elements of each tile, including internal padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 * @param[in] get_rankof
 *          A function that return the MPI rank of the tile A(m,n).
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Create_OOC_User(CHAM_desc_t **descptr, cham_flttype_t dtyp, int mb, int nb, int bsiz,
                                   int lm, int ln, int i, int j, int m, int n, int p, int q,
                                   blkrankof_fct_t get_rankof )
{
#if !defined (CHAMELEON_SCHED_STARPU)
    (void)descptr; (void)dtyp; (void)mb; (void)nb; (void)bsiz;
    (void)lm; (void)ln; (void)i; (void)j; (void)m; (void)n; (void)p; (void)q;
    (void)get_rankof;

    chameleon_error("CHAMELEON_Desc_Create_OOC_User", "Only StarPU supports on-demand tile allocation");
    return CHAMELEON_ERR_NOT_SUPPORTED;
#else
    int rc;
    rc = CHAMELEON_Desc_Create_User( descptr, CHAMELEON_MAT_OOC, dtyp, mb, nb, bsiz,
                                     lm, ln, i, j, m, n, p, q,
                                     chameleon_getaddr_null, NULL, get_rankof );
    return rc;
#endif
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Create_OOC - Create matrix descriptor for tiled matrix which may
 *  not fit memory.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg ChamRealFloat:     single precision real (S),
 *          @arg ChamRealDouble:    double precision real (D),
 *          @arg ChamComplexFloat:  single precision complex (C),
 *          @arg ChamComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in number of elements of each tile, including internal padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Create_OOC(CHAM_desc_t **descptr, cham_flttype_t dtyp, int mb, int nb, int bsiz,
                              int lm, int ln, int i, int j, int m, int n, int p, int q)
{
    return CHAMELEON_Desc_Create_User( descptr, CHAMELEON_MAT_OOC, dtyp, mb, nb, bsiz,
                                       lm, ln, i, j, m, n, p, q,
                                       chameleon_getaddr_null, NULL, NULL );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Creates a new descriptor with the same properties as the one given as
 * input.
 *
 * @warning This function copies the descriptor structure, but does not copy the
 * matrix data.
 *
 ******************************************************************************
 *
 * @param[in] descin
 *          The descriptor structure to duplicate.
 *
 * @param[in] mat
 *          Memory location for the copy. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 ******************************************************************************
 *
 * @retval The new matrix descriptor.
 *
 */
CHAM_desc_t *CHAMELEON_Desc_Copy( const CHAM_desc_t *descin, void *mat )
{
    CHAM_desc_t *descout = NULL;
    CHAMELEON_Desc_Create_User( &descout, mat,
                                descin->dtyp, descin->mb, descin->nb, descin->bsiz,
                                descin->lm, descin->ln, descin->i, descin->j, descin->m, descin->n, descin->p, descin->q,
                                NULL, NULL, descin->get_rankof );
    return descout;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Creates a new descriptor with the same properties as the one given as
 * input and restricted on node 0.
 *
 * @warning This function copies the descriptor structure, but does not copy the
 * matrix data.
 *
 ******************************************************************************
 *
 * @param[in] descin
 *          The descriptor structure to duplicate.
 *
 * @param[in] mat
 *          Memory location for the copy. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 ******************************************************************************
 *
 * @retval The new matrix descriptor.
 *
 */
CHAM_desc_t *CHAMELEON_Desc_CopyOnZero( const CHAM_desc_t *descin, void *mat )
{
    CHAM_desc_t *descout = NULL;
    CHAMELEON_Desc_Create_User( &descout, mat,
                                descin->dtyp, descin->mb, descin->nb, descin->bsiz,
                                descin->lm, descin->ln, descin->i, descin->j, descin->m, descin->n, 1, 1,
                                NULL, NULL, descin->get_rankof );
    return descout;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Destroy - Destroys matrix descriptor.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Destroy(CHAM_desc_t **descptr)
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t *desc;
    int m, n;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Destroy", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ((descptr == NULL) || (*descptr == NULL)) {
        chameleon_error("CHAMELEON_Desc_Destroy", "attempting to destroy a NULL descriptor");
        return CHAMELEON_ERR_UNALLOCATED;
    }

    desc = *descptr;
    for ( n=0; n<desc->nt; n++ ) {
        for ( m=0; m<desc->mt; m++ ) {
            CHAM_tile_t *tile;

            tile = desc->get_blktile( desc, m, n );

            if ( tile->format == CHAMELEON_TILE_DESC ) {
                CHAM_desc_t *tiledesc = tile->mat;

                /* Recursive names are allocated internally, we need to free them */
                if ( tiledesc->name ) {
                    free( (void*)(tiledesc->name) );
                }
                CHAMELEON_Desc_Destroy( &tiledesc );
                assert( tiledesc == NULL );
            }
        }
    }

    chameleon_desc_destroy( desc );
    free(desc);
    *descptr = NULL;
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Acquire - Ensures that all data of the descriptor are
 *  up-to-date.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Acquire( const CHAM_desc_t *desc ) {
    return RUNTIME_desc_acquire( desc );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Release - Release the data of the descriptor acquired by the
 *  application. Should be called if CHAMELEON_Desc_Acquire has been called on the
 *  descriptor and if you do not need to access to its data anymore.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Release( const CHAM_desc_t *desc ) {
    return RUNTIME_desc_release( desc );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_Desc_Flush - Flushes the data in the sequence when they won't be
 *  reused. This calls cleans up the distributed communication caches, and
 *  transfer the data back to the CPU.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc_Flush( const CHAM_desc_t        *desc,
                          const RUNTIME_sequence_t *sequence )
{
    RUNTIME_desc_flush( desc, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  CHAMELEON_user_tag_size - Set the sizes for the MPI tags
 *  Default value: tag_width=31, tag_sep=24, meaning that the MPI tag is stored
 *  in 31 bits, with 24 bits for the tile tag and 7 for the descriptor.  This
 *  function must be called before any descriptor creation.
 *
 ******************************************************************************
 *
 * @param[in] user_tag_width
 *          The new value for tag_width.
 *
 * @param[in] user_tag_sep
 *          The new value for tag_sep.
 *
 */
void CHAMELEON_user_tag_size(int user_tag_width, int user_tag_sep) {
    RUNTIME_comm_set_tag_sizes( user_tag_width, user_tag_sep );
    return;
}

static void
chameleon_desc_print( const CHAM_desc_t *desc, int shift )
{
    intptr_t base = (intptr_t)desc->mat;
    int m, n, rank;
    CHAM_context_t *chamctxt = chameleon_context_self();

    rank = CHAMELEON_Comm_rank();

    for ( n=0; n<desc->nt; n++ ) {
        for ( m=0; m<desc->mt; m++ ) {
            const CHAM_tile_t *tile;
            const CHAM_desc_t *tiledesc;
            intptr_t ptr;
            int      trank;

            trank    = desc->get_rankof( desc, m, n );
            tile     = desc->get_blktile( desc, m, n );
            tiledesc = tile->mat;

            ptr = ( tile->format == CHAMELEON_TILE_DESC ) ? (intptr_t)(tiledesc->mat) : (intptr_t)(tile->mat);

            if ( trank == rank ) {
                fprintf( stdout, "[%2d]%*s%s(%3d,%3d): %d * %d / ld = %d / offset= %ld\n",
                         rank, shift, " ", desc->name, m, n, tile->m, tile->n, tile->ld, ptr - base );

                if ( tile->format == CHAMELEON_TILE_DESC ) {
                    chameleon_desc_print( tiledesc, shift+2 );
                }
            }
            else {
                assert( ptr == 0 );
            }

            if ( chamctxt->scheduler != RUNTIME_SCHED_OPENMP ) {
                RUNTIME_barrier(chamctxt);
            }
        }
    }
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  @brief Print descriptor structure for debug purpose
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          The input desc for which to describe to print the tile structure
 */
void
CHAMELEON_Desc_Print( const CHAM_desc_t *desc )
{
    chameleon_desc_print( desc, 2 );
}
