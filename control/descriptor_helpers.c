/**
 *
 * @file descriptor_helpers.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @date 2022-02-22
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

/**
 * @brief Return the rank of the tile A( m, n ) in a classic 2D Block Cyclic
 * distribution PxQ.
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The rank of the tile A( m, n )
 *
 */
int chameleon_getrankof_2d( const CHAM_desc_t *A, int m, int n )
{
    int mm = m + A->i / A->mb;
    int nn = n + A->j / A->nb;
    return (mm % A->p) * A->q + (nn % A->q);
}

/**
 * @brief Return the rank associated to the diagonal tile ( m, m ) of a classic
 * 2D Block Cyclic distribution PxQ.
 *
 * @param[in] A
 *        A specific matrix descriptor that holds only diagonal tiles. Thus, n is never used.
 *
 * @param[in] m
 *        The row and column index of the tile.
 *
 * @param[in] n
 *        Unused
 *
 * @return The rank of the tile A( m, m )
 *
 */
int chameleon_getrankof_2d_diag( const CHAM_desc_t *A, int m, int n )
{
    int mm = m + A->i / A->mb;
    (void)n;
    return (mm % A->p) * A->q + (mm % A->q);
}

/**
 * @brief Return the address of the tile A( m, n ) in a tile storage.
 *
 * @WARNING The only mapping valid with this data storage is the 2D block cyclic
 * storage as soon as multiple nodes are involved (see chameleon_getrankof_2d()).
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The address of the tile A( m, n )
 *
 */
void* chameleon_getaddr_ccrb( const CHAM_desc_t *A, int m, int n )
{
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = CHAMELEON_Element_Size(A->dtyp);
    size_t offset  = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn ) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    if (mm < (size_t)(A->llm1)) {
        if (nn < (size_t)(A->lln1))
            offset = (size_t)(A->bsiz) * (mm + (size_t)(A->llm1) * nn );
        else
            offset = A->A12 + ((size_t)(A->mb * (A->lln%A->nb)) * mm );
    }
    else {
        if (nn < (size_t)(A->lln1))
            offset = A->A21 + ((size_t)((A->llm%A->mb) * A->nb) * nn );
        else
            offset = A->A22;
    }

    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/**
 * @brief Return the address of the tile A( m, n ) in a column major storage.
 *
 * @WARNING The only mapping valid with this data storage is the 2D block cyclic
 * storage as soon as multiple nodes are involved (see chameleon_getrankof_2d()).
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The address of the tile A( m, n )
 *
 */
void *chameleon_getaddr_cm( const CHAM_desc_t *A, int m, int n )
{
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = CHAMELEON_Element_Size(A->dtyp);
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn ) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    offset = (size_t)(A->llm * A->nb) * nn + (size_t)(A->mb) * mm;
    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/**
 * @brief Return the address of the tile A( m, m ) in a tile storage.
 *
 * @WARNING The only mapping valid with this data storage is the diagonal 2D block cyclic
 * storage as soon as multiple nodes are involved (see chameleon_getrankof_2d_diag()).
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row and column index of the tile.
 *
 * @param[in] n
 *        Unused.
 *
 * @return The address of the tile A( m, n )
 *
 */
void *chameleon_getaddr_diag( const CHAM_desc_t *A, int m, int n )
{
    (void)n;
    return chameleon_getaddr_ccrb( A, m, 0 );
}

/**
 * @brief Return the address of the tile A( m, m ) in a dynamic storage.
 *
 * Template that returns a NULL address for all tiles. This is used jointly with
 * the StarPU runtime system to allocate data on the fly and it can be used with
 * any data mapping.
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile (unused).
 *
 * @param[in] m
 *        The row index of the tile (unused).
 *
 * @param[in] n
 *        The column index of the tile (unused).
 *
 * @return A null pointer
 *
 */
void *chameleon_getaddr_null( const CHAM_desc_t *A, int m, int n )
{
    (void)A; (void)m; (void)n;
    return NULL;
}

/**
 * @brief Return the leading dimension of the tile A( m, m ) stored in a tiled storage.
 *
 * This functions returns the leading dimension of the tile A( m, n ) and since
 * the tiled storage is compatible with any mapping strategy, this function is
 * also valid for any mapping.
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The leading dimension of the tile A( m, n ).
 *
 */
int chameleon_getblkldd_ccrb( const CHAM_desc_t *A, int m )
{
    int mm = m + A->i / A->mb;
    return ( ((mm+1) == A->lmt) && ((A->lm % A->mb) != 0)) ? A->lm % A->mb : A->mb;
}

/**
 * @brief Return the leading dimension of the tile A( m, m ) stored in a column major storage.
 *
 * This functions returns the leading dimension of the tile A( m, n ) gfor
 * column major storage. It can be only be used jointly with
 * chameleon_getaddr_cm() and chameleon_getrankof_2d().
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The leading dimension of the tile A( m, n ).
 *
 */
int chameleon_getblkldd_cm( const CHAM_desc_t *A, int m )
{
    (void)m;
    return A->llm;
}
