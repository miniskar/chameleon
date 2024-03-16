/**
 *
 * @file starpu/runtime_descriptor_ipiv.c
 *
 * @copyright 2022-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU descriptor routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2024-03-16
 *
 */
#include "chameleon_starpu.h"

/**
 *  Create ws_pivot runtime structures
 */
void RUNTIME_ipiv_create( CHAM_ipiv_t *ipiv )
{
    assert( ipiv );
    starpu_data_handle_t *handles = calloc( 5 * ipiv->mt, sizeof(starpu_data_handle_t) );
    ipiv->ipiv    = handles;
    handles += ipiv->mt;
    ipiv->nextpiv = handles;
    handles += ipiv->mt;
    ipiv->prevpiv = handles;
    handles += ipiv->mt;
    ipiv->perm    = handles;
    handles += ipiv->mt;
    ipiv->invp    = handles;
#if defined(CHAMELEON_USE_MPI)
    /*
     * Book the number of tags required to describe pivot structure
     * One per handle type
     */
    {
        chameleon_starpu_tag_init();
        ipiv->mpitag_ipiv = chameleon_starpu_tag_book( (int64_t)(ipiv->mt) * 5 );
        if ( ipiv->mpitag_ipiv == -1 ) {
            chameleon_fatal_error("RUNTIME_ipiv_create", "Can't pursue computation since no more tags are available for ipiv structure");
            return;
        }
        ipiv->mpitag_nextpiv = ipiv->mpitag_ipiv    + ipiv->mt;
        ipiv->mpitag_prevpiv = ipiv->mpitag_nextpiv + ipiv->mt;
        ipiv->mpitag_perm    = ipiv->mpitag_prevpiv + ipiv->mt;
        ipiv->mpitag_invp    = ipiv->mpitag_perm    + ipiv->mt;
    }
#endif
}

/**
 *  Destroy ws_pivot runtime structures
 */
void RUNTIME_ipiv_destroy( CHAM_ipiv_t *ipiv )
{
    int                   i;
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->ipiv);

    for(i=0; i<(5 * ipiv->mt); i++) {
        if ( *handle != NULL ) {
            starpu_data_unregister( *handle );
            *handle = NULL;
        }
        handle++;
    }

    free( ipiv->ipiv    );
    ipiv->ipiv    = NULL;
    ipiv->nextpiv = NULL;
    ipiv->prevpiv = NULL;
    ipiv->perm    = NULL;
    ipiv->invp    = NULL;
    chameleon_starpu_tag_release( ipiv->mpitag_ipiv );
}

void *RUNTIME_ipiv_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->ipiv);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    handle += mm;
    assert( handle );

    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    const CHAM_desc_t *A = ipiv->desc;
    int owner = A->get_rankof( A, m, m );
    int ncols = (mm == (ipiv->mt-1)) ? ipiv->m - mm * ipiv->mb : ipiv->mb;

    starpu_vector_data_register( handle, -1, (uintptr_t)NULL, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
    {
        int64_t tag = ipiv->mpitag_ipiv + mm;
        starpu_mpi_data_register( *handle, tag, owner );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

    assert( *handle );
    return (void*)(*handle);
}

void *RUNTIME_nextpiv_getaddr( const CHAM_ipiv_t *ipiv, int m, int h )
{
    starpu_data_handle_t *nextpiv = (starpu_data_handle_t*)(ipiv->nextpiv);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    nextpiv += mm;
    assert( nextpiv );

    if ( *nextpiv != NULL ) {
        return (void*)(*nextpiv);
    }

    const CHAM_desc_t *A = ipiv->desc;
    int     owner = A->get_rankof( A, m, m );
    int     ncols = (mm == (A->nt-1)) ? A->n - mm * A->nb : A->nb;
    int64_t tag   = ipiv->mpitag_nextpiv + mm;

    cppi_register( nextpiv, A->dtyp, ncols, tag, owner );

    assert( *nextpiv );
    return (void*)(*nextpiv);
}

void *RUNTIME_prevpiv_getaddr( const CHAM_ipiv_t *ipiv, int m, int h )
{
    starpu_data_handle_t *prevpiv = (starpu_data_handle_t*)(ipiv->prevpiv);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    prevpiv += mm;
    assert( prevpiv );

    if ( *prevpiv != NULL ) {
        return (void*)(*prevpiv);
    }

    const CHAM_desc_t *A = ipiv->desc;
    int     owner = A->get_rankof( A, m, m );
    int     ncols = (mm == (A->nt-1)) ? A->n - mm * A->nb : A->nb;
    int64_t tag   = ipiv->mpitag_prevpiv + mm;

    cppi_register( prevpiv, A->dtyp, ncols, tag, owner );

    assert( *prevpiv );
    return (void*)(*prevpiv);
}

void *RUNTIME_perm_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->perm);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    handle += mm;
    assert( handle );

    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    const CHAM_desc_t *A = ipiv->desc;
    int owner = A->get_rankof( A, m, m );
    int ncols = ipiv->mb;

    starpu_vector_data_register( handle, -1, (uintptr_t)NULL, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
    {
        int64_t tag = ipiv->mpitag_perm + mm;
        starpu_mpi_data_register( *handle, tag, owner );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

    assert( *handle );
    return (void*)(*handle);
}

void *RUNTIME_invp_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->invp);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    handle += mm;
    assert( handle );

    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    const CHAM_desc_t *A = ipiv->desc;
    int owner = A->get_rankof( A, m, m );
    int ncols = ipiv->mb;

    starpu_vector_data_register( handle, -1, (uintptr_t)NULL, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
    {
        int64_t tag = ipiv->mpitag_invp + mm;
        starpu_mpi_data_register( *handle, tag, owner );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

    assert( *handle );
    return (void*)(*handle);
}

void RUNTIME_ipiv_flushk( const RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle;
    const CHAM_desc_t *A = ipiv->desc;
    int64_t mm = m + ( ipiv->i / ipiv->mb );

    handle = (starpu_data_handle_t*)(ipiv->nextpiv);
    handle += mm;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == A->myrank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    handle = (starpu_data_handle_t*)(ipiv->prevpiv);
    handle += mm;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == A->myrank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    (void)sequence;
    (void)ipiv;
    (void)m;
}

void RUNTIME_ipiv_flush( const RUNTIME_sequence_t *sequence,
                         const CHAM_ipiv_t        *ipiv )
{
    int m;

    for (m = 0; m < ipiv->mt; m++)
    {
        RUNTIME_ipiv_flushk( sequence, ipiv, m );
    }
}

void RUNTIME_perm_flushk( const RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle;
    const CHAM_desc_t *A = ipiv->desc;
    int64_t mm = m + ( ipiv->i / ipiv->mb );

    handle = (starpu_data_handle_t*)(ipiv->perm);
    handle += mm;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == A->myrank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    handle = (starpu_data_handle_t*)(ipiv->invp);
    handle += mm;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == A->myrank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    (void)sequence;
    (void)ipiv;
    (void)m;
}

void RUNTIME_ipiv_gather( const RUNTIME_sequence_t *sequence,
                          CHAM_ipiv_t *desc, int *ipiv, int node )
{
    int64_t mt   = desc->mt;
    int64_t mb   = desc->mb;
    int64_t tag  = chameleon_starpu_tag_book( (int64_t)(desc->mt) );
    int     rank = CHAMELEON_Comm_rank();
    int     owner = rank;
    int     m;

    for (m = 0; m < mt; m++, ipiv += mb) {
        starpu_data_handle_t ipiv_src = RUNTIME_ipiv_getaddr( desc, m );

#if defined(CHAMELEON_USE_MPI)
        owner = starpu_mpi_data_get_rank( ipiv_src );
        if ( node != owner ) {
            starpu_mpi_tag_t tag = starpu_mpi_data_get_tag( ipiv_src );

            if ( rank == node )
            {
                /* Need to receive the data */
                int already_received = starpu_mpi_cached_receive_set( ipiv_src );
                if (already_received == 0)
                {
                    MPI_Status status;
                    starpu_mpi_recv( ipiv_src, owner, tag, sequence->comm, &status );
                }
            }
            else if ( rank == owner )
            {
                /* Need to send the data */
                int already_sent = starpu_mpi_cached_send_set( ipiv_src, node );
                if (already_sent == 0)
                {
                    starpu_mpi_send( ipiv_src, node, tag, sequence->comm );
                }
            }
        }
        if ( rank == node )
#endif
        {
            starpu_data_handle_t ipiv_dst;
            int       ncols     = (m == (mt-1)) ? desc->m - m * mb : mb;
            uintptr_t ipivptr   = (rank == node) ? (uintptr_t)ipiv : 0;
            int       home_node = (rank == node) ? STARPU_MAIN_RAM : -1;

            starpu_vector_data_register( &ipiv_dst, home_node, ipivptr, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
            starpu_mpi_data_register( ipiv_dst, tag + m, 0 );
#endif /* defined(CHAMELEON_USE_MPI) */

            assert( ipiv_dst );

            starpu_data_cpy( ipiv_dst, ipiv_src, 0, NULL, NULL );
            starpu_data_unregister( ipiv_dst );
        }
    }
}
