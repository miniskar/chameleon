/**
 *
 * @file parsec/codelet_map.c
 *
 * @copyright 2018-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map PaRSEC codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2024-03-11
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks.h"

struct parsec_map_args_s {
    cham_uplo_t          uplo;
    int                  m, n;
    cham_map_operator_t *op_fcts;
    void                *op_args;
    const CHAM_desc_t   *desc[1];
};

static inline int
CORE_map_one_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    struct parsec_map_args_s *pargs = NULL;
    CHAM_tile_t              *tileA;

    parsec_dtd_unpack_args( this_task, &pargs, &tileA );
    pargs->op_fcts->cpufunc( pargs->op_args, pargs->uplo, pargs->m, pargs->n, 1,
                             pargs->desc[0], tileA );

    free( pargs );
}

static inline int
CORE_map_two_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    struct parsec_map_args_s *pargs = NULL;
    CHAM_tile_t              *tileA;
    CHAM_tile_t              *tileB;

    parsec_dtd_unpack_args( this_task, &pargs, &tileA, &tileB );
    pargs->op_fcts->cpufunc( pargs->op_args, pargs->uplo, pargs->m, pargs->n, 2,
                             pargs->desc[0], tileA, pargs->desc[1], tileB );

    free( pargs );
}

static inline int
CORE_map_three_parsec( parsec_execution_stream_t *context,
                       parsec_task_t             *this_task )
{
    struct parsec_map_args_s *pargs = NULL;
    CHAM_tile_t              *tileA;
    CHAM_tile_t              *tileB;
    CHAM_tile_t              *tileC;

    parsec_dtd_unpack_args( this_task, &pargs, &tileA, &tileB, &tileC );
    pargs->op_fcts->cpufunc( pargs->op_args, pargs->uplo, pargs->m, pargs->n, 3,
                             pargs->desc[0], tileA, pargs->desc[1], tileB,
                             pargs->desc[2], tileC );

    free( pargs );
}

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args )
{
    parsec_taskpool_t        *PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    struct parsec_map_args_s *pargs      = NULL;
    size_t                    pargs_size = 0;
    int                       i;

    if ( ( ndata < 0 ) || ( ndata > 3 ) ) {
        fprintf( stderr, "INSERT_TASK_map() can handle only 1 to 3 parameters\n" );
        return;
    }

    pargs_size = sizeof( struct parsec_map_args_s ) + (ndata - 1) * sizeof( CHAM_desc_t * );
    pargs = malloc( pargs_size );
    pargs->uplo    = uplo;
    pargs->m       = m;
    pargs->n       = n;
    pargs->op_fcts = op_fcts;
    pargs->op_args = op_args;
    for( i=0; i<ndata; i++ ) {
        pargs->desc[i] = data[i].desc;
    }

    switch( ndata ) {
    case 1:
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_map_one_parsec, options->priority, op_fcts->name,
            sizeof(struct parsec_map_args_s*), &pargs, VALUE,
            PASSED_BY_REF, RTBLKADDR( data[0].desc, void, m, n ), chameleon_parsec_get_arena_index( data[0].desc ) | cham_to_parsec_access( data[0].access ),
            PARSEC_DTD_ARG_END );
        break;

    case 2:
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_map_two_parsec, options->priority, op_fcts->name,
            sizeof(struct parsec_map_args_s*), &pargs, VALUE,
            PASSED_BY_REF, RTBLKADDR( data[0].desc, void, m, n ), chameleon_parsec_get_arena_index( data[0].desc ) | cham_to_parsec_access( data[0].access ),
            PASSED_BY_REF, RTBLKADDR( data[1].desc, void, m, n ), chameleon_parsec_get_arena_index( data[1].desc ) | cham_to_parsec_access( data[1].access ),
            PARSEC_DTD_ARG_END );
        break;

    case 3:
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_map_three_parsec, options->priority, op_fcts->name,
            sizeof(struct parsec_map_args_s*), &pargs, VALUE,
            PASSED_BY_REF, RTBLKADDR( data[0].desc, void, m, n ), chameleon_parsec_get_arena_index( data[0].desc ) | cham_to_parsec_access( data[0].access ),
            PASSED_BY_REF, RTBLKADDR( data[1].desc, void, m, n ), chameleon_parsec_get_arena_index( data[1].desc ) | cham_to_parsec_access( data[1].access ),
            PASSED_BY_REF, RTBLKADDR( data[2].desc, void, m, n ), chameleon_parsec_get_arena_index( data[2].desc ) | cham_to_parsec_access( data[2].access ),
            PARSEC_DTD_ARG_END );
        break;
    }
}
