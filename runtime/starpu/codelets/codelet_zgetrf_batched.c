/**
 *
 * @file starpu/codelet_zgetrf_batched.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel batched StarPU codelets
 *
 * @version 1.2.0
 * @comment Codelets to perform batched panel factorization with partial pivoting
 *
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @date 2024-01-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

struct cl_getrf_batched_args_t {
    char                    *cl_name;
    int                      tasks_nbr;
    int                      diag;
    int                      h;
    int                      ib;
    int                      m[CHAMELEON_BATCH_SIZE];
    int                      n[CHAMELEON_BATCH_SIZE];
    int                      m0[CHAMELEON_BATCH_SIZE];
    struct starpu_data_descr handle_mode[CHAMELEON_BATCH_SIZE];
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgetrf_panel_offdiag_batched_cpu_func( void *descr[],
                                          void *cl_arg )
{
    struct cl_getrf_batched_args_t *clargs  = (struct cl_getrf_batched_args_t *) cl_arg;
    cppi_interface_t               *nextpiv = (cppi_interface_t*) descr[0];
    cppi_interface_t               *prevpiv = (cppi_interface_t*) descr[1];
    int                             i, m, n, h, m0, lda;
    CHAM_tile_t                    *tileA;

    nextpiv->h = clargs->h;

    for ( i = 0; i < clargs->tasks_nbr; i++ ) {
        tileA = cti_interface_get( descr[ i + 2 ] );
        lda   = tileA->ld;
        m     = clargs->m[ i ];
        n     = clargs->n[ i ];
        h     = clargs->h;
        m0    = clargs->m0[ i ];
        CORE_zgetrf_panel_offdiag( m, n, h, m0, n, CHAM_tile_get_ptr(tileA), lda,
                                   NULL, -1, &( nextpiv->pivot ), &( prevpiv->pivot ) );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

CODELETS_CPU( zgetrf_panel_offdiag_batched, cl_zgetrf_panel_offdiag_batched_cpu_func )

void
INSERT_TASK_zgetrf_panel_offdiag_batched( const RUNTIME_option_t *options,
                                          int m, int n, int h, int m0,
                                          void *ws,
                                          CHAM_desc_t *A, int Am, int An,
                                          void **clargs_ptr,
                                          CHAM_ipiv_t *ipiv )
{
    CHAM_tile_t *tileA      = A->get_blktile( A, Am, An );
    int          task_num   = 0;
    int          exec       = 0;
    int          batch_size = ((struct chameleon_pzgetrf_s *)ws)->batch_size;
    void (*callback)(void*) = NULL;
    struct cl_getrf_batched_args_t *clargs = *clargs_ptr;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_getrf_batched_args_t ) ) ;
        clargs->tasks_nbr   = 0;
        clargs->h           = h;
        clargs->cl_name     = "zgetrf_panel_offdiag_batched";

        *clargs_ptr = clargs;
    }

    task_num               = clargs->tasks_nbr;
    clargs->m[ task_num ]  = m;
    clargs->n[ task_num ]  = n;
    clargs->m0[ task_num ] = m0;
    clargs->handle_mode[ task_num ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    clargs->handle_mode[ task_num ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;
    /* Refine name */
    clargs->cl_name = chameleon_codelet_name( clargs->cl_name, 1,
                                              A->get_blktile( A, Am, An ) );

    if ( clargs->tasks_nbr == batch_size ) {
        rt_starpu_insert_task(
            &cl_zgetrf_panel_offdiag_batched,
            /* Task codelet arguments */
            STARPU_CL_ARGS,           clargs, sizeof(struct cl_getrf_batched_args_t),
            STARPU_REDUX,             RUNTIME_pivot_getaddr( ipiv, An, h   ),
            STARPU_R,                 RUNTIME_pivot_getaddr( ipiv, An, h-1 ),
            STARPU_DATA_MODE_ARRAY,   clargs->handle_mode, clargs->tasks_nbr,
            STARPU_PRIORITY,          options->priority,
            STARPU_CALLBACK,          callback,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME,              clargs->cl_name,
#endif
            0);

        /* clargs is freed by starpu. */
        *clargs_ptr = NULL;
    }
}

void
INSERT_TASK_zgetrf_panel_offdiag_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An,
                                                void **clargs_ptr,
                                                CHAM_ipiv_t *ipiv )
{
    void (*callback)(void*) = NULL;
    struct cl_getrf_batched_args_t *clargs = *clargs_ptr;

    if ( clargs == NULL ) {
        return;
    }

    rt_starpu_insert_task(
        &cl_zgetrf_panel_offdiag_batched,
        /* Task codelet arguments */
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_getrf_batched_args_t),
        STARPU_REDUX,             RUNTIME_pivot_getaddr( ipiv, An, clargs->h   ),
        STARPU_R,                 RUNTIME_pivot_getaddr( ipiv, An, clargs->h-1 ),
        STARPU_DATA_MODE_ARRAY,   clargs->handle_mode, clargs->tasks_nbr,
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              clargs->cl_name,
#endif
        0);

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}
