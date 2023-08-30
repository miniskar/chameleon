/**
 *
 * @file starpu/codelet_zgetrf_percol.c
 *
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel StarPU codelets
 *
 * @version 1.3.0
 * @comment Codelets to perform panel factorization with partial pivoting
 *
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2023-08-22
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

CHAMELEON_CL_CB( zgetrf_percol_diag,    cti_handle_get_m(task->handles[0]), 0, 0, M );
CHAMELEON_CL_CB( zgetrf_percol_offdiag, cti_handle_get_m(task->handles[0]), 0, 0, M );

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_percol_diag_cpu_func(void *descr[], void *cl_arg)
{
    int                 h, m0;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t  *request;
    CHAM_tile_t        *tileA;
    int                *ipiv;
    cppi_interface_t   *nextpiv;
    cppi_interface_t   *prevpiv;

    starpu_codelet_unpack_args( cl_arg, &h, &m0,
                                &sequence, &request );

    tileA   = cti_interface_get(descr[0]);
    ipiv    = (int *)STARPU_VECTOR_GET_PTR(descr[1]);
    nextpiv = (cppi_interface_t*) descr[2];
    prevpiv = (cppi_interface_t*) descr[3];

    if ( h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv before call: " );
    }
    if ( h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv before call: " );
    }

    /*
     * Make sure the nextpiv interface store the right information about the
     * column and diagonal row for the reduction
     */
    nextpiv->h        = h;
    nextpiv->has_diag = 1;

    CORE_zgetrf_panel_diag( tileA->m, tileA->n, h, m0,
                            CHAM_tile_get_ptr( tileA ), tileA->ld,
                            ipiv, &(nextpiv->pivot), &(prevpiv->pivot) );

    if ( h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv after call: " );
    }
    if ( h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv after call: " );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_percol_diag, cl_zgetrf_percol_diag_cpu_func );

void INSERT_TASK_zgetrf_percol_diag( const RUNTIME_option_t *options,
                                     int h, int m0,
                                     CHAM_desc_t *A, int Am, int An,
                                     CHAM_ipiv_t *ipiv )
{
    struct starpu_codelet *codelet = &cl_zgetrf_percol_diag;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_percol_diag_callback : NULL;

    int access_ipiv = ( h == 0 )       ? STARPU_W    : STARPU_RW;
    int access_npiv = ( h == ipiv->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )       ? STARPU_NONE : STARPU_R;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &h,                   sizeof(int),
        STARPU_VALUE,             &m0,                  sizeof(int),
        STARPU_VALUE,             &(options->sequence), sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,             &(options->request),  sizeof(RUNTIME_request_t*),
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        access_ipiv,              RUNTIME_ipiv_getaddr( ipiv, An ),
        access_npiv,              RUNTIME_pivot_getaddr( ipiv, An, h   ),
        access_ppiv,              RUNTIME_pivot_getaddr( ipiv, An, h-1 ),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf_percol_diag",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_percol_offdiag_cpu_func(void *descr[], void *cl_arg)
{
    int                 h, m0;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t  *request;
    CHAM_tile_t        *tileA;
    cppi_interface_t   *nextpiv;
    cppi_interface_t   *prevpiv;

    starpu_codelet_unpack_args( cl_arg, &h, &m0, &sequence, &request );

    tileA   = cti_interface_get(descr[0]);
    nextpiv = (cppi_interface_t*) descr[1];
    prevpiv = (cppi_interface_t*) descr[2];

    nextpiv->h = h; /* Initialize in case it uses a copy */

    CORE_zgetrf_panel_offdiag( tileA->m, tileA->n, h, m0,
                               CHAM_tile_get_ptr(tileA), tileA->ld,
                               &(nextpiv->pivot), &(prevpiv->pivot) );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_percol_offdiag, cl_zgetrf_percol_offdiag_cpu_func)

void INSERT_TASK_zgetrf_percol_offdiag( const RUNTIME_option_t *options,
                                        int h, int m0,
                                        CHAM_desc_t *A, int Am, int An,
                                        CHAM_ipiv_t *ipiv )
{
    struct starpu_codelet *codelet = &cl_zgetrf_percol_offdiag;

    void (*callback)(void*) = options->profiling ? cl_zgetrf_percol_offdiag_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &h,                   sizeof(int),
        STARPU_VALUE,    &m0,                  sizeof(int),
        STARPU_VALUE,    &(options->sequence), sizeof(RUNTIME_sequence_t *),
        STARPU_VALUE,    &(options->request),  sizeof(RUNTIME_request_t *),
        STARPU_RW,       RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_REDUX,    RUNTIME_pivot_getaddr( ipiv, An, h   ),
        STARPU_R,        RUNTIME_pivot_getaddr( ipiv, An, h-1 ),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf_percol_offdiag",
#endif
        0);
}
