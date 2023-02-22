/**
 *
 * @file starpu/codelet_zpanel.c
 *
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel StarPU codelets
 *
 * @version 1.2.0
 * @comment Codelets to perform panel factorization with partial pivoting
 *
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2023-02-21
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

static const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
static const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_panel_nopiv_percol_diag_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t           *tileA, *tileU;
    int                    m, n, k, lda, iinfo;
    RUNTIME_sequence_t    *sequence;
    RUNTIME_request_t     *request;
    CHAMELEON_Complex64_t *A, *row, pivot;

    tileA = cti_interface_get( descr[0] );
    tileU = cti_interface_get( descr[1] );

    starpu_codelet_unpack_args( cl_arg, &m, &n, &k, &iinfo, &sequence, &request );

    A   = tileA->mat;
    lda = tileA->ld;
    row = tileU->mat;

    /* Shift to the diagonal element */
    A += k * (lda + 1);

    /* Extract row into buffer */
    cblas_zcopy( n-k, A, lda, row, 1 );

    /* Perform update on current diagonal block directly here */
    if ( *row == 0. ) {
        if ( sequence->status == CHAMELEON_SUCCESS ) {
            RUNTIME_sequence_flush( NULL, sequence, request, iinfo+k+1 );
        }
        return;
    }

    pivot = 1. / *row;
    cblas_zscal( m-k-1, CBLAS_SADDR( pivot ), A + 1, 1 );

    CORE_zgemm( ChamNoTrans, ChamNoTrans,
                m-k-1, n-k-1, 1,
                mzone, A   + 1,       lda,
                       row + 1,       1,
                zone,  A   + 1 + lda, lda );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_panel_nopiv_percol_diag, cl_zgetrf_panel_nopiv_percol_diag_cpu_func );

void INSERT_TASK_zgetrf_panel_nopiv_percol_diag( const RUNTIME_option_t *options,
                                                 int m, int n, int k,
                                                 const CHAM_desc_t *A, int Am, int An,
                                                 const CHAM_desc_t *U, int Um, int Un,
                                                 int iinfo )
{
    struct starpu_codelet *codelet = &cl_zgetrf_panel_nopiv_percol_diag;
    // void (*callback)(void*) = options->profiling ? cl_zgetrf_panel_nopiv_percol_diag_callback : NULL;
    void (*callback)(void*) = NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    CHAMELEON_ACCESS_W( U, Um, Un );
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m,                   sizeof(int),
        STARPU_VALUE,             &n,                   sizeof(int),
        STARPU_VALUE,             &k,                   sizeof(int),
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,                 RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_VALUE,             &iinfo,               sizeof(int),
        STARPU_VALUE,             &(options->sequence), sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,             &(options->request),  sizeof(RUNTIME_request_t*),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              "zgetrf_panel_nopiv_percol_diag",
#endif
        0);
}

/*
 * Update column blocs
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_panel_nopiv_percol_trsm_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t           *tileA, *tileU;
    int                    m, n, k, lda;
    CHAMELEON_Complex64_t *A, *row, pivot;

    tileA = cti_interface_get( descr[0] );
    tileU = cti_interface_get( descr[1] );

    starpu_codelet_unpack_args( cl_arg, &m, &n, &k );

    A   = tileA->mat;
    lda = tileA->ld;
    row = tileU->mat;

    /* Shift to the right column */
    A += k * lda;

    pivot = 1. / *row;
    cblas_zscal( m, CBLAS_SADDR( pivot ), A, 1 );

    /* Update trailing matrix from k+1 to n */
    CORE_zgemm( ChamNoTrans, ChamNoTrans,
                m, n-k-1, 1,
                mzone, A,       lda,
                       row + 1, 1,
                zone,  A + lda, lda );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_panel_nopiv_percol_trsm, cl_zgetrf_panel_nopiv_percol_trsm_cpu_func );

void INSERT_TASK_zgetrf_panel_nopiv_percol_trsm( const RUNTIME_option_t *options,
                                                 int m, int n, int k,
                                                 const CHAM_desc_t *A, int Am, int An,
                                                 const CHAM_desc_t *U, int Um, int Un )
{
    struct starpu_codelet *codelet = &cl_zgetrf_panel_nopiv_percol_trsm;
    // void (*callback)(void*) = options->profiling ? cl_zgetrf_panel_nopiv_percol_trsm_callback : NULL;
    void (*callback)(void*) = NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_R(U, Um, Un);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m, sizeof(int),
        STARPU_VALUE,             &n, sizeof(int),
        STARPU_VALUE,             &k, sizeof(int),
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,                 RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf_panel_nopiv_percol_trsm",
#endif
        0);
}

