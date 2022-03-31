/**
 *
 * @file pzgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf parallel algorithm
 *
 * @version 1.2.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Matthieu Kuhn
 * @date 2022-02-22
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n

/*
 * All the functions below are panel factorization variant.
 * The parameters are:
 *   @param[inout] ws
 *      The data structure associated to the algorithm that holds all extra
 *      information that may be needed for LU factorization
 *
 *   @param[inout] A
 *      The descriptor of the full matrix A (not just the panel)
 *
 *   @param[in] k
 *      The index of the column to factorize
 *
 *   @param[in] ib
 *      The index of the column to factorize
 *
 *   @param[inout] options
 *      The runtime options data structure to pass through all insert_task calls.
 */
static inline void
chameleon_pzgetrf_panel_facto_nopiv( void             *ws,
                                     CHAM_desc_t      *A,
                                     int               k,
                                     RUNTIME_option_t *options )
{
    const CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t) 1.0;
    int m, tempkm, tempkn, tempmm;

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

    /*
     * Algorithm per block without pivoting
     */
    INSERT_TASK_zgetrf_nopiv(
        options,
        tempkm, tempkn, 32, A->mb,
         A(k, k), 0);

    for (m = k+1; m < A->mt; m++) {
        tempmm = (m == (A->mt - 1)) ? A->m - m * A->mb : A->mb;
        INSERT_TASK_ztrsm(
            options,
            ChamRight, ChamUpper, ChamNoTrans, ChamNonUnit,
            tempmm, tempkn, A->mb,
            zone, A(k, k),
                  A(m, k) );
    }
}
static inline void
chameleon_pzgetrf_panel_facto( void             *ws,
                               CHAM_desc_t      *A,
                               int               k,
                               RUNTIME_option_t *options )
{
    chameleon_pzgetrf_panel_facto_nopiv( ws, A, k, options );
}

/**
 *  Permutation of the panel n at step k
 */
static inline void
chameleon_pzgetrf_panel_permute( void             *ws,
                                 CHAM_desc_t      *A,
                                 int               k,
                                 int               n,
                                 RUNTIME_option_t *options )
{
    (void)ws;
    (void)A;
    (void)k;
    (void)n;
    (void)options;
}

static inline void
chameleon_pzgetrf_panel_update( void             *ws,
                                CHAM_desc_t      *A,
                                int               k,
                                int               n,
                                RUNTIME_option_t *options )
{
    const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    int m, tempkm, tempmm, tempnn;

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

    chameleon_pzgetrf_panel_permute( ws, A, k, n, options );

    INSERT_TASK_ztrsm(
        options,
        ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
        tempkm, tempnn, A->mb,
        zone, A(k, k),
              A(k, n) );

    for (m = k+1; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

        INSERT_TASK_zgemm(
            options,
            ChamNoTrans, ChamNoTrans,
            tempmm, tempnn, A->mb, A->mb,
            mzone, A(m, k),
                   A(k, n),
            zone,  A(m, n) );
    }

    RUNTIME_data_flush( options->sequence, A(k, n) );
}

/**
 *  Parallel tile LU factorization with no pivoting - dynamic scheduling
 */
void chameleon_pzgetrf( void               *ws,
                        CHAM_desc_t        *A,
                        RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t  *request )
{
    CHAM_context_t  *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int min_mnt = chameleon_min( A->mt, A->nt );

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    for (k = 0; k < min_mnt; k++) {
        RUNTIME_iteration_push( chamctxt, k );

        options.priority = A->nt;
        chameleon_pzgetrf_panel_facto( ws, A, k, &options );

        for (n = k+1; n < A->nt; n++) {
            options.priority = A->nt-n;
            chameleon_pzgetrf_panel_update( ws, A, k, n, &options );
        }

        /* Flush panel k */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_flush( sequence, A(m, k) );
        }

        RUNTIME_iteration_pop( chamctxt );
    }

    /* Backward pivoting */
    for (k = 1; k < min_mnt; k++) {
        for (n = 0; n < k; n++) {
            chameleon_pzgetrf_panel_permute( ws, A, k, n, &options );
        }
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
