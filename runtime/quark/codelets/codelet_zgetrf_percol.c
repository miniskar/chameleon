/**
 *
 * @file quark/codelet_zgetrf_percol.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_percol Quark codelets
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
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zgetrf_percol_diag( const RUNTIME_option_t *options,
                                     int h, int m0,
                                     CHAM_desc_t *A, int Am, int An,
                                     CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)h;
    (void)m0;
    (void)A;
    (void)Am;
    (void)An;
    (void)ipiv;
}

void INSERT_TASK_zgetrf_percol_offdiag( const RUNTIME_option_t *options,
                                        int h, int m0,
                                        CHAM_desc_t *A, int Am, int An,
                                        CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)h;
    (void)m0;
    (void)A;
    (void)Am;
    (void)An;
    (void)ipiv;
}
