/**
 *
 * @file quark/codelet_zgetrf_nopiv_percol.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv_percol Quark codelets
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

void INSERT_TASK_zgetrf_nopiv_percol_diag( const RUNTIME_option_t *options,
                                           int m, int n, int k,
                                           const CHAM_desc_t *A, int Am, int An,
                                           const CHAM_desc_t *U, int Um, int Un,
                                           int iinfo )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)k;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)iinfo;
}

void INSERT_TASK_zgetrf_nopiv_percol_trsm( const RUNTIME_option_t *options,
                                           int m, int n, int k,
                                           const CHAM_desc_t *A, int Am, int An,
                                           const CHAM_desc_t *U, int Um, int Un )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)k;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
}
