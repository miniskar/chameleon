/**
 *
 * @file openmp/codelet_zherk.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk OpenMP codelet
 *
 * @version 1.2.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2022-02-22
 * @precisions normal z -> c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, cham_trans_t trans,
                      int n, int k, int nb,
                      double alpha, const CHAM_desc_t *A, int Am, int An,
                      double beta, const CHAM_desc_t *C, int Cm, int Cn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );
#pragma omp task firstprivate( uplo, trans, n, k, alpha, tileA, beta, tileC ) depend( in:tileA[0] ) depend( inout:tileC[0] )
    TCORE_zherk( uplo, trans,
        n, k,
        alpha, tileA,
        beta, tileC );

    (void)options;
    (void)nb;
}
