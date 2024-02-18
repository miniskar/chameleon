/**
 *
 * @file openmp/codelet_zhemm.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm OpenMP codelet
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

void
INSERT_TASK_zhemm( const RUNTIME_option_t *options,
                   cham_side_t side, cham_uplo_t uplo,
                   int m, int n, int nb,
                   CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                const CHAM_desc_t *B, int Bm, int Bn,
                   CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );
#pragma omp task firstprivate( side, uplo, m, n, alpha, tileA, tileB, beta, tileC ) depend( in:tileA[0], tileB[0] ) depend( inout:tileC[0] )
    TCORE_zhemm( side, uplo,
                 m, n,
                 alpha, tileA,
                 tileB,
                 beta, tileC );

    (void)options;
    (void)nb;
}

void
INSERT_TASK_zhemm_Astat( const RUNTIME_option_t *options,
                         cham_side_t side, cham_uplo_t uplo,
                         int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                      const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    INSERT_TASK_zhemm( options, side, uplo, m, n, nb,
                       alpha, A, Am, An, B, Bm, Bn,
                       beta, C, Cm, Cn );
}
