/**
 *
 * @file openmp/codelet_ztrtri.c
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri OpenMP codelet
 *
 * @version 0.9.2
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2019-11-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_diag_t diag,
                       int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       int iinfo )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
#pragma omp task firstprivate( uplo, diag, n, tileA, iinfo ) depend( inout:tileA[0] )
    TCORE_ztrtri( uplo, diag, n, tileA, &iinfo );
}
