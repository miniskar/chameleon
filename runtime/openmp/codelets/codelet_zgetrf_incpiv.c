/**
 *
 * @file openmp/codelet_zgetrf_incpiv.c
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv OpenMP codelet
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

void INSERT_TASK_zgetrf_incpiv( const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *L, int Lm, int Ln,
                              int *IPIV,
                              cham_bool_t check_info, int iinfo )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    int info = 0;
#pragma omp task firstprivate( m, n, ib, tileA, IPIV, info ) depend( out:IPIV[0] ) depend( inout:tileA[0] )
    TCORE_zgetrf_incpiv( m, n, ib, tileA, IPIV, &info );
}
