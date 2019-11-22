/**
 *
 * @file openmp/codelet_map.c
 *
 * @copyright 2018-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map OpenMP codelet
 *
 * @version 0.9.2
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2019-11-19
 *
 */
#include "chameleon_openmp.h"

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

#pragma omp task depend( inout: tileA[0] )
    {
        op_fct( A, uplo, Am, An, tileA, op_args );
    }

}
