/**
 *
 * @file quark/codelet_map.c
 *
 * @copyright 2018-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map Quark codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2018-09-24
 *
 */
#include "chameleon_openmp.h"

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t operator, void *op_args )
{
    char *ptrA = RTBLKADDR( A, char, Am, An );

#pragma omp task depend(inout: ptrA[0])
    {
        operator( A, uplo, Am, An, ptrA, op_args );
    }

}
