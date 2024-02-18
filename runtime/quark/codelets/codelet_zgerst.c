/**
 *
 * @file quark/codelet_zgerst.c
 *
 * @copyright 2023-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgerst Quark codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2023-07-06
 * @precisions normal z -> d
 *
 */
#include "chameleon_quark.h"

void INSERT_TASK_zgerst( const RUNTIME_option_t *options,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An )
{
    fprintf( stderr, "WARNING: gerst kernel is not available with Quark\n" );

    (void)options;
    (void)m;
    (void)n;
    (void)A;
    (void)Am;
    (void)An;
}
