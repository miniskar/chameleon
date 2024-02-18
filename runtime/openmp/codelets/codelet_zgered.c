/**
 *
 * @file openmp/codelet_zgered.c
 *
 * @copyright 2023-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgered OpenMP codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2023-07-06
 * @precisions normal z -> d
 *
 */
#include "chameleon_openmp.h"

void INSERT_TASK_zgered( const RUNTIME_option_t *options,
                         double threshold, double Anorm, int m, int n,
                         const CHAM_desc_t *A, int Am, int An )
{
    fprintf( stderr, "WARNING: gered kernel is not available with OpenMP\n" );

    (void)options;
    (void)threshold;
    (void)Anorm;
    (void)m;
    (void)n;
    (void)A;
    (void)Am;
    (void)An;
}
