/**
 *
 * @file openmp/codelet_zplgsy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplgsy StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Philippe Virouleau
 * @date 2018-06-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An, int lda,
                         int bigM, int m0, int n0, unsigned long long int seed )
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
#pragma omp task firstprivate(bump, m, n, ptrA, lda, bigM, m0, n0, seed) depend(inout:ptrA[0])
    CORE_zplgsy( bump, m, n, ptrA, lda, bigM, m0, n0, seed );
}
