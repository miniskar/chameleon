/**
 *
 * @file codelet_zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Philippe Virouleau
 * @date 2018-06-20
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void INSERT_TASK_zgetrf(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       int *IPIV,
                       cham_bool_t check_info, int iinfo)
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    int info = 0;
#pragma omp task firstprivate(m, n, ptrA, lda, IPIV, info) depend(inout:ptrA[0])
    CORE_zgetrf( m, n, ptrA, lda, IPIV, &info );
}
