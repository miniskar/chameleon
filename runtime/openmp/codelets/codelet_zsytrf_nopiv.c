/**
 *
 * @file openmp/codelet_zsytrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf_nopiv StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Marc Sergent
 * @date 2011-10-09
 * @precisions normal z -> c
 *
 */

#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
void INSERT_TASK_zsytrf_nopiv(const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An, int lda,
                             int iinfo)
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
#pragma omp task firstprivate(uplo, n, ptrA, lda) depend(inout:ptrA[0:Am*An])
    CORE_zsytf2_nopiv(uplo, n, ptrA, lda);
}
