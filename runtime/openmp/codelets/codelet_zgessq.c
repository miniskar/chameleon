/**
 *
 * @file openmp/codelet_zgessq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessq StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2018-06-20
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void INSERT_TASK_zgessq( const RUNTIME_option_t *options,
                        int m, int n,
                        const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    double *ptrScaleSum = RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn);
#pragma omp task firstprivate(m, n, ptrA, lda, ptrScaleSum) depend(in:ptrA[0]) depend(inout:ptrScaleSum[0])
    CORE_zgessq( m, n, ptrA, lda, &ptrScaleSum[0], &ptrScaleSum[1] );
}
