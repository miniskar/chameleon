/**
 *
 * @file openmp/codelet_zsyssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyssq StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, int n,
                        const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    double *ptrSCALESUMSQ = RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn);
#pragma omp task firstprivate(uplo, n, ptrA, lda, ptrSCALESUMSQ) depend(in:ptrA[0]) depend(inout:ptrSCALESUMSQ[0])
    CORE_zsyssq( uplo, n, ptrA, lda, &ptrSCALESUMSQ[0], &ptrSCALESUMSQ[1] );
}
