/**
 *
 * @file openmp/codelet_zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2018-06-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void INSERT_TASK_zlange(const RUNTIME_option_t *options,
                       cham_normtype_t norm, int M, int N, int NB,
                       const CHAM_desc_t *A, int Am, int An, int LDA,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    double *ptrB = RTBLKADDR(B, double, Bm, Bn);
    int ws_size = options->ws_wsize;
#pragma omp task firstprivate(ws_size, M, N, ptrA, LDA, ptrB, options) depend(in:ptrA[0]) depend(inout:ptrB[0])
    {
      double work[ws_size];
      CORE_zlange( norm, M, N, ptrA, LDA, work, ptrB);
    }
}

void INSERT_TASK_zlange_max(const RUNTIME_option_t *options,
                           const CHAM_desc_t *A, int Am, int An,
                           const CHAM_desc_t *B, int Bm, int Bn)
{
    double *ptrA = RTBLKADDR(A, double, Am, An);
    double *ptrB = RTBLKADDR(B, double, Bm, Bn);

#pragma omp task firstprivate(ptrA, ptrB) depend(in:ptrA[0]) depend(inout:ptrB[0])
    {
        if ( *ptrA > *ptrB )
            *ptrB = *ptrA;
    }
}
