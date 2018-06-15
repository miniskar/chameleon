/**
 *
 * @file codelet_zlanhe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlanhe StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2018-06-20
 * @precisions normal z -> c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void INSERT_TASK_zlanhe(const RUNTIME_option_t *options,
                       cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                       const CHAM_desc_t *A, int Am, int An, int LDA,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    double *work = options->ws_worker;
    double *normA = RTBLKADDR(B, double, Bm, Bn);
#pragma omp task firstprivate(norm, uplo, N, ptrA, LDA, work, normA) depend(in:ptrA[0:Am*An]) depend(inout:normA[0:Bm*Bn])
    CORE_zlanhe( norm, uplo, N, ptrA, LDA, work, normA);
}
