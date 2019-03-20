/**
 *
 * @file openmp/codelet_ztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpqrt StarPU codelet
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @date 2018-06-15
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    CHAMELEON_Complex64_t *ptrB = RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn);
    CHAMELEON_Complex64_t *ptrT = RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn);
    int ws_size = options->ws_wsize;

#pragma omp task firstprivate(ws_size, M, N, L, ib, ptrT, ldt, ptrA, lda, ptrB, ldb) depend(inout:ptrA[0], ptrB[0]) depend(out:ptrT[0])
    {
      CHAMELEON_Complex64_t tmp[ws_size];

      CORE_zlaset( ChamUpperLower, ib, N, 0., 0., ptrT, ldt );
      CORE_ztpqrt( M, N, L, ib,
                   ptrA, lda, ptrB, ldb, ptrT, ldt, tmp );
    }
}
