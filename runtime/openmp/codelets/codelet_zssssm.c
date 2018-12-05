/**
 *
 * @file openmp/codelet_zssssm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */

#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zssssm applies the LU factorization update from a complex
 *  matrix formed by a lower triangular IB-by-K tile L1 on top of a
 *  M2-by-K tile L2 to a second complex matrix formed by a M1-by-N1
 *  tile A1 on top of a M2-by-N2 tile A2 (N1 == N2).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M1
 *         The number of rows of the tile A1.  M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1.  N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2 and of the tile L2.
 *         M2 >= 0.
 *
 * @param[in] N2
 *         The number of columns of the tile A2.  N2 >= 0.
 *
 * @param[in] K
 *         The number of columns of the tiles L1 and L2.  K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is updated by the application of L (L1 L2).
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is updated by the application of L (L1 L2).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M2).
 *
 * @param[in] L1
 *         The IB-by-K lower triangular tile as returned by
 *         CORE_ztstrf.
 *
 * @param[in] LDL1
 *         The leading dimension of the array L1.  LDL1 >= max(1,IB).
 *
 * @param[in] L2
 *         The M2-by-K tile as returned by CORE_ztstrf.
 *
 * @param[in] LDL2
 *         The leading dimension of the array L2.  LDL2 >= max(1,M2).
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_ztstrf.
 *
 *******************************************************************************
 *
 * @return
 *         \retval CHAMELEON_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */

void INSERT_TASK_zssssm(const RUNTIME_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                       const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                       const CHAM_desc_t *L1, int L1m, int L1n, int ldl1,
                       const CHAM_desc_t *L2, int L2m, int L2n, int ldl2,
                       const int *IPIV)
{
    CHAMELEON_Complex64_t *ptrA1 = RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n);
    CHAMELEON_Complex64_t *ptrA2 = RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n);
    CHAMELEON_Complex64_t *ptrL1 = RTBLKADDR(L1, CHAMELEON_Complex64_t, L1m, L1n);
    CHAMELEON_Complex64_t *ptrL2 = RTBLKADDR(L2, CHAMELEON_Complex64_t, L2m, L2n);
#pragma omp task firstprivate(m1, n1, m2, n2, k, ib, ptrA1, ptrA2, ptrL1, ptrL2, lda1, lda2, ldl1, ldl2, IPIV)\
    depend(inout:ptrA1[0:A1m*A1n])\
    depend(inout:ptrA2[0:A2m*A2n])\
    depend(in:ptrL1[0:L1m*L1n])\
    depend(in:ptrL2[0:L2m*L2n])
    CORE_zssssm(m1, n1, m2, n2, k, ib, ptrA1, lda1, ptrA2, lda2, ptrL1, ldl1, ptrL2, ldl2, IPIV);
}
