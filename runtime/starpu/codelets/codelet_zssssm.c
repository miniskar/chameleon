/**
 *
 * @file starpu/codelet_zssssm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zssssm_cpu_func(void *descr[], void *cl_arg)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A1;
    int ldA1;
    CHAMELEON_Complex64_t *A2;
    int ldA2;
    CHAMELEON_Complex64_t *L1;
    int ldL1;
    CHAMELEON_Complex64_t *L2;
    int ldL2;
    int *IPIV;

    A1 = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2 = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    L1 = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    L2 = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);
    ldA1 = STARPU_MATRIX_GET_LD( descr[0] );
    ldA2 = STARPU_MATRIX_GET_LD( descr[1] );
    ldL1 = STARPU_MATRIX_GET_LD( descr[2] );
    ldL2 = STARPU_MATRIX_GET_LD( descr[3] );
    starpu_codelet_unpack_args(cl_arg, &m1, &n1, &m2, &n2, &k, &ib, &IPIV);
    CORE_zssssm(m1, n1, m2, n2, k, ib, A1, ldA1, A2, ldA2, L1, ldL1, L2, ldL2, IPIV);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zssssm, 4, cl_zssssm_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
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
 * @param[in] ldA1
 *         The leading dimension of the array A1.  ldA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is updated by the application of L (L1 L2).
 *
 * @param[in] ldA2
 *         The leading dimension of the array A2.  ldA2 >= max(1,M2).
 *
 * @param[in] L1
 *         The IB-by-K lower triangular tile as returned by
 *         CORE_ztstrf.
 *
 * @param[in] ldL1
 *         The leading dimension of the array L1.  ldL1 >= max(1,IB).
 *
 * @param[in] L2
 *         The M2-by-K tile as returned by CORE_ztstrf.
 *
 * @param[in] ldL2
 *         The leading dimension of the array L2.  ldL2 >= max(1,M2).
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_ztstrf.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */
void INSERT_TASK_zssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n, int ldA1,
                         const CHAM_desc_t *A2, int A2m, int A2n, int ldA2,
                         const CHAM_desc_t *L1, int L1m, int L1n, int ldL1,
                         const CHAM_desc_t *L2, int L2m, int L2n, int ldL2,
                         const int *IPIV )
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zssssm;
    void (*callback)(void*) = options->profiling ? cl_zssssm_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A1, A1m, A1n);
    CHAMELEON_ACCESS_RW(A2, A2m, A2n);
    CHAMELEON_ACCESS_R(L1, L1m, L1n);
    CHAMELEON_ACCESS_R(L2, L2m, L2n);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m1,                        sizeof(int),
        STARPU_VALUE,    &n1,                        sizeof(int),
        STARPU_VALUE,    &m2,                        sizeof(int),
        STARPU_VALUE,    &n2,                        sizeof(int),
        STARPU_VALUE,     &k,                        sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_RW,            RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n),
        STARPU_RW,            RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n),
        STARPU_R,             RTBLKADDR(L1, CHAMELEON_Complex64_t, L1m, L1n),
        STARPU_R,             RTBLKADDR(L2, CHAMELEON_Complex64_t, L2m, L2n),
        STARPU_VALUE,          &IPIV,                      sizeof(int*),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zssssm",
#endif
        0);
}
