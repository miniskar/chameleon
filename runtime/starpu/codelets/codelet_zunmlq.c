/**
 *
 * @file starpu/codelet_zunmlq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zunmlq_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    const CHAMELEON_Complex64_t *A;
    int lda;
    const CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *C;
    int ldc;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    A    = (const CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    T    = (const CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    WORK = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args(cl_arg, &side, &trans, &m, &n, &k, &ib,
                               &lda, &ldt, &ldc, &ldwork);

    CORE_zunmlq(side, trans, m, n, k, ib,
                A, lda, T, ldt, C, ldc, WORK, ldwork);
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zunmlq_cuda_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    const cuDoubleComplex *A, *T;
    cuDoubleComplex *C, *WORK;
    int lda, ldt, ldc, ldwork;

    starpu_codelet_unpack_args(cl_arg, &side, &trans, &m, &n, &k, &ib,
                               &lda, &ldt, &ldc, &ldwork);

    A    = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    T    = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C    = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    WORK = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[3]); /* ib * nb */

    RUNTIME_getStream(stream);

    CUDA_zunmlqt(
            side, trans, m, n, k, ib,
            A, lda, T, ldt, C, ldc, WORK, ldwork, stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zunmlq, 4, cl_zunmlq_cpu_func, cl_zunmlq_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 *  CORE_zunmlq overwrites the general complex M-by-N tile C with
 *
 *                    SIDE = 'L'     SIDE = 'R'
 *    TRANS = 'N':      Q * C          C * Q
 *    TRANS = 'C':      Q^H * C       C * Q^H
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(k) . . . H(2) H(1)
 *
 *  as returned by CORE_zgelqt. Q is of order M if SIDE = 'L' and of order N
 *  if SIDE = 'R'.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg ChamLeft  : apply Q or Q^H from the Left;
 *         @arg ChamRight : apply Q or Q^H from the Right.
 *
 * @param[in] trans
 *         @arg ChamNoTrans   :  No transpose, apply Q;
 *         @arg ChamConjTrans :  Transpose, apply Q^H.
 *
 * @param[in] M
 *         The number of rows of the tile C.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile C.  N >= 0.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *         If SIDE = ChamLeft,  M >= K >= 0;
 *         if SIDE = ChamRight, N >= K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] A
 *         Dimension:  (LDA,M) if SIDE = ChamLeft,
 *                     (LDA,N) if SIDE = ChamRight,
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_zgelqt in the first k rows of its array argument A.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the M-by-N tile C.
 *         On exit, C is overwritten by Q*C or Q^T*C or C*Q^T or C*Q.
 *
 * @param[in] LDC
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *         If SIDE = ChamLeft,  LDWORK >= max(1,N);
 *         if SIDE = ChamRight, LDWORK >= max(1,M).
 *
 *******************************************************************************
 *
 *          @retval CHAMELEON_SUCCESS successful exit
 *          @retval <0 if -i, the i-th argument had an illegal value
 *
 */
void INSERT_TASK_zunmlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                         const CHAM_desc_t *C, int Cm, int Cn, int ldc )
{
    struct starpu_codelet *codelet = &cl_zunmlq;
    void (*callback)(void*) = options->profiling ? cl_zunmlq_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &side,              sizeof(int),
        STARPU_VALUE,    &trans,             sizeof(int),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_R,         RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_R,         RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        STARPU_VALUE,    &ldt,               sizeof(int),
        STARPU_RW,        RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_VALUE,    &ldc,               sizeof(int),
        /* ib * nb */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &nb,                sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zunmlq",
#endif
        0);
}
