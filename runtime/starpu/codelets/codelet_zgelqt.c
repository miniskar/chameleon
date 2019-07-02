/**
 *
 * @file starpu/codelet_zgelqt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqt StarPU codelet
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
static void cl_zgelqt_cpu_func(void *descr[], void *cl_arg)
{
    CHAMELEON_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    CHAMELEON_Complex64_t *A;
    int ldA;
    CHAMELEON_Complex64_t *T;
    int ldT;
    CHAMELEON_Complex64_t *TAU, *WORK;

    A   = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    T   = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    TAU = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]); /* max(m,n) + ib*n */
    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldT = STARPU_MATRIX_GET_LD( descr[1] );

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &h_work);

    WORK = TAU + chameleon_max( m, n );
    CORE_zlaset( ChamUpperLower, ib, m, 0., 0., T, ldT );
    CORE_zgelqt(m, n, ib, A, ldA, T, ldT, TAU, WORK);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgelqt, 3, cl_zgelqt_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 *  CORE_zgelqt - computes a LQ factorization of a complex M-by-N tile A: A = L * Q.
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
 *  A(i,i+1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-min(M,N) lower trapezoidal tile L (L is
 *         lower triangular if M <= N); the elements above the diagonal,
 *         with the array TAU, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] ldA
 *         The leading dimension of the array A.  ldA >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldT
 *         The leading dimension of the array T. ldT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 *          @retval CHAMELEON_SUCCESS successful exit
 *          @retval <0 if -i, the i-th argument had an illegal value
 *
 */
void INSERT_TASK_zgelqt(const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An, int ldA,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldT)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgelqt;
    void (*callback)(void*) = options->profiling ? cl_zgelqt_callback : NULL;
    CHAMELEON_starpu_ws_t *h_work = (CHAMELEON_starpu_ws_t*)(options->ws_host);

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_W(T, Tm, Tn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,         RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        /* max( nb * (ib+1), ib * (ib+nb) ) */
        STARPU_SCRATCH,   options->ws_worker,
        /* /\* ib*n + 3*ib*ib + max(m,n) *\/ */
        STARPU_VALUE,    &h_work,            sizeof(CHAMELEON_starpu_ws_t *),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgelqt",
#endif
        0);
    (void)ldT;
    (void)ldA;
}
