/**
 *
 * @file starpu/codelet_ztsmlq_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmlq_hetra1 StarPU codelet
 *
 * @version 0.9.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @author Lucas Barros de Assis
 * @date 2016-12-09
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztsmlq_hetra1_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    int nb;
    CHAMELEON_Complex64_t *A1;
    int ldA1;
    CHAMELEON_Complex64_t *A2;
    int ldA2;
    CHAMELEON_Complex64_t *V;
    int ldV;
    CHAMELEON_Complex64_t *T;
    int ldT;

    CHAMELEON_Complex64_t *WORK;
    int ldWORK;

    A1    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    V     = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    T     = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);
    WORK  = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[4]); /* ib * nb */
    ldA1 = STARPU_MATRIX_GET_LD( descr[0] );
    ldA2 = STARPU_MATRIX_GET_LD( descr[1] );
    ldV = STARPU_MATRIX_GET_LD( descr[2] );
    ldT = STARPU_MATRIX_GET_LD( descr[3] );
    starpu_codelet_unpack_args( cl_arg, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib, &nb, &ldWORK);
    CORE_ztsmlq_hetra1(side, trans, m1, n1, m2, n2, k,
                       ib, A1, ldA1, A2, ldA2, V, ldV, T, ldT, WORK, ldWORK);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztsmlq_hetra1, 5, cl_ztsmlq_hetra1_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_ztsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n, int ldA1,
                                const CHAM_desc_t *A2, int A2m, int A2n, int ldA2,
                                const CHAM_desc_t *V,  int Vm,  int Vn,  int ldV,
                                const CHAM_desc_t *T,  int Tm,  int Tn,  int ldT )
{
    struct starpu_codelet *codelet = &cl_ztsmlq_hetra1;
    void (*callback)(void*) = options->profiling ? cl_ztsmlq_hetra1_callback : NULL;

    int ldWORK = side == ChamLeft ? ib : nb;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A1, A1m, A1n);
    CHAMELEON_ACCESS_RW(A2, A2m, A2n);
    CHAMELEON_ACCESS_R(V, Vm, Vn);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &side,              sizeof(int),
        STARPU_VALUE,    &trans,             sizeof(int),
        STARPU_VALUE,    &m1,                sizeof(int),
        STARPU_VALUE,    &n1,                sizeof(int),
        STARPU_VALUE,    &m2,                sizeof(int),
        STARPU_VALUE,    &n2,                sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_VALUE,    &nb,                sizeof(int),
        STARPU_RW,        RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n),
        STARPU_RW,        RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n),
        STARPU_R,         RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn),
        STARPU_R,         RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &ldWORK,            sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztsmlq_hetra1",
#endif
        0);
}
