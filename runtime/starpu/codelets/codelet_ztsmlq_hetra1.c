/**
 *
 * @file starpu/codelet_ztsmlq_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmlq_hetra1 StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 */
void INSERT_TASK_ztsmlq_hetra1(const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                              const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                              const CHAM_desc_t *V,  int Vm,  int Vn,  int ldv,
                              const CHAM_desc_t *T,  int Tm,  int Tn,  int ldt)
{
    struct starpu_codelet *codelet = &cl_ztsmlq_hetra1;
    void (*callback)(void*) = options->profiling ? cl_ztsmlq_hetra1_callback : NULL;

    int ldwork = side == ChamLeft ? ib : nb;

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
        STARPU_VALUE,    &lda1,              sizeof(int),
        STARPU_RW,        RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n),
        STARPU_VALUE,    &lda2,              sizeof(int),
        STARPU_R,         RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn),
        STARPU_VALUE,    &ldv,               sizeof(int),
        STARPU_R,         RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        STARPU_VALUE,    &ldt,               sizeof(int),
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &ldwork,            sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "ztsmlq_hetra1",
#endif
        0);
}

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
    int lda1;
    CHAMELEON_Complex64_t *A2;
    int lda2;
    CHAMELEON_Complex64_t *V;
    int ldv;
    CHAMELEON_Complex64_t *T;
    int ldt;

    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    A1    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    V     = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    T     = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]);
    WORK  = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[4]); /* ib * nb */

    starpu_codelet_unpack_args(cl_arg, &side, &trans, &m1, &n1, &m2, &n2, &k,
                               &ib, &nb, &lda1, &lda2, &ldv, &ldt, &ldwork);
    CORE_ztsmlq_hetra1(side, trans, m1, n1, m2, n2, k,
                       ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztsmlq_hetra1, 5, cl_ztsmlq_hetra1_cpu_func)
