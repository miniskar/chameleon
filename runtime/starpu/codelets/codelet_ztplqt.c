/**
 *
 * @file starpu/codelet_ztplqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztplqt StarPU codelet
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @date 2018-01-31
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztplqt_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int L;
    int ib;
    CHAMELEON_Complex64_t *A;
    int ldA;
    CHAMELEON_Complex64_t *B;
    int ldB;
    CHAMELEON_Complex64_t *T;
    int ldT;
    CHAMELEON_Complex64_t *WORK;

    A    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    T    = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    WORK = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]); /* ib * nb */
    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldB = STARPU_MATRIX_GET_LD( descr[1] );
    ldT = STARPU_MATRIX_GET_LD( descr[2] );
    starpu_codelet_unpack_args( cl_arg, &M, &N, &L, &ib );

    CORE_zlaset( ChamUpperLower, ib, M, 0., 0., T, ldT );
    CORE_ztplqt( M, N, L, ib,
                 A, ldA, B, ldB, T, ldT, WORK );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztplqt, 4, cl_ztplqt_cpu_func)

void INSERT_TASK_ztplqt( const RUNTIME_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int ldA,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldB,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldT )
{
    struct starpu_codelet *codelet = &cl_ztplqt;
    void (*callback)(void*) = options->profiling ? cl_ztplqt_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_ACCESS_W(T, Tm, Tn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,    sizeof(int),
        STARPU_RW,     RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_RW,     RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_W,      RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        /* Other options */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, B->get_rankof(B, Bm, Bn),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, (L == 0) ? "ztplqs" : "ztplqt",
#endif
        0);
    (void)ldB;
    (void)ldA;

    (void)ib; (void)nb;
}
