/**
 *
 * @file starpu/codelet_zgeqrt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrt StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgeqrt_cpu_func(void *descr[], void *cl_arg)
{
    CHAMELEON_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileW;
    CHAMELEON_Complex64_t *TAU;
    CHAMELEON_Complex64_t *WORK;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileW = cti_interface_get(descr[2]); /* max(m,n) + ib * n */

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &h_work);

    TAU  = tileW->mat;
    WORK = TAU + chameleon_max( m, n );

    TCORE_zlaset( ChamUpperLower, ib, n, 0., 0., tileT );
    TCORE_zgeqrt(m, n, ib, tileA, tileT, TAU, WORK );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgeqrt, 3, cl_zgeqrt_cpu_func)

void INSERT_TASK_zgeqrt(const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgeqrt;
    void (*callback)(void*) = options->profiling ? cl_zgeqrt_callback : NULL;
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
        /* ib * (m+3*ib) + max(m,n) */
        STARPU_VALUE,    &h_work,            sizeof(CHAMELEON_starpu_ws_t *),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgeqrt",
#endif
        0);
}
