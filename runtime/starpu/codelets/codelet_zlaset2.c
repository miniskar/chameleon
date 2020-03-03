/**
 *
 * @file starpu/codelet_zlaset2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset2 StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
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
static void cl_zlaset2_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &alpha);
    TCORE_zlaset2(uplo, M, N, alpha, tileA);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlaset2, 1, cl_zlaset2_cpu_func)

void INSERT_TASK_zlaset2(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int M, int N,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An)
{

    struct starpu_codelet *codelet = &cl_zlaset2;
    void (*callback)(void*) = options->profiling ? cl_zlaset2_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_W(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,  &uplo,                sizeof(int),
        STARPU_VALUE,     &M,                        sizeof(int),
        STARPU_VALUE,     &N,                        sizeof(int),
        STARPU_VALUE, &alpha,         sizeof(CHAMELEON_Complex64_t),
        STARPU_W,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zlaset2",
#endif
        0);
}
