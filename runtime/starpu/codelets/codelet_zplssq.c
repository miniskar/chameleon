/**
 *
 * @file starpu/codelet_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplssq_cpu_func(void *descr[], void *cl_arg)
{
    cham_store_t storev;
    int M;
    int N;
    double *SCLSSQ_IN;
    double *SCLSSQ_OUT;

    starpu_codelet_unpack_args(cl_arg, &storev, &M, &N);
    SCLSSQ_IN  = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    SCLSSQ_OUT = (double *)STARPU_MATRIX_GET_PTR(descr[1]);

    CORE_zplssq(storev, M, N, SCLSSQ_IN, SCLSSQ_OUT);

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq, 2, cl_zplssq_cpu_func)

void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCLSSQ_IN,  int SCLSSQ_INm,  int SCLSSQ_INn,
                         const CHAM_desc_t *SCLSSQ_OUT, int SCLSSQ_OUTm, int SCLSSQ_OUTn )
{
    struct starpu_codelet *codelet = &cl_zplssq;
    void (*callback)(void*) = options->profiling ? cl_zplssq_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(  SCLSSQ_IN,  SCLSSQ_INm,  SCLSSQ_INn  );
    CHAMELEON_ACCESS_RW( SCLSSQ_OUT, SCLSSQ_OUTm, SCLSSQ_OUTn );
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &storev,            sizeof(int),
        STARPU_VALUE,    &M,                 sizeof(int),
        STARPU_VALUE,    &N,                 sizeof(int),
        STARPU_R,  RTBLKADDR( SCLSSQ_IN,  double, SCLSSQ_INm,  SCLSSQ_INn  ),
        STARPU_RW, RTBLKADDR( SCLSSQ_OUT, double, SCLSSQ_OUTm, SCLSSQ_OUTn ),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zplssq",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplssq2_cpu_func(void *descr[], void *cl_arg)
{
    int N;
    double *RESULT;

    starpu_codelet_unpack_args(cl_arg, &N);
    RESULT = (double *)STARPU_MATRIX_GET_PTR(descr[0]);

    CORE_zplssq2(N, RESULT);

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq2, 1, cl_zplssq2_cpu_func)

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    struct starpu_codelet *codelet = &cl_zplssq2;
    void (*callback)(void*) = options->profiling ? cl_zplssq2_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( RESULT, RESULTm, RESULTn );
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &N,                 sizeof(int),
        STARPU_RW, RTBLKADDR(RESULT, double, RESULTm, RESULTn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zplssq2",
#endif
        0);
}
