/**
 *
 * @file starpu/codelet_zhe2ge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhe2ge StarPU codelet
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @date 2016-12-09
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zhe2ge_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int M;
    int N;
    const CHAMELEON_Complex64_t *A;
    int ldA;
    CHAMELEON_Complex64_t *B;
    int ldB;

    A = (const CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    ldA = STARPU_MATRIX_GET_LD( descr[0] );
    ldB = STARPU_MATRIX_GET_LD( descr[1] );

    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N);
    CORE_zhe2ge(uplo, M, N, A, ldA, B, ldB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhe2ge, 2, cl_zhe2ge_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zhe2ge(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int m, int n, int mb,
                       const CHAM_desc_t *A, int Am, int An, int ldA,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldB)
{
    (void)mb;
    struct starpu_codelet *codelet = &cl_zhe2ge;
    void (*callback)(void*) = options->profiling ? cl_zhe2ge_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,  &uplo,                sizeof(int),
        STARPU_VALUE,     &m,                        sizeof(int),
        STARPU_VALUE,     &n,                        sizeof(int),
        STARPU_R,             RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,             RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zhe2ge",
#endif
        0);
    (void)ldB;
    (void)ldA;
}
