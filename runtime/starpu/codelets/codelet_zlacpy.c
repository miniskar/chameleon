/**
 *
 * @file starpu/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zlacpyx(const RUNTIME_option_t *options,
                        cham_uplo_t uplo, int m, int n, int nb,
                        int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                        int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zlacpy;
    void (*callback)(void*) = options->profiling ? cl_zlacpy_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,   &uplo,                sizeof(int),
        STARPU_VALUE,   &m,                   sizeof(int),
        STARPU_VALUE,   &n,                   sizeof(int),
        STARPU_VALUE,   &displA,              sizeof(int),
        STARPU_R,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,   &lda,                 sizeof(int),
        STARPU_VALUE,   &displB,              sizeof(int),
        STARPU_W,        RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_VALUE,   &ldb,                 sizeof(int),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zlacpy",
#endif
        0);
}

void INSERT_TASK_zlacpy(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    INSERT_TASK_zlacpyx( options, uplo, m, n, nb,
                        0, A, Am, An, lda,
                        0, B, Bm, Bn, ldb );
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlacpy_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int M;
    int N;
    int displA;
    int displB;
    const CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;

    A = (const CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (CHAMELEON_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &displA, &LDA, &displB, &LDB);
    CORE_zlacpy(uplo, M, N, A + displA, LDA, B + displB, LDB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlacpy, 2, cl_zlacpy_cpu_func)
