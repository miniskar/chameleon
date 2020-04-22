/**
 *
 * @file starpu/codelet_zlatro.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
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
static void cl_zlatro_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int M;
    int N;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &M, &N);
    TCORE_zlatro(uplo, trans, M, N, tileA, tileB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlatro, 2, cl_zlatro_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    struct starpu_codelet *codelet = &cl_zlatro;
    void (*callback)(void*) = NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,   &uplo,    sizeof(int),
        STARPU_VALUE,   &trans,   sizeof(int),
        STARPU_VALUE,   &m,       sizeof(int),
        STARPU_VALUE,   &n,       sizeof(int),
        STARPU_R,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,        RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zlatro",
#endif
        0);
    (void)mb;
}
