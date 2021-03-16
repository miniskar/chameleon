/**
 *
 * @file starpu/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy StarPU codelet
 *
 * @version 1.0.0
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

struct cl_zlacpy_args_s {
    cham_uplo_t uplo;
    int m;
    int n;
    int displA;
    int displB;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlacpy_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s clargs;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &clargs );
    assert( clargs.displA == 0 );
    assert( clargs.displB == 0 );
    /* A = tileA->mat; */
    /* B = tileB->mat; */
    /* CORE_zlacpy( uplo, M, N, A + displA, tileA->ld, B + displB, tileB->ld ); */
    TCORE_zlacpy( clargs.uplo, clargs.m, clargs.n, tileA, tileB );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlacpy, cl_zlacpy_cpu_func )

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int nb,
                          int displA, const CHAM_desc_t *A, int Am, int An,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn )
{
    struct cl_zlacpy_args_s clargs = {
        .uplo   = uplo,
        .m      = m,
        .n      = n,
        .displA = displA,
        .displB = displB,
        .tileA  = A->get_blktile( A, Am, An ),
        .tileB  = B->get_blktile( B, Bm, Bn ),
    };
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    starpu_option_request_t *schedopt = (starpu_option_request_t *)(request->schedopt);
    int                      workerid;
    char                    *cl_name = "zlacpy";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlacpy_callback : NULL;

    /* Fix the worker id */
    workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlacpy,
        /* Task codelet arguments */
        STARPU_VALUE, &clargs, sizeof(struct cl_zlacpy_args_s),
        STARPU_R,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,      RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );

    (void)nb;
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    INSERT_TASK_zlacpyx( options, uplo, m, n, nb,
                         0, A, Am, An,
                         0, B, Bm, Bn );
}
