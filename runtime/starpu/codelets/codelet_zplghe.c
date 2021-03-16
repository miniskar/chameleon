/**
 *
 * @file starpu/codelet_zplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe StarPU codelet
 *
 * @version 1.0.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

struct cl_zplghe_args_s {
    double bump;
    int m;
    int n;
    CHAM_tile_t *tileA;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplghe_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zplghe_args_s clargs;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    starpu_codelet_unpack_args( cl_arg, &clargs );
    TCORE_zplghe( clargs.bump, clargs.m, clargs.n, tileA,
                  clargs.bigM, clargs.m0, clargs.n0, clargs.seed );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zplghe, cl_zplghe_cpu_func )

void INSERT_TASK_zplghe( const RUNTIME_option_t *options,
                         double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed )
{
    struct cl_zplghe_args_s clargs = {
        .bump  = bump,
        .m     = m,
        .n     = n,
        .tileA = A->get_blktile( A, Am, An ),
        .bigM  = bigM,
        .m0    = m0,
        .n0    = n0,
        .seed  = seed,
    };
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    starpu_option_request_t *schedopt = (starpu_option_request_t *)(request->schedopt);
    int                      workerid;
    char                    *cl_name = "zplghe";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_W(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zplghe_callback : NULL;

    /* Fix the worker id */
    workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zplghe,
        /* Task codelet arguments */
        STARPU_VALUE, &clargs, sizeof(struct cl_zplghe_args_s),
        STARPU_W,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );

}
