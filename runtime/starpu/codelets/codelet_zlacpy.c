/**
 *
 * @file starpu/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy StarPU codelet
 *
 * @version 1.2.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Alycia Lisito
 * @date 2022-02-22
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
    int lda;
    int ldb;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlacpy_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( clargs->displA == 0 );
    assert( clargs->displB == 0 );

    TCORE_zlacpy( clargs->uplo, clargs->m, clargs->n, tileA, tileB );
}

static void
cl_zlacpyx_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_zlacpyx( clargs->uplo, clargs->m, clargs->n, clargs->displA,
                   tileA, clargs->lda, clargs->displB, tileB, clargs->ldb );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlacpy,  cl_zlacpy_cpu_func  )
CODELETS_CPU( zlacpyx, cl_zlacpyx_cpu_func )

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    struct cl_zlacpy_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      exec = 0;
    char                    *cl_name = "zlacpyx";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zlacpy_args_s ) );
        clargs->uplo   = uplo;
        clargs->m      = m;
        clargs->n      = n;
        clargs->displA = displA;
        clargs->displB = displB;
        clargs->tileA  = A->get_blktile( A, Am, An );
        clargs->tileB  = B->get_blktile( B, Bm, Bn );
        clargs->lda    = lda;
        clargs->ldb    = ldb;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlacpyx_callback : NULL;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlacpyx,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlacpy_args_s),
        STARPU_R,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,      RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    struct cl_zlacpy_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      exec    = 0;
    char                    *cl_name = "zlacpy";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zlacpy_args_s ) );
        clargs->uplo   = uplo;
        clargs->m      = m;
        clargs->n      = n;
        clargs->displA = 0;
        clargs->displB = 0;
        clargs->tileA  = A->get_blktile( A, Am, An );
        clargs->tileB  = B->get_blktile( B, Bm, Bn );
        clargs->lda    = clargs->tileA->ld;
        clargs->ldb    = clargs->tileB->ld;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlacpy_callback : NULL;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlacpy,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlacpy_args_s),
        STARPU_R,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,      RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );
}
