/**
 *
 * @file starpu/codelet_zlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum StarPU codelet
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
 * @author Gwenole Lucas
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

struct cl_zlauum_args_s {
    cham_uplo_t uplo;
    int n;
    CHAM_tile_t *tileA;
};

#if defined(CHAMELEON_USE_BUBBLE)
static inline int
cl_zlauum_is_bubble( struct starpu_task *t, void *_args )
{
    struct cl_zlauum_args_s *clargs = (struct cl_zlauum_args_s *)(t->cl_arg);
    (void)_args;

    return( clargs->tileA->format & CHAMELEON_TILE_DESC );
}

static void
cl_zlauum_bubble_func( struct starpu_task *t, void *_args )
{
    struct cl_zlauum_args_s *clargs  = (struct cl_zlauum_args_s *)(t->cl_arg);
    bubble_args_t           *b_args  = (bubble_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    /* Register the task parent */
    request.parent = t;

#if defined(CHAMELEON_BUBBLE_PARALLEL_INSERT)
    request.dependency = t;
    starpu_task_end_dep_add( t, 1 );
#endif

    chameleon_pzlauum( clargs->uplo, clargs->tileA->mat,
                       b_args->sequence, &request );

    free( _args );
}
#endif /* defined(CHAMELEON_USE_BUBBLE) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlauum_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlauum_args_s *clargs = (struct cl_zlauum_args_s *)cl_arg;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zlauum( clargs->uplo, clargs->n, tileA );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlauum, cl_zlauum_cpu_func )

void INSERT_TASK_zlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An )
{
    struct cl_zlauum_args_s *clargs = NULL;
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    int                      is_bubble;
    int                      exec = 0;
    char                    *cl_name = "zlauum";
    bubble_args_t           *b_args = NULL;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zlauum_args_s ) );
        clargs->uplo  = uplo;
        clargs->n     = n;
        clargs->tileA = A->get_blktile( A, Am, An );
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlauum_callback : NULL;

    /* Check if this is a bubble */
    is_bubble = ( clargs->tileA->format & CHAMELEON_TILE_DESC );
    if ( is_bubble ) {
        b_args = malloc( sizeof(bubble_args_t) + sizeof(struct cl_zlauum_args_s) );
        b_args->sequence = options->sequence;
        b_args->parent   = request->parent;
        memcpy( &(b_args->clargs), clargs, sizeof(struct cl_zlauum_args_s) );
        cl_name = "zlauum_bubble";
    }

#if defined(CHAMELEON_KERNELS_TRACE)
    {
        char *cl_fullname;
        chameleon_asprintf( &cl_fullname, "%s( %s )", cl_name, clargs->tileA->name );
        cl_name = cl_fullname;
    }
#endif

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlauum,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlauum_args_s),
        STARPU_RW,     RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        /* Bubble management */
#if defined(CHAMELEON_USE_BUBBLE)
        STARPU_BUBBLE_FUNC,             is_bubble_func,
        STARPU_BUBBLE_FUNC_ARG,         b_args,
        STARPU_BUBBLE_GEN_DAG_FUNC,     cl_zlauum_bubble_func,
        STARPU_BUBBLE_GEN_DAG_FUNC_ARG, b_args,

#if defined(CHAMELEON_BUBBLE_PROFILE)
        STARPU_BUBBLE_PARENT, request->parent,
#endif

#if defined(CHAMELEON_BUBBLE_PARALLEL_INSERT)
        STARPU_CALLBACK_WITH_ARG_NFREE, callback_end_dep_release, request->dependency,
#endif
#endif
        0 );

    /* Dependency is used only by the first submitted task and should not be reused */
    request->dependency = NULL;
    (void)nb;
}
