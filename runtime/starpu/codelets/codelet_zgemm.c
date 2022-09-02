/**
 *
 * @file starpu/codelet_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm StarPU codelet
 *
 * @version 1.2.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Gwenole Lucas
 * @author Philippe Swartvagher
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

struct cl_zgemm_args_s {
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileC;
};

#if defined(CHAMELEON_USE_BUBBLE)
static inline int
cl_zgemm_is_bubble( struct starpu_task *t, void *_args )
{
    struct cl_zgemm_args_s *clargs = (struct cl_zgemm_args_s *)(t->cl_arg);
    (void)_args;

    return( ( clargs->tileA->format & CHAMELEON_TILE_DESC ) &&
            ( clargs->tileB->format & CHAMELEON_TILE_DESC ) &&
            ( clargs->tileC->format & CHAMELEON_TILE_DESC ) );
}

static void
cl_zgemm_bubble_func( struct starpu_task *t, void *_args )
{
    struct cl_zgemm_args_s *clargs = (struct cl_zgemm_args_s *)(t->cl_arg);
    bubble_args_t          *b_args  = (bubble_args_t *)_args;
    RUNTIME_request_t       request = RUNTIME_REQUEST_INITIALIZER;
    (void)_args;

    /* Register the task parent */
    request.parent = t;

#if defined(CHAMELEON_BUBBLE_PARALLEL_INSERT)
    request.dependency = t;
    starpu_task_end_dep_add( t, 1 );
#endif

    chameleon_pzgemm( NULL, clargs->transA, clargs->transB,
                      clargs->alpha, clargs->tileA->mat, clargs->tileB->mat,
                      clargs->beta,  clargs->tileC->mat,
                      b_args->sequence, &request );

    free( _args );
}
#endif /* defined(CHAMELEON_USE_BUBBLE) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgemm_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zgemm_args_s *clargs = (struct cl_zgemm_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    TCORE_zgemm( clargs->transA, clargs->transB,
                 clargs->m, clargs->n, clargs->k,
                 clargs->alpha, tileA, tileB,
                 clargs->beta,  tileC );
}

#ifdef CHAMELEON_USE_CUDA
static void
cl_zgemm_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_zgemm_args_s *clargs = (struct cl_zgemm_args_s *)cl_arg;
    cublasHandle_t          handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    assert( tileC->format & CHAMELEON_TILE_FULLRANK );

    CUDA_zgemm(
        clargs->transA, clargs->transB,
        clargs->m, clargs->n, clargs->k,
        (cuDoubleComplex*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        (cuDoubleComplex*)&(clargs->beta),
        tileC->mat, tileC->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zgemm, cl_zgemm_cpu_func, cl_zgemm_cuda_func, STARPU_CUDA_ASYNC )

void INSERT_TASK_zgemm_Astat( const RUNTIME_option_t *options,
                              cham_trans_t transA, cham_trans_t transB,
                              int m, int n, int k, int nb,
                              CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                           const CHAM_desc_t *B, int Bm, int Bn,
                              CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                                    beta, C, Cm, Cn );
    }

    struct cl_zgemm_args_s  *clargs = NULL;
    void (*callback)(void*);
    int                      accessC;
    int                      exec    = 0;
    char                    *cl_name = "zgemm_Astat";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
     /* Check A as write, since it will be the owner of the computation */
    CHAMELEON_ACCESS_W(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
     /* Check C as read, since it will be used in a reduction */
    CHAMELEON_ACCESS_R(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zgemm_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->tileA  = A->get_blktile( A, Am, An );
        clargs->tileB  = B->get_blktile( B, Bm, Bn );
        clargs->beta   = beta;
        clargs->tileC  = C->get_blktile( C, Cm, Cn );
    }

    /* Callback for profiling information */
    callback = options->profiling ? cl_zgemm_callback : NULL;

    /* Reduce the C access if needed */
    if ( beta == 0. ) {
        accessC = STARPU_W;
    }
#if defined(HAVE_STARPU_MPI_REDUX)
    else if ( beta == 1. ) {
        accessC = STARPU_MPI_REDUX;
    }
#endif
    else {
        accessC = STARPU_RW;
    }

#if defined(CHAMELEON_KERNELS_TRACE)
    {
        char *cl_fullname;
        chameleon_asprintf( &cl_fullname, "%s( %s, %s, %s )", cl_name, clargs->tileA->name, clargs->tileB->name, clargs->tileC->name );
        cl_name = cl_fullname;
    }
#endif

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zgemm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgemm_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R, RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        accessC,  RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_NODE,   A->get_rankof(A, Am, An),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif
        0 );
}

void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                                    beta, C, Cm, Cn );
    }

    struct cl_zgemm_args_s  *clargs = NULL;
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    int                      is_bubble, accessC;
    int                      exec = 0;
    char                    *cl_name = "zgemm";
    bubble_args_t           *b_args = NULL;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zgemm_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->tileA  = A->get_blktile( A, Am, An );
        clargs->tileB  = B->get_blktile( B, Bm, Bn );
        clargs->beta   = beta;
        clargs->tileC  = C->get_blktile( C, Cm, Cn );
    }

    /* Callback for profiling information */
    callback = options->profiling ? cl_zgemm_callback : NULL;

    /* Reduce the C access if needed */
    accessC = ( beta == 0. ) ? STARPU_W : (STARPU_RW | ((beta == 1.) ? STARPU_COMMUTE : 0));

    /* Check if this is a bubble */
    is_bubble = ( ( clargs->tileA->format & CHAMELEON_TILE_DESC ) &&
                  ( clargs->tileB->format & CHAMELEON_TILE_DESC ) &&
                  ( clargs->tileC->format & CHAMELEON_TILE_DESC ) );
    if ( is_bubble ) {
        b_args = malloc( sizeof(bubble_args_t) + sizeof(struct cl_zgemm_args_s) );
        b_args->sequence = options->sequence;
        b_args->parent   = request->parent;
        cl_name = "zgemm_bubble";
    }

#if defined(CHAMELEON_KERNELS_TRACE)
    {
        char *cl_fullname;
        chameleon_asprintf( &cl_fullname, "%s( %s, %s, %s )", cl_name, clargs->tileA->name, clargs->tileB->name, clargs->tileC->name );
        cl_name = cl_fullname;
    }
#endif

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zgemm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgemm_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R, RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        accessC,  RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),

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
        STARPU_BUBBLE_GEN_DAG_FUNC,     cl_zgemm_bubble_func,
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
}
