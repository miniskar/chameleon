/**
 *
 * @file starpu/codelet_gemmex.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon gemmex StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2023-07-06
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelets.h"

CHAMELEON_CL_CB( gemmex, cti_handle_get_m(task->handles[2]), cti_handle_get_n(task->handles[2]), cti_handle_get_n(task->handles[0]), 2. *M*N*K) /* If A^t, computation is wrong */

struct cl_gemmex_args_s {
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    double alpha;
    double beta;
};

#if !defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA)
static void
cl_gemmex_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_gemmex_args_s *clargs = (struct cl_gemmex_args_s *)cl_arg;
    cublasHandle_t          handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;
    void *ptrAlpha, *ptrBeta;

    switch( tileC->flttype ) {
    case ChamRealHalf:
    {
        CHAMELEON_Real16_t halpha = clargs->alpha;
        CHAMELEON_Real16_t hbeta  = clargs->beta;
        ptrAlpha = &halpha;
        ptrBeta  = &hbeta;
    }
    break;
    case ChamRealFloat:
    {
        float salpha = clargs->alpha;
        float sbeta  = clargs->beta;
        ptrAlpha = &salpha;
        ptrBeta  = &sbeta;
    }
    break;
    case ChamRealDouble:
    {
        double dalpha = clargs->alpha;
        double dbeta  = clargs->beta;
        ptrAlpha = &dalpha;
        ptrBeta  = &dbeta;
    }
    break;
    case ChamComplexFloat:
    {
        CHAMELEON_Complex32_t calpha = clargs->alpha;
        CHAMELEON_Complex32_t cbeta  = clargs->beta;
        ptrAlpha = &calpha;
        ptrBeta  = &cbeta;
    }
    break;
    case ChamComplexDouble:
    {
        CHAMELEON_Complex64_t zalpha = clargs->alpha;
        CHAMELEON_Complex64_t zbeta  = clargs->beta;
        ptrAlpha = &zalpha;
        ptrBeta  = &zbeta;
    }
    break;
    default:
        fprintf( stderr, "cl_gemmex: Unknown C datatype\n" );
        return;
    }

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    assert( tileC->format & CHAMELEON_TILE_FULLRANK );

    CUDA_gemmex(
        clargs->transA, clargs->transB,
        clargs->m, clargs->n, clargs->k,
        ptrAlpha,
        tileA->mat, tileA->ld, tileA->flttype,
        tileB->mat, tileB->ld, tileB->flttype,
        ptrBeta,
        tileC->mat, tileC->ld, tileC->flttype,
        handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( gemmex, NULL, cl_gemmex_cuda_func, STARPU_CUDA_ASYNC )

void INSERT_TASK_gemmex_Astat( const RUNTIME_option_t *options,
                               cham_trans_t transA, cham_trans_t transB,
                               int m, int n, int k, int nb,
                               double alpha, const CHAM_desc_t *A, int Am, int An,
                                             const CHAM_desc_t *B, int Bm, int Bn,
                               double beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    /* if ( alpha == 0. ) { */
    /*     INSERT_TASK_hlascal( options, ChamUpperLower, m, n, nb, */
    /*                          beta, C, Cm, Cn ); */
    /*     return; */
    /* } */

    struct cl_gemmex_args_s  *clargs = NULL;
    void (*callback)(void*);
    int                      accessC;
    int                      exec    = 0;
    char                    *cl_name = "gemmex_Astat";

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
        clargs = malloc( sizeof( struct cl_gemmex_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->beta   = beta;
    }

    /* Callback for profiling information */
    callback = options->profiling ? cl_gemmex_callback : NULL;

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

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 3,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ),
                                      C->get_blktile( C, Cm, Cn ) );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_gemmex,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_gemmex_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamRealHalf, Am, An),
        STARPU_R, RTBLKADDR(B, ChamRealHalf, Bm, Bn),
        accessC,  RTBLKADDR(C, ChamRealHalf, Cm, Cn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_NODE,   A->get_rankof(A, Am, An),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif
        0 );
}

void INSERT_TASK_gemmex( const RUNTIME_option_t *options,
                         cham_trans_t transA, cham_trans_t transB,
                         int m, int n, int k, int nb,
                         double alpha, const CHAM_desc_t *A, int Am, int An,
                                       const CHAM_desc_t *B, int Bm, int Bn,
                         double beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    /* if ( alpha == 0. ) { */
    /*     INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb, */
    /*                          beta, C, Cm, Cn ); */
    /*     return; */
    /* } */

    struct cl_gemmex_args_s  *clargs = NULL;
    void (*callback)(void*);
    int                      accessC;
    int                      exec = 0;
    char                    *cl_name = "gemmex";

    if ( !(options->withcuda) ) {
        /* Fallback to cpu version */
        INSERT_TASK_gemm( options, transA, transB, m, n, k, nb,
                          alpha, A, Am, An, B, Bm, Bn, beta, C, Cm, Cn );
        return;
    }

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_gemmex_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->beta   = beta;
    }

    /* Callback for profiling information */
    callback = options->profiling ? cl_gemmex_callback : NULL;

    /* Reduce the C access if needed */
    accessC = ( beta == 0. ) ? STARPU_W : (STARPU_RW | ((beta == 1.) ? STARPU_COMMUTE : 0));

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 3,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ),
                                      C->get_blktile( C, Cm, Cn ) );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_gemmex,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_gemmex_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamRealHalf, Am, An),
        STARPU_R, RTBLKADDR(B, ChamRealHalf, Bm, Bn),
        accessC,  RTBLKADDR(C, ChamRealHalf, Cm, Cn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_POSSIBLY_PARALLEL, options->parallel,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif
        0 );
}
