/**
 *
 * @file starpu/codelet_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 *  INSERT_TASK_zplssq returns: scl * sqrt(ssq)
 *
 * with scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( 2*i )**2 * A( 2*i+1 ) )
 *                      i
 *
 * The values of A(2*i+1) are assumed to be at least unity.
 * The values of A(2*i) are assumed to be non-negative and scl is
 *
 *    scl = max( A( 2*i ) ),
 *           i
 *
 * The routine makes only one pass through the matrix A.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of couple (scale, sumsq) in the matrix A.
 *
 *  @param[in] A
 *          The 2-by-M matrix.
 *
 *  @param[out] result
 *          On exit, result contains scl * sqrt( ssq )
 *
 */
void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
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
static void cl_zplssq_cpu_func(void *descr[], void *cl_arg)
{
    double *SCLSSQ_IN;
    double *SCLSSQ_OUT;

    SCLSSQ_IN  = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    SCLSSQ_OUT = (double *)STARPU_MATRIX_GET_PTR(descr[1]);

    assert( SCLSSQ_OUT[0] >= 0. );
    if( SCLSSQ_OUT[0] < SCLSSQ_IN[0] ) {
        SCLSSQ_OUT[1] = SCLSSQ_IN[1]  + (SCLSSQ_OUT[1] * (( SCLSSQ_OUT[0] / SCLSSQ_IN[0] ) * ( SCLSSQ_OUT[0] / SCLSSQ_IN[0] )));
        SCLSSQ_OUT[0] = SCLSSQ_IN[0];
    } else {
        if ( SCLSSQ_OUT[0] > 0 ) {
            SCLSSQ_OUT[1] = SCLSSQ_OUT[1] + (SCLSSQ_IN[1]  * (( SCLSSQ_IN[0] / SCLSSQ_OUT[0] ) * ( SCLSSQ_IN[0] / SCLSSQ_OUT[0] )));
        }
    }

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq, 2, cl_zplssq_cpu_func)

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options,
                         const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    struct starpu_codelet *codelet = &cl_zplssq2;
    void (*callback)(void*) = options->profiling ? cl_zplssq2_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( RESULT, RESULTm, RESULTn );
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_RW, RTBLKADDR(RESULT, double, RESULTm, RESULTn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zplssq2",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zplssq2_cpu_func(void *descr[], void *cl_arg)
{
    double *RESULT;

    RESULT = (double *)STARPU_MATRIX_GET_PTR(descr[0]);

    RESULT[0] = RESULT[0] * sqrt( RESULT[1] );

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq2, 1, cl_zplssq2_cpu_func)
