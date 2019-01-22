/**
 *
 * @file parsec/codelet_zplssq.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

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
static inline int
CORE_zplssq_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    double *SCLSSQ_IN;
    double *SCLSSQ_OUT;

    parsec_dtd_unpack_args(
        this_task, &SCLSSQ_IN, &SCLSSQ_OUT );

    assert( SCLSSQ_OUT[0] >= 0. );
    if( SCLSSQ_OUT[0] < SCLSSQ_IN[0] ) {
        SCLSSQ_OUT[1] = SCLSSQ_IN[1]  + (SCLSSQ_OUT[1] * (( SCLSSQ_OUT[0] / SCLSSQ_IN[0] ) * ( SCLSSQ_OUT[0] / SCLSSQ_IN[0] )));
        SCLSSQ_OUT[0] = SCLSSQ_IN[0];
    } else {
        if ( SCLSSQ_OUT[0] > 0 ) {
            SCLSSQ_OUT[1] = SCLSSQ_OUT[1] + (SCLSSQ_IN[1]  * (( SCLSSQ_IN[0] / SCLSSQ_OUT[0] ) * ( SCLSSQ_IN[0] / SCLSSQ_OUT[0] )));
        }
    }

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                        const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zplssq_parsec, options->priority, "plssq",
        PASSED_BY_REF,         RTBLKADDR( SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn ),    INPUT,
        PASSED_BY_REF,         RTBLKADDR( SCLSSQ, double, SCLSSQm, SCLSSQn ),                INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}

static inline int
CORE_zplssq2_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    double *RESULT;

    parsec_dtd_unpack_args(
        this_task, &RESULT );

    RESULT[0] = RESULT[0] * sqrt( RESULT[1] );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options,
                         const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zplssq2_parsec, options->priority, "plssq2",
        PASSED_BY_REF,         RTBLKADDR( RESULT, double, RESULTm, RESULTn ),     INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
