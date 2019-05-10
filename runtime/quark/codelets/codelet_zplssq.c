/**
 *
 * @file quark/codelet_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq Quark codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zplssq_quark(Quark *quark)
{
    cham_store_t storev;
    int M;
    int N;
    double *SCLSSQ_IN;
    double *SCLSSQ_OUT;

    quark_unpack_args_5( quark, storev, M, N, SCLSSQ_IN, SCLSSQ_OUT );

    CORE_zplssq(storev, M, N, SCLSSQ_IN, SCLSSQ_OUT);
}

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 * @brief Compute sum( a_ij ^ 2 ) = scl * sqrt(ssq)
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
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn )
{
    int sizein = 2*M*N;
    int sizeout;

    if ( storev == ChamColumnwise ) {
        sizeout = 2*N;
    } else if ( storev == ChamRowwise ) {
        sizeout = 2*M;
    } else {
        sizeout = 2;
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_PLSSQ;
    QUARK_Insert_Task(opt->quark, CORE_zplssq_quark, (Quark_Task_Flags*)opt,
        sizeof(int),            &storev,    VALUE,
        sizeof(int),            &M,         VALUE,
        sizeof(int),            &N,         VALUE,
        sizeof(double)*sizein,  RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn), INPUT,
        sizeof(double)*sizeout, RTBLKADDR(SCLSSQ,     double, SCLSSQm,     SCLSSQn),     INOUT,
        0);
}


void CORE_zplssq2_quark(Quark *quark)
{
    int N;
    double *RESULT;

    quark_unpack_args_2( quark, N, RESULT );

    CORE_zplssq2(N, RESULT);
}

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_PLSSQ2;
    QUARK_Insert_Task(opt->quark, CORE_zplssq2_quark, (Quark_Task_Flags*)opt,
        sizeof(int),        &N,         VALUE,
        sizeof(double)*2*N, RTBLKADDR(RESULT, double, RESULTm, RESULTn), INOUT,
        0);
}
