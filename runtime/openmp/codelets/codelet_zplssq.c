/**
 *
 * @file openmp/codelet_zplssq.c
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
 * @author Philippe Virouleau
 * @date 2018-06-20
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
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
                         const CHAM_desc_t *IN,  int INm,  int INn,
                         const CHAM_desc_t *OUT, int OUTm, int OUTn )
{
    double *sclssq_in  = RTBLKADDR(IN,  double, INm,  INn );
    double *sclssq_out = RTBLKADDR(OUT, double, OUTm, OUTn);
#pragma omp task depend(in: sclssq_in[0]) depend(inout: sclssq_out[0])
    {
        if( sclssq_out[0] < sclssq_in[0] ) {
            sclssq_out[1] = sclssq_in[1]  + (sclssq_out[1] * (( sclssq_out[0] / sclssq_in[0] ) * ( sclssq_out[0] / sclssq_in[0] )));
            sclssq_out[0] = sclssq_in[0];
        } else {
            sclssq_out[1] = sclssq_out[1] + (sclssq_in[1]  * (( sclssq_in[0] / sclssq_out[0] ) * ( sclssq_in[0] / sclssq_out[0] )));
        }
    }
}

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    CHAMELEON_Complex64_t *res = RTBLKADDR(RESULT, CHAMELEON_Complex64_t, RESULTm, RESULTn);

#pragma omp task depend(inout: res[0])
    {
        res[0] = res[0] * sqrt( res[1] );
    }
}
