/**
 *
 * @file openmp/codelet_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2018-06-15
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *IN,  int INm,  int INn,
                         const CHAM_desc_t *OUT, int OUTm, int OUTn )
{
    double *sclssq_in  = RTBLKADDR(IN,  double, INm,  INn );
    double *sclssq_out = RTBLKADDR(OUT, double, OUTm, OUTn);
#pragma omp task firstprivate(storev, M, N) depend(in: sclssq_in[0]) depend(inout: sclssq_out[0])
    CORE_zplssq(storev, M, N, sclssq_in, sclssq_out);
}

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    double *res = RTBLKADDR(RESULT, double, RESULTm, RESULTn);

#pragma omp task firstprivate(N) depend(inout: res[0])
    CORE_zplssq2(N, res);
}
