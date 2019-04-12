/**
 *
 * @file quark/codelet_zhessq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhessq Quark codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @date 2014-11-16
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zhessq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    INSERT_TASK_zsyssq( options, storev, uplo, n,
                        A, Am, An, lda,
                        SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
