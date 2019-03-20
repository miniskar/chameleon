/**
 *
 * @file pzlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c parallel algorithm
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2014-11-16
 * @precisions mixed zc -> ds
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define SA(m,n) SA,  m,  n
#define SB(m,n) SB,  m,  n
/**
 *
 */
/**
 *
 */
void chameleon_pclag2z(CHAM_desc_t *SA, CHAM_desc_t *B,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int X, Y;
    int m, n;
    int ldam, ldbm;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for(m = 0; m < SA->mt; m++) {
        X = m == SA->mt-1 ? SA->m-m*SA->mb : SA->mb;
        ldam = BLKLDD(SA, m);
        ldbm = BLKLDD(B, m);
        for(n = 0; n < SA->nt; n++) {
            Y = n == SA->nt-1 ? SA->n-n*SA->nb : SA->nb;
            INSERT_TASK_clag2z(
                &options,
                X, Y, SA->mb,
                SA(m, n), ldam,
                B(m, n), ldbm);
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
