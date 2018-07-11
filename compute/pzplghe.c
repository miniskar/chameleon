/**
 *
 * @file pzplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe parallel algorithm
 *
 * @version 1.0.0
 * @comment This file is a copy from pzplghe.c
 *          wich has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Rade Mathis
 * @author Florent Pruvost
 * @date 2016-08-01
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  morse_pzplghe - Generate a random hermitian (positive definite if 'bump' is large enough) half-matrix by tiles.
 */
void morse_pzplghe( double bump, cham_uplo_t uplo, CHAM_desc_t *A,
                    unsigned long long int seed,
                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *morse;
    RUNTIME_option_t options;

    int m, n;
    int ldam;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < A->mt; m++) {
    	tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
    	ldam = BLKLDD(A, m);

        /*
         *  ChamLower
         */
        if (uplo == ChamLower) {
            for (n = 0; n <= m; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n), ldam,
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
        /*
         * ChamUpper
         */
        else if (uplo == ChamUpper) {
            for (n = m; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n), ldam,
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
        /*
         * ChamUpperLower
         */
        else {
            for (n = 0; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n), ldam,
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
