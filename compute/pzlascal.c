/**
 *
 * @file pzlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal parallel algorithm
 *
 * @version 0.9.2
 * @author Dalal Sukkari
 * @date 2016-11-30
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
/**
 *  Parallel scale of a matrix A
 */
void chameleon_pzlascal(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int tempmm, tempnn, tempmn, tempnm;
    int m, n;
    int ldam, ldan;
    int minmnt = chameleon_min(A->mt, A->nt);

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    switch(uplo) {
    case ChamLower:
        for (n = 0; n < minmnt; n++) {
            tempnm = n == A->mt-1 ? A->m-n*A->mb : A->mb;
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            ldan = BLKLDD(A, n);

            INSERT_TASK_zlascal(
                &options,
                ChamLower, tempnm, tempnn, A->mb,
                alpha, A(n, n), ldan);

            for (m = n+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-A->mb*m : A->nb;
                ldam = BLKLDD(A, m);

                INSERT_TASK_zlascal(
                    &options,
                    ChamUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n), ldam);
            }
        }
        break;

    case ChamUpper:
        for (m = 0; m < minmnt; m++) {
            tempmm = m == A->mt-1 ? A->m-A->mb*m : A->nb;
            tempmn = m == A->nt-1 ? A->n-m*A->nb : A->nb;
            ldam = BLKLDD(A, m);

            INSERT_TASK_zlascal(
                &options,
                ChamUpper, tempmm, tempmn, A->mb,
                alpha, A(m, m), ldam);

            for (n = m+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                INSERT_TASK_zlascal(
                    &options,
                    ChamUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n), ldam);
            }
        }
        break;

    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-A->mb*m : A->nb;
            ldam = BLKLDD(A, m);

            for (n = 0; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                INSERT_TASK_zlascal(
                    &options,
                    ChamUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n), ldam);
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
