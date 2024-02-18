/**
 *
 * @file pzgerst.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgerst parallel algorithm
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2023-07-06
 * @precisions normal z -> d
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n

void chameleon_pzgerst( cham_uplo_t         uplo,
                        CHAM_desc_t        *A,
                        RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for(m = 0; m < A->mt; m++) {
        int tempmm = ( m == (A->mt-1) ) ? A->m - m * A->mb : A->mb;
        int nmin   = ( uplo == ChamUpper ) ? m                         : 0;
        int nmax   = ( uplo == ChamLower ) ? chameleon_min(m+1, A->nt) : A->nt;

        for(n = nmin; n < nmax; n++) {
            CHAM_tile_t *tile = A->get_blktile( A, m, n );

            if (( tile->rank == A->myrank ) &&
                ( tile->flttype != ChamComplexDouble ) )
            {
                int tempnn = ( n == (A->nt-1) ) ? A->n - n * A->nb : A->nb;

                INSERT_TASK_zgerst( &options,
                                     tempmm, tempnn, A( m, n ) );
            }
        }
    }

    RUNTIME_options_finalize(&options, chamctxt);
}
