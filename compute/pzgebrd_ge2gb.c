/**
 *
 * @file pzgebrd_ge2gb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgebrd_ge2gb parallel algorithm
 *
 * @version 0.9.2
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2016-12-09
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

void chameleon_pzgebrd_ge2gb( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    int k;
    int tempkm, tempkn;
    CHAM_desc_t *A1, *A2, *T1, *D1 = NULL;

    if (A->m >= A->n){
        for (k = 0; k < A->nt; k++) {
            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

            A1 = chameleon_desc_submatrix(A, k*A->mb,     k*A->nb, A->m-k*A->mb, tempkn);
            A2 = chameleon_desc_submatrix(A, k*A->mb, (k+1)*A->nb, A->m-k*A->mb, A->n-(k+1)*A->nb);
            T1 = chameleon_desc_submatrix(T, k*T->mb,     k*T->nb, T->m-k*T->mb, T->nb );
            if ( D != NULL ) {
                D1 = chameleon_desc_submatrix(D, k*D->mb, k*D->nb, D->m-k*D->mb, tempkn);
            }

            chameleon_pzgeqrf( genD, A1, T1, D1,
                               sequence, request);

            chameleon_pzunmqr( 0, ChamLeft, ChamConjTrans,
                               A1, A2, T1, D1,
                               sequence, request);

            if (k+1 < A->nt){
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                A1 = chameleon_desc_submatrix(A,     k*A->mb, (k+1)*A->nb, tempkm,           A->n-(k+1)*A->nb);
                A2 = chameleon_desc_submatrix(A, (k+1)*A->mb, (k+1)*A->nb, A->m-(k+1)*A->mb, A->n-(k+1)*A->nb);
                T1 = chameleon_desc_submatrix(T,     k*T->mb, (k+1)*T->nb, T->mb,            T->n-(k+1)*T->nb);
                if ( D != NULL ) {
                    D1 = chameleon_desc_submatrix(D, k*D->mb, (k+1)*D->nb, tempkm,           D->n-(k+1)*D->nb);
                }

                chameleon_pzgelqf( genD, A1, T1, D1,
                                   sequence, request);

                chameleon_pzunmlq( 0, ChamRight, ChamConjTrans,
                                   A1, A2, T1, D1,
                                   sequence, request);
            }
        }
    }
    else{
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

            A1 = chameleon_desc_submatrix(A,     k*A->mb, k*A->nb, tempkm,           A->n-k*A->nb);
            A2 = chameleon_desc_submatrix(A, (k+1)*A->mb, k*A->nb, A->m-(k+1)*A->mb, A->n-k*A->nb);
            T1 = chameleon_desc_submatrix(T,     k*T->mb, k*T->nb, T->mb,            T->n-k*T->nb);
            if ( D != NULL ) {
                D1 = chameleon_desc_submatrix(D, k*D->mb, k*D->nb, tempkm,           D->n-k*D->nb);
            }
            chameleon_pzgelqf( genD, A1, T1, D1,
                               sequence, request);

            chameleon_pzunmlq( 0, ChamRight, ChamConjTrans,
                               A1, A2, T1, D1,
                               sequence, request);

            if (k+1 < A->mt){
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                A1 = chameleon_desc_submatrix(A, (k+1)*A->mb,     k*A->nb, A->m-(k+1)*A->mb, tempkn);
                A2 = chameleon_desc_submatrix(A, (k+1)*A->mb, (k+1)*A->nb, A->m-(k+1)*A->mb, A->n-(k+1)*A->nb);
                T1 = chameleon_desc_submatrix(T, (k+1)*T->mb,     k*T->nb, T->m-(k+1)*T->mb, T->nb );
                if ( D != NULL ) {
                    D1 = chameleon_desc_submatrix(D, (k+1)*D->mb, k*D->nb, D->m-(k+1)*D->mb, tempkn);
                }

                chameleon_pzgeqrf( genD, A1, T1, D1,
                                   sequence, request);

                chameleon_pzunmqr( 0, ChamLeft, ChamConjTrans,
                                   A1, A2, T1, D1,
                                   sequence, request);
            }
        }
    }
}
