/**
 *
 * @file pzgebrd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgebrd parallel algorithm
 *
 * @version 1.2.0
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-22
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/lapacke.h"
#endif

void chameleon_pzgebrd_ge2gb( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    int k;
    int tempkm, tempkn;
    CHAM_desc_t *A1, *A2, *T1, *D1 = NULL;

    if ( A->m >= A->n ){
        for ( k = 0; k < A->nt; k++ ) {
            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

            A1 = chameleon_desc_submatrix( A, k*A->mb,     k*A->nb, A->m-k*A->mb, tempkn );
            A2 = chameleon_desc_submatrix( A, k*A->mb, (k+1)*A->nb, A->m-k*A->mb, A->n-(k+1)*A->nb );
            T1 = chameleon_desc_submatrix( T, k*T->mb,     k*T->nb, T->m-k*T->mb, T->nb );
            if ( D != NULL ) {
                D1 = chameleon_desc_submatrix( D, k*D->mb, k*D->nb, D->m-k*D->mb, tempkn );
            }

            chameleon_pzgeqrf( genD, A1, T1, D1, sequence, request );
            chameleon_pzunmqr( 0, ChamLeft, ChamConjTrans, A1, A2, T1, D1, sequence, request );

            free( A1 );
            free( A2 );
            free( T1 );
            if ( D != NULL ) {
                free( D1 );
            }

            if ( k+1 < A->nt ) {
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                A1 = chameleon_desc_submatrix( A,     k*A->mb, (k+1)*A->nb, tempkm,           A->n-(k+1)*A->nb );
                A2 = chameleon_desc_submatrix( A, (k+1)*A->mb, (k+1)*A->nb, A->m-(k+1)*A->mb, A->n-(k+1)*A->nb );
                T1 = chameleon_desc_submatrix( T,     k*T->mb, (k+1)*T->nb, T->mb,            T->n-(k+1)*T->nb );
                if ( D != NULL ) {
                    D1 = chameleon_desc_submatrix( D, k*D->mb, (k+1)*D->nb, tempkm, D->n-(k+1)*D->nb );
                }

                chameleon_pzgelqf( genD, A1, T1, D1, sequence, request );
                chameleon_pzunmlq( 0, ChamRight, ChamConjTrans, A1, A2, T1, D1, sequence, request );

                free( A1 );
                free( A2 );
                free( T1 );
                if ( D != NULL ) {
                    free( D1 );
                }
            }
        }
    }
    else {
        for ( k = 0; k < A->mt; k++ ) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

            A1 = chameleon_desc_submatrix( A,     k*A->mb, k*A->nb, tempkm,           A->n-k*A->nb );
            A2 = chameleon_desc_submatrix( A, (k+1)*A->mb, k*A->nb, A->m-(k+1)*A->mb, A->n-k*A->nb );
            T1 = chameleon_desc_submatrix( T,     k*T->mb, k*T->nb, T->mb,            T->n-k*T->nb );
            if ( D != NULL ) {
                D1 = chameleon_desc_submatrix( D, k*D->mb, k*D->nb, tempkm,           D->n-k*D->nb );
            }

            chameleon_pzgelqf( genD, A1, T1, D1, sequence, request );
            chameleon_pzunmlq( 0, ChamRight, ChamConjTrans, A1, A2, T1, D1, sequence, request );

            free( A1 );
            free( A2 );
            free( T1 );
            if ( D != NULL ) {
                free( D1 );
            }

            if ( k+1 < A->mt ) {
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                A1 = chameleon_desc_submatrix( A, (k+1)*A->mb,     k*A->nb, A->m-(k+1)*A->mb, tempkn );
                A2 = chameleon_desc_submatrix( A, (k+1)*A->mb, (k+1)*A->nb, A->m-(k+1)*A->mb, A->n-(k+1)*A->nb );
                T1 = chameleon_desc_submatrix( T, (k+1)*T->mb,     k*T->nb, T->m-(k+1)*T->mb, T->nb );
                if ( D != NULL ) {
                    D1 = chameleon_desc_submatrix( D, (k+1)*D->mb, k*D->nb, D->m-(k+1)*D->mb, tempkn );
                }

                chameleon_pzgeqrf( genD, A1, T1, D1, sequence, request );
                chameleon_pzunmqr( 0, ChamLeft, ChamConjTrans, A1, A2, T1, D1, sequence, request );

                free( A1 );
                free( A2 );
                free( T1 );
                if ( D != NULL ) {
                    free( D1 );
                }
            }
        }
    }

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( T, sequence );
    if ( D != NULL ) {
        CHAMELEON_Desc_Flush( D, sequence );
    }
}

int chameleon_pzgebrd_gb2bd( cham_job_t jobu, cham_job_t jobvt, CHAM_desc_t *A,
                             CHAMELEON_Complex64_t *U, int LDU,
                             CHAMELEON_Complex64_t *VT, int LDVT,
                             double *E, double *S,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t descAB;
    cham_uplo_t uplo;
    int M, N, MINMN, NB, LDAB, ABn;
    int info;
    int KL, KU;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return sequence->status;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    M     = A->m;
    N     = A->n;
    MINMN = chameleon_min(M, N);
    NB    = A->mb;
    LDAB  = NB + 1;
    uplo  = M >= N ? ChamUpper : ChamLower;
    ABn   = MINMN;

    /* Allocate band structure */
    chameleon_zdesc_alloc( descAB, LDAB, NB, /* mb, nb */
                           LDAB, ABn,        /* lm, ln */
                           0, 0,             /* i,  j  */
                           LDAB, ABn,        /* m,  n  */
                           NULL );

    /* Convert matrix to band form */
    chameleon_pztile2band( uplo, A, &descAB, sequence, request );

    /* NCC = 0, C = NULL, we do not update any matrix with new singular vectors */
    /* On exit, AB = U (S +~ E) VT */
    KL = uplo == ChamUpper ? 0  : NB;
    KU = uplo == ChamUpper ? NB : 0;

    /* Manage the case where only singular values are required */
    char gbbrd_vect;
    if ( jobu == ChamNoVec ) {
        if ( jobvt == ChamNoVec ) {
            gbbrd_vect = 'N';
        }
        else {
            gbbrd_vect = 'P';
        }
    }
    else {
        if ( jobvt == ChamNoVec ) {
            gbbrd_vect = 'Q';
        }
        else {
            gbbrd_vect = 'B';
        }
    }

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( &descAB, sequence );
    chameleon_sequence_wait( chamctxt, sequence );

#if !defined(CHAMELEON_SIMULATION)
    info = LAPACKE_zgbbrd( LAPACK_COL_MAJOR, gbbrd_vect, M, N, 0, KL, KU,
                           (CHAMELEON_Complex64_t *) descAB.mat, LDAB, S, E,
                           U, LDU, VT, LDVT, NULL, 1 );
    if ( info != 0 ) {
        fprintf( stderr, "CHAMELEON_zgesvd_Tile_Async: LAPACKE_zgbbrd = %d\n", info );
    }
    assert( info == 0 );
#endif /* !defined(CHAMELEON_SIMULATION) */

    chameleon_desc_destroy( &descAB );

    RUNTIME_options_finalize( &options, chamctxt );
}

int chameleon_pzgebrd( int genD, cham_job_t jobu, cham_job_t jobvt,
                       CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                       CHAMELEON_Complex64_t *U, int LDU,
                       CHAMELEON_Complex64_t *VT, int LDVT,
                       double *E, double *S,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *subA, *subT, *subUVT, *subD;
    CHAM_desc_t descUl, descUt;
    CHAM_desc_t descVTl, descVTt;
    int M, N, NB;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return sequence->status;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    chameleon_pzgebrd_ge2gb( genD, A, T, D, sequence, request );
    chameleon_pzgebrd_gb2bd( jobu, jobvt, A, U, LDU, VT, LDVT, E, S, sequence, request );

    /* Update U and Vt according to jobu and jobvt */
    subA   = NULL;
    subT   = NULL;
    subUVT = NULL;
    subD   = NULL;
    M      = A->m;
    N      = A->n;
    NB     = A->mb;

    if ( jobu != ChamNoVec ) {
        chameleon_zlap2tile( chamctxt, &descUl, &descUt, ChamDescInout, ChamUpperLower,
                             U, NB, NB, LDU, M, M, M, sequence, request );

        if ( M < N ) {
            subA   = chameleon_desc_submatrix( A,       chameleon_min(A->mb, A->m),             0,
                                                        chameleon_max(0, A->m - A->mb),         A->n );
            subUVT = chameleon_desc_submatrix( &descUt, chameleon_min(descUt.mb, descUt.m),     0,
                                                        chameleon_max(0, descUt.m - descUt.mb), descUt.n);
            subT   = chameleon_desc_submatrix( T,       chameleon_min(T->mb, T->m),             0,
                                                        chameleon_max(0, T->m - T->mb),         T->n );
            if ( D != NULL ) {
                subD = chameleon_desc_submatrix( D,     chameleon_min(D->mb, D->m),             0,
                                                        chameleon_max(0, D->m - D->mb),         D->n );
            }
            chameleon_pzunmqr( 0, ChamLeft, ChamNoTrans, subA, subUVT, subT, subD, sequence, request );

            free( subA );
            free( subUVT );
            free( subT );
            if ( D != NULL ) {
                free( subD );
            }
        }
        else {
            chameleon_pzunmqr( 0, ChamLeft, ChamNoTrans, A, &descUt, T, D, sequence, request );
        }

        chameleon_ztile2lap( chamctxt, &descUl, &descUt, ChamDescInout, ChamUpperLower, sequence, request );
    }
    if ( jobvt != ChamNoVec ) {
        chameleon_zlap2tile( chamctxt, &descVTl, &descVTt, ChamDescInout, ChamUpperLower,
                             VT, NB, NB, LDVT, N, N, N, sequence, request );

        if ( M < N ){
            chameleon_pzunmlq( 0, ChamRight, ChamNoTrans, A, &descVTt, T, D, sequence, request );
        }
        else {
            subA   = chameleon_desc_submatrix( A,        0,         chameleon_min(A->nb, A->n),
                                                         A->m,      chameleon_max(0, A->n - A->nb) );
            subUVT = chameleon_desc_submatrix( &descVTt, 0,         chameleon_min(descVTt.nb, descVTt.n),
                                                         descVTt.m, chameleon_max(0, descVTt.n - descVTt.nb) );
            subT   = chameleon_desc_submatrix( T,        0,         chameleon_min(T->nb, T->n),
                                                         T->m,      chameleon_max(0, T->n - T->nb) );
            if ( D != NULL ) {
                subD = chameleon_desc_submatrix( D,      0,         chameleon_min(D->nb, D->n),
                                                         D->m,      chameleon_max(0, D->n - D->nb) );
            }

            chameleon_pzunmlq( 0, ChamRight, ChamNoTrans, subA, subUVT, subT, subD, sequence, request );

            free( subA );
            free( subUVT );
            free( subT );
            if ( D != NULL ) {
                free( subD );
            }
        }

        chameleon_ztile2lap( chamctxt, &descVTl, &descVTt,
                             ChamDescInout, ChamUpperLower, sequence, request );
    }
    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( T, sequence );
    if ( D != NULL ) {
        CHAMELEON_Desc_Flush( D, sequence );
    }
    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    if ( jobu != ChamNoVec ) {
        chameleon_ztile2lap_cleanup( chamctxt, &descUl,  &descUt  );
    }
    if ( jobvt != ChamNoVec ) {
        chameleon_ztile2lap_cleanup( chamctxt, &descVTl, &descVTt );
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
