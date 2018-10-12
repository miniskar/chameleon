/**
 *
 * @file pzlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2014-07-21
 * @precisions normal z -> s d c
 *
 */
//ALLOC_WS :  A->mb
//ALLOC_WS :  A->nb
//WS_ADD :  A->mb + A->nb
#include "control/common.h"

#define A(m, n)    A,    (m), (n)
#define Wcol(m, n) Wcol, (m), (n)
#define Welt(m, n) Welt, (m), (n)

static inline void
chameleon_pzlange_one( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                       CHAM_desc_t *Wcol, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     * Step 1:
     *  For j in [1,P], W(i, n) = reduce( A(i+k*P, n) )
     */
    for(n = 0; n < NT; n++) {
        int mmin = ( uplo == ChamLower ) ? n                      : 0;
        int mmax = ( uplo == ChamUpper ) ? chameleon_min(n+1, MT) : MT;

        int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

        for(m = mmin; m < mmax; m++) {
            int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;
            int ldam = BLKLDD( A, m );

            if ( (n == m)  && (uplo != ChamUpperLower) ) {
                INSERT_TASK_ztrasm(
                    options,
                    ChamColumnwise, uplo, diag, tempmm, tempnn,
                    A(m, n), ldam, Wcol(m, n) );
            }
            else {
                INSERT_TASK_dzasum(
                    options,
                    ChamColumnwise, ChamUpperLower, tempmm, tempnn,
                    A(m, n), ldam, Wcol(m, n) );
            }

            if ( m >= P ) {
                INSERT_TASK_dgeadd(
                    options,
                    ChamNoTrans, tempnn, 1, A->nb,
                    1.0, Wcol(m,   n), tempnn,
                    1.0, Wcol(m%P, n), tempnn );
            }
        }

        /**
         * Step 2:
         *  For each i, W(i, n) = reduce( W(0..P-1, n) )
         */
        for(m = 1; m < P; n++) {
            INSERT_TASK_dgeadd(
                options,
                ChamNoTrans, tempnn, 1, A->nb,
                1.0, Wcol(m, n), tempnn,
                1.0, Wcol(0, n), tempnn );
        }

        INSERT_TASK_dlange(
            options,
            ChamMaxNorm, tempnn, 1, A->nb,
            Wcol(0, n), tempnn, Welt(0, n));
    }

    /**
     * Step 3:
     *  For n in 0..Q-1, W(m, n) = max( W(m, n..nt[Q] ) )
     */
    for(n = Q; n < NT; n++) {
        INSERT_TASK_dlange_max(
            options,
            Welt(0, n), Welt(0, n%Q) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(n = 1; n < Q; n++) {
        INSERT_TASK_dlange_max(
            options,
            Welt(0, n), Welt(0, 0) );
    }
}

static inline void
chameleon_pzlange_inf( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                       CHAM_desc_t *Wcol, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     * Step 1:
     *  For j in [1,Q], Wcol(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;
        int ldam = BLKLDD( A, m );

        for(n = nmin; n < nmax; n++) {
            int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

            if ( (n == m)  && (uplo != ChamUpperLower) ) {
                INSERT_TASK_ztrasm(
                    options,
                    ChamRowwise, uplo, diag, tempmm, tempnn,
                    A(m, n), ldam, Wcol(m, n) );
            }
            else {
                INSERT_TASK_dzasum(
                    options,
                    ChamRowwise, ChamUpperLower, tempmm, tempnn,
                    A(m, n), ldam, Wcol(m, n) );
            }

            if ( n >= Q ) {
                INSERT_TASK_dgeadd(
                    options,
                    ChamNoTrans, tempmm, 1, A->mb,
                    1.0, Wcol(m, n  ), tempmm,
                    1.0, Wcol(m, n%Q), tempmm );
            }
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Wcol(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_dgeadd(
                options,
                ChamNoTrans, tempmm, 1, A->mb,
                1.0, Wcol(m, n), tempmm,
                1.0, Wcol(m, 0), tempmm );
        }

        INSERT_TASK_dlange(
            options,
            ChamMaxNorm, tempmm, 1, A->nb,
            Wcol(m, 0), 1, Welt(m, 0));
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Wcol(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dlange_max(
            options,
            Welt(m, 0), Welt(m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dlange_max(
            options,
            Welt(m, 0), Welt(0, 0) );
    }
}

static inline void
chameleon_pzlange_max( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;
        int ldam = BLKLDD( A, m );

        for(n = nmin; n < nmax; n++) {
            int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

            if ( (n == m)  && (uplo != ChamUpperLower) ) {
                INSERT_TASK_zlantr(
                    options,
                    ChamMaxNorm, uplo, diag, tempmm, tempnn, A->nb,
                    A(m, n), ldam, Welt(m, n));
            }
            else {
                INSERT_TASK_zlange(
                    options,
                    ChamMaxNorm, tempmm, tempnn, A->nb,
                    A(m, n), ldam, Welt(m, n));
            }

            if ( n >= Q ) {
                INSERT_TASK_dlange_max(
                    options,
                    Welt(m, n), Welt(m, n%Q) );
            }
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Welt(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_dlange_max(
                options,
                Welt(m, n), Welt(m, 0) );
        }
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Welt(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dlange_max(
            options,
            Welt(m, 0), Welt(m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dlange_max(
            options,
            Welt(m, 0), Welt(0, 0) );
    }
}

static inline void
chameleon_pzlange_frb( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;
        int ldam = BLKLDD( A, m );

        for(n = nmin; n < nmax; n++) {
            int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

            if ( (n == m) && (uplo != ChamUpperLower) ) {
                INSERT_TASK_ztrssq(
                    options,
                    uplo, diag, tempmm, tempnn,
                    A(m, n), ldam, Welt(m, n) );
            }
            else {
                INSERT_TASK_zgessq(
                    options,
                    tempmm, tempnn,
                    A(m, n), ldam, Welt(m, n) );
            }

            if ( n >= Q ) {
                INSERT_TASK_dplssq(
                    options, Welt(m, n), Welt(m, n%Q) );
            }
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Welt(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_dplssq(
                options, Welt(m, n), Welt(m, 0) );
        }
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Welt(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dplssq(
            options, Welt(m, 0), Welt(m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dplssq(
            options, Welt(m, 0), Welt(0, 0) );
    }

    INSERT_TASK_dplssq2(
        options, Welt(0, 0) );
}

/**
 *
 */
void chameleon_pzlange_generic( cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                                CHAM_desc_t *A, double *result,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *Wcol = NULL;
    CHAM_desc_t *Welt = NULL;
    double alpha = 0.0;
    double beta  = 0.0;

    int workn, workmt, worknt;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    *result = 0.0;

    workmt = chameleon_max( A->mt, A->p );
    worknt = chameleon_max( A->nt, A->q );
    workn  = chameleon_max( A->n,  A->q );

    switch ( norm ) {
    case ChamOneNorm:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        CHAMELEON_Desc_Create( &Wcol, NULL, ChamRealDouble, 1, A->nb, A->nb,
                               workmt, worknt * A->nb, 0, 0, workmt, worknt * A->nb, A->p, A->q );

        CHAMELEON_Desc_Create( &Welt, NULL, ChamRealDouble, 1, 1, 1,
                               A->p, worknt, 0, 0, A->p, worknt, A->p, A->q );

        break;

        /*
         *  ChamInfNorm
         */
    case ChamInfNorm:
        RUNTIME_options_ws_alloc( &options, A->mb, 0 );

        CHAMELEON_Desc_Create( &Wcol, NULL, ChamRealDouble, A->mb, 1, A->mb,
                               workmt * A->mb, worknt, 0, 0, workmt * A->mb, worknt, A->p, A->q );

        CHAMELEON_Desc_Create( &Welt, NULL, ChamRealDouble, 1, 1, 1,
                               workmt, A->q, 0, 0, workmt, A->q, A->p, A->q );
        break;

        /*
         *  ChamFrobeniusNorm
         */
    case ChamFrobeniusNorm:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        alpha = 1.;
        CHAMELEON_Desc_Create( &Welt, NULL, ChamRealDouble, 2, 1, 2,
                               workmt*2, workn, 0, 0, workmt*2, workn, A->p, A->q );
        break;

        /*
         *  ChamMaxNorm
         */
    case ChamMaxNorm:
    default:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        CHAMELEON_Desc_Create( &Welt, NULL, ChamRealDouble, 1, 1, 1,
                               workmt, workn, 0, 0, workmt, workn, A->p, A->q );
    }

    /* Initialize workspaces */
    if ( (norm == ChamInfNorm) ||
         (norm == ChamOneNorm) )
    {
        /* Initialize Wcol tile */
        for(m = 0; m < Wcol->mt; m++) {
            for(n = 0; n < Wcol->nt; n++) {
                INSERT_TASK_dlaset(
                    &options,
                    ChamUpperLower, Wcol->mb, Wcol->nb,
                    alpha, beta,
                    Wcol(m,n), Wcol->mb );
            }
        }
    }
    for(m = 0; m < Welt->mt; m++) {
        for(n = 0; n < Welt->nt; n++) {
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, Welt->mb, Welt->nb,
                alpha, beta,
                Welt(m,n), Welt->mb );
        }
    }

    switch ( norm ) {
    case ChamOneNorm:
        chameleon_pzlange_one( uplo, diag, A, Wcol, Welt, &options );
        CHAMELEON_Desc_Flush( Wcol, sequence );
        break;

    case ChamInfNorm:
        chameleon_pzlange_inf( uplo, diag, A, Wcol, Welt, &options );
        CHAMELEON_Desc_Flush( Wcol, sequence );
        break;

    case ChamFrobeniusNorm:
        chameleon_pzlange_frb( uplo, diag, A, Welt, &options );
        break;

    case ChamMaxNorm:
    default:
        chameleon_pzlange_max( uplo, diag, A, Welt, &options );
    }

    /**
     * Broadcast the result
     */
    for(m = 0; m < A->p; m++) {
        for(n = 0; n < A->q; n++) {
            if ( (m != 0) && (n != 0) ) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    Welt(0,0), 1, Welt(m, n), 1);
            }
        }
    }

    CHAMELEON_Desc_Flush( Welt, sequence );
    RUNTIME_sequence_wait(chamctxt, sequence);

    *result = *(double *)Welt->get_blkaddr(Welt, A->myrank / A->q, A->myrank % A->q );

    if ( Wcol != NULL ) {
        CHAMELEON_Desc_Destroy( &Wcol );
    }
    CHAMELEON_Desc_Destroy( &Welt );

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
