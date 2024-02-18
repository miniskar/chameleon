/**
 *
 * @file pzcesca.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca parallel algorithm
 *
 * @version 1.2.0
 * @author Florent Pruvost
 * @date 2022-02-22
 * @precisions normal z -> s d c z
 *
 */
#include "control/common.h"

#define A( m, n )        A,     (m), (n)
#define W( desc, m, n ) (desc), (m), (n)

static inline void
chameleon_pzcesca_internal( int center,
                            int scale,
                            cham_store_t axis,
                            CHAM_desc_t *A,
                            CHAM_desc_t *Wgcol,
                            CHAM_desc_t *Wgrow,
                            CHAM_desc_t *Wgelt,
                            CHAM_desc_t *Wdcol,
                            CHAM_desc_t *Wdrow,
                            RUNTIME_option_t *options )
{
    int m, n;
    int MT = A->mt;
    int NT = A->nt;
    int M  = A->m;
    int N  = A->n;
    int P  = A->p;
    int Q  = A->q;

    /**
     *  1) compute sums and sum-square (scl,ssq) in each tile
     */
    for(n = 0; n < NT; n++) {
        int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;
        for(m = 0; m < MT; m++) {
            int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;
            if ( (center == 1) && ( (axis == ChamColumnwise) || (axis == ChamEltwise) ) ) {
                INSERT_TASK_zgesum(
                    options, ChamColumnwise, tempmm, tempnn,
                    A(m, n), W( Wgcol, m, n) );
            }
            if ( (center == 1) && ( (axis == ChamRowwise) || (axis == ChamEltwise) ) ) {
                INSERT_TASK_zgesum(
                    options, ChamRowwise, tempmm, tempnn,
                    A(m, n), W( Wgrow, m, n) );
            }
            if ( (scale == 1) && (axis == ChamColumnwise) ) {
                INSERT_TASK_dgessq(
                    options, ChamColumnwise, tempmm, tempnn,
                    A(m, n), W( Wdcol, m, n) );
            }
            if ( (scale == 1) && (axis == ChamRowwise) ) {
                INSERT_TASK_dgessq(
                    options, ChamRowwise, tempmm, tempnn,
                    A(m, n), W( Wdrow, m, n) );
            }
        }
    }

    for(n = 0; n < NT; n++) {
        int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

        if ( (center == 1) && ( (axis == ChamColumnwise) || (axis == ChamEltwise) ) ) {
            /**
             *  2) reduce sums of columns between tiles per processus (between lines)
             */
            for(m = P; m < MT; m++) {
                INSERT_TASK_zgeadd(
                    options, ChamNoTrans, 1, tempnn, Wgcol->nb,
                    1., W( Wgcol, m,   n ),
                    1., W( Wgcol, m%P, n ) );
            }
            /**
             *  3) reduce sums of columns between tiles of different processus on the first line of tiles
             */
            for(m = 1; m < P; m++) {
                INSERT_TASK_zgeadd(
                    options, ChamNoTrans, 1, tempnn, Wgcol->nb,
                    1., W( Wgcol, m, n ),
                    1., W( Wgcol, 0, n ) );
            }
            if ( axis == ChamEltwise ) {
                /* 4) reduce sums inside each tile of the first line of tiles for the global sum */
                INSERT_TASK_zgesum(
                    options, ChamEltwise, 1, tempnn,
                    W( Wgcol, 0, n ),
                    W( Wgelt, 0, n ) );
            }
        }

        if ( (scale == 1) && (axis == ChamColumnwise) ) {
            /**
             *  2) reduce columns (scl,ssq) tiles per processus (between lines)
             */
            for(m = P; m < MT; m++) {
                INSERT_TASK_dplssq(
                    options, ChamColumnwise, 1, tempnn,
                    W( Wdcol, m,   n ),
                    W( Wdcol, m%P, n ) );
            }
            /**
             *  3) reduce columns (scl,ssq) tiles on the first line of tiles
             */
            for(m = 1; m < P; m++) {
                INSERT_TASK_dplssq(
                    options, ChamColumnwise, 1, tempnn,
                    W( Wdcol, m, n ),
                    W( Wdcol, 0, n ) );
            }
            /* 5) deduce the sum square for each column from the pairs (scl,ssq) -> sqrt(sum) = scl*sqrt(ssq) */
            INSERT_TASK_dplssq2( options, tempnn, W( Wdcol, 0, n ) );
        }
    }

    for(m = 0; m < MT; m++) {
        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

        if ( (center == 1) && ( (axis == ChamRowwise) || (axis == ChamEltwise) ) ) {
            /**
             *  2) reduce sums of rows between tiles per processus (between columns)
             */
            for(n = Q; n < NT; n++) {
                INSERT_TASK_zgeadd(
                    options, ChamNoTrans, tempmm, 1, Wgrow->mb,
                    1., W( Wgrow, m,   n ),
                    1., W( Wgrow, m, n%Q ) );
            }
            /**
             *  3) reduce sums of rows between tiles of different processus on the first columns of tiles
             */
            for(n = 1; n < Q; n++) {
                INSERT_TASK_zgeadd(
                    options, ChamNoTrans, tempmm, 1, Wgcol->nb,
                    1., W( Wgrow, m, n ),
                    1., W( Wgrow, m, 0 ) );
            }
        }

        if ( (scale == 1) && (axis == ChamRowwise) ) {
            /**
             *  2) reduce rows (scl,ssq) tiles per processus (between lines)
             */
            for(n = Q; n < NT; n++) {
                INSERT_TASK_dplssq(
                    options, ChamRowwise, tempmm, 1,
                    W( Wdrow, m,   n ),
                    W( Wdrow, m, n%Q ) );
            }
            /**
             *  3) reduce rows (scl,ssq) tiles on the first column of tiles
             */
            for(n = 1; n < Q; n++) {
                INSERT_TASK_dplssq(
                    options, ChamRowwise, tempmm, 1,
                    W( Wdrow, m, n ),
                    W( Wdrow, m, 0 ) );
            }
            /* 5) deduce the sum square for each column from the pairs (scl,ssq) -> sqrt(sum) = scl*sqrt(ssq) */
            INSERT_TASK_dplssq2( options, tempmm, W( Wdrow, m, 0 ) );
        }
    }

    if ( (center == 1) && (axis == ChamEltwise) ) {
        /* 6) reduce global sum on each processus (between columns) */
        for(n = Q; n < NT; n++) {
            INSERT_TASK_zgeadd(
                options, ChamNoTrans, 1, 1, 1,
                1., W( Wgelt, 0, n),
                1., W( Wgelt, 0, n%Q) );
        }

        /* 7) reduce global sum on the first tile (index 0, 0) */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_zgeadd(
                options, ChamNoTrans, 1, 1, 1,
                1., W( Wgelt, 0, n),
                1., W( Wgelt, 0, 0) );
        }
    }

    /* Finally compute Centered-Scaled matrix coefficients inplace */
    for(n = 0; n < NT; n++) {
        int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

        for(m = 0; m < MT; m++) {
            int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

            INSERT_TASK_zcesca(
                options,
                center, scale, axis,
                A->m, A->n, tempmm, tempnn,
                W( Wgrow, m, 0 ),
                W( Wgcol, 0, n ),
                W( Wgelt, 0, 0 ),
                W( Wdrow, m, 0 ),
                W( Wdcol, 0, n ),
                A( m, n ) );
        }
    }
}

/**
 *
 */
void chameleon_pzcesca( struct chameleon_pzcesca_s *ws, int center, int scale, cham_store_t axis, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *Wgcol = &(ws->Wgcol);
    CHAM_desc_t *Wgrow = &(ws->Wgrow);
    CHAM_desc_t *Wgelt = &(ws->Wgelt);
    CHAM_desc_t *Wdcol = &(ws->Wdcol);
    CHAM_desc_t *Wdrow = &(ws->Wdrow);
    int m, n, tempmm, tempnn;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    /* Initialize Wgcol */
    for(m = 0; m < Wgcol->mt; m++) {
        tempmm = m == Wgcol->mt-1 ? Wgcol->m-m*Wgcol->mb : Wgcol->mb;
        for(n = 0; n < Wgcol->nt; n++) {
            tempnn = n == Wgcol->nt-1 ? Wgcol->n-n*Wgcol->nb : Wgcol->nb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempmm, tempnn,
                0., 0.,
                W( Wgcol, m, n ) );
        }
    }
    /* Initialize Wgrow */
    for(m = 0; m < Wgrow->mt; m++) {
        tempmm = m == Wgrow->mt-1 ? Wgrow->m-m*Wgrow->mb : Wgrow->mb;
        for(n = 0; n < Wgrow->nt; n++) {
            tempnn = n == Wgrow->nt-1 ? Wgrow->n-n*Wgrow->nb : Wgrow->nb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempmm, tempnn,
                0., 0.,
                W( Wgrow, m, n ) );
        }
    }
    /* Initialize Wgelt */
    for(m = 0; m < Wgelt->mt; m++) {
        tempmm = m == Wgelt->mt-1 ? Wgelt->m-m*Wgelt->mb : Wgelt->mb;
        for(n = 0; n < Wgelt->nt; n++) {
            tempnn = n == Wgelt->nt-1 ? Wgelt->n-n*Wgelt->nb : Wgelt->nb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempmm, tempnn,
                0., 0.,
                W( Wgelt, m, n ) );
        }
    }
    /* Initialize Wdcol */
    for(m = 0; m < Wdcol->mt; m++) {
        tempmm = m == Wdcol->mt-1 ? Wdcol->m-m*Wdcol->mb : Wdcol->mb;
        for(n = 0; n < Wdcol->nt; n++) {
            tempnn = n == Wdcol->nt-1 ? Wdcol->n-n*Wdcol->nb : Wdcol->nb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempmm, tempnn,
                -1., -1.,
                W( Wdcol, m, n ) );
        }
    }
    /* Initialize Wdrow */
    for(m = 0; m < Wdrow->mt; m++) {
        tempmm = m == Wdrow->mt-1 ? Wdrow->m-m*Wdrow->mb : Wdrow->mb;
        for(n = 0; n < Wdrow->nt; n++) {
            tempnn = n == Wdrow->nt-1 ? Wdrow->n-n*Wdrow->nb : Wdrow->nb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempmm, tempnn,
                -1., -1.,
                W( Wdrow, m, n ) );
        }
    }

    chameleon_pzcesca_internal( center, scale, axis, A, Wgcol, Wgrow, Wgelt, Wdcol, Wdrow, &options );

    CHAMELEON_Desc_Flush( Wgcol, sequence );
    CHAMELEON_Desc_Flush( Wgrow, sequence );
    CHAMELEON_Desc_Flush( Wgelt, sequence );
    CHAMELEON_Desc_Flush( Wdcol, sequence );
    CHAMELEON_Desc_Flush( Wdrow, sequence );
    RUNTIME_options_finalize(&options, chamctxt);
}
