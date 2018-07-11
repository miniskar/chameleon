/**
 *
 * @file pzlansy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlansy parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
//ALLOC_WS :  A->mb
#include <stdlib.h>
#include <math.h>
//WS_ADD :  A->mb
#include "control/common.h"

#define A(m, n) A, m, n
#define VECNORMS_STEP1(m, n) VECNORMS_STEP1, m, n
#define VECNORMS_STEP2(m, n) VECNORMS_STEP2, m, n
#define RESULT(m, n) RESULT, m, n
/**
 *
 */
/**
 *
 */
void morse_pzlansy(cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A, double *result,
                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_desc_t *VECNORMS_STEP1 = NULL;
    CHAM_desc_t *VECNORMS_STEP2 = NULL;
    CHAM_desc_t *RESULT         = NULL;
    CHAM_context_t *morse;
    RUNTIME_option_t options;

    int workm, workn;
    int tempkm, tempkn;
    int ldam;
    int m, n;
    /* int part_p, part_q; */

    /* part_p = A->myrank / A->q; */
    /* part_q = A->myrank % A->q; */

    morse = morse_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    *result = 0.0;
    switch ( norm ) {
    /*
     *  ChamOneNorm / ChamInfNorm
     */
    case ChamOneNorm:
    case ChamInfNorm:
        /* Init workspace handle for the call to zlange */
        RUNTIME_options_ws_alloc( &options, A->mb, 0 );

        workm = A->m;
        workn = chameleon_max( A->nt, A->q );
        CHAMELEON_Desc_Create(&(VECNORMS_STEP1), NULL, ChamRealDouble, A->mb, 1, A->mb,
                          workm, workn, 0, 0, workm, workn, A->p, A->q);

        CHAMELEON_Desc_Create(&(VECNORMS_STEP2), NULL, ChamRealDouble, A->mb, 1, A->mb,
                          workm, 1, 0, 0, workm, 1, A->p, A->q);

        CHAMELEON_Desc_Create(&(RESULT), NULL, ChamRealDouble, 1, 1, 1,
                          1, 1, 0, 0, 1, 1, 1, 1);

        /* Zeroes my intermediate vectors */
        for(m = 0; m < A->mt; m++) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            for(n = 0; n < workn; n++) {
                INSERT_TASK_dlaset(
                    &options,
                    ChamUpperLower, tempkm, 1,
                    0., 0.,
                    VECNORMS_STEP1(m, n), 1);
            }
        }

        for(m = 0; m < A->mt; m++) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);

            /* compute sums of absolute values on diagonal tile m */
            INSERT_TASK_dzasum(
                &options,
                ChamRowwise, uplo, tempkm, tempkm,
                A(m, m), ldam, VECNORMS_STEP1(m, m));

            /*
             *  ChamLower
             */
            if (uplo == ChamLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n < m; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* compute sums of absolute values on rows of tile m */
                    INSERT_TASK_dzasum(
                        &options,
                        ChamRowwise, ChamUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    INSERT_TASK_dzasum(
                        &options,
                        ChamColumnwise, ChamUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(n, m));
                }
            }
            /*
             *  ChamUpper
             */
            else {
                //                for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
                //                    n < A->mt; n+=A->q) {
                for(n = m+1; n < A->mt; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* compute sums of absolute values on rows of tile m */
                    INSERT_TASK_dzasum(
                        &options,
                        ChamRowwise, ChamUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    INSERT_TASK_dzasum(
                        &options,
                        ChamColumnwise, ChamUpperLower, tempkm, tempkn,
                        A(m, n), ldam, VECNORMS_STEP1(n, m));
                }
            }
        }

        /* compute vector sum between tiles in rows */
        for(m = 0; m < A->mt; m++) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempkm, 1,
                0., 0.,
                VECNORMS_STEP2(m, 0), 1);
            for(n = 0; n < A->nt; n++) {
                INSERT_TASK_dgeadd(
                    &options,
                    ChamNoTrans, tempkm, 1, A->mb,
                    1.0, VECNORMS_STEP1(m, n), tempkm,
                    1.0, VECNORMS_STEP2(m, 0), tempkm);
            }
            /*
             * Compute max norm of each segment of the final vector in the
             * previous workspace
             */
            INSERT_TASK_dlange(
                &options,
                ChamMaxNorm, tempkm, 1, A->nb,
                VECNORMS_STEP2(m, 0), tempkm,
                VECNORMS_STEP1(m, 0));
        }

        /* Initialize RESULT array */
        INSERT_TASK_dlaset(
            &options,
            ChamUpperLower, 1, 1,
            0., 0.,
            RESULT(0,0), 1);

        /* compute max norm between tiles in the column */
        if (A->myrank % A->q == 0) {
            for(m = 0; m < A->mt; m++) {
                INSERT_TASK_dlange_max(
                    &options,
                    VECNORMS_STEP1(m, 0),
                    RESULT(0,0));
            }
        }

        /* Scatter norm over processus */
        for(m = 0; m < A->p; m++) {
            for(n = 0; n < A->q; n++) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    RESULT(0,0), 1,
                    VECNORMS_STEP1(m, n), 1 );
            }
        }
        CHAMELEON_Desc_Flush( VECNORMS_STEP2, sequence );
        CHAMELEON_Desc_Flush( VECNORMS_STEP1, sequence );
        CHAMELEON_Desc_Flush( RESULT, sequence );
        RUNTIME_sequence_wait(morse, sequence);
        CHAMELEON_Desc_Destroy( &(VECNORMS_STEP2) );
        break;
    /*
     *  ChamFrobeniusNorm
     */
    case ChamFrobeniusNorm:
        workm = chameleon_max( A->mt, A->p );
        workn = chameleon_max( A->nt, A->q );

        CHAMELEON_Desc_Create(&(VECNORMS_STEP1), NULL, ChamRealDouble, 1, 2, 2,
                          workm, 2*workn, 0, 0, workm, 2*workn, A->p, A->q);
        CHAMELEON_Desc_Create(&(RESULT), NULL, ChamRealDouble, 1, 2, 2,
                          1, 2, 0, 0, 1, 2, 1, 1);

        /* Compute local norm to each tile */
        for(m = (A->myrank / A->q); m < A->mt; m+=A->p) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);

            /* Zeroes my intermediate vector */
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, 1, 2,
                1., 0.,
                VECNORMS_STEP1(m, m), 1);
            /* compute norm on diagonal tile m */
            INSERT_TASK_zsyssq(
                &options,
                uplo, tempkm,
                A(m, m), ldam,
                VECNORMS_STEP1(m, m));

            /*
             *  ChamLower
             */
            if (uplo == ChamLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n < m; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* Zeroes my intermediate vector */
                    INSERT_TASK_dlaset(
                        &options,
                        ChamUpperLower, 1, 2,
                        1., 0.,
                        VECNORMS_STEP1(m, n), 1);
                    /* compute norm on the lower part */
                    INSERT_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    INSERT_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                }
            }
            /*
             *  ChamUpper
             */
            else {
                //                for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
                //                    n < A->mt; n+=A->q) {
                for(n = m+1; n < A->mt; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    /* Zeroes my intermediate vector */
                    INSERT_TASK_dlaset(
                        &options,
                        ChamUpperLower, 1, 2,
                        1., 0.,
                        VECNORMS_STEP1(m, n), 1);
                    /* compute norm on the lower part */
                    INSERT_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                    /* same operation on the symmetric part */
                    INSERT_TASK_zgessq(
                        &options,
                        tempkm, tempkn,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                }
            }
        }

        /* Zeroes my intermediate vector */
        INSERT_TASK_dlaset(
            &options,
            ChamUpperLower, 1, 2,
            1., 0.,
            RESULT(0,0), 1);

        /* Compute accumulation of scl and ssq */
        for(m = (A->myrank / A->q); m < A->mt; m+=A->p) {
            /*
             *  ChamLower
             */
            if (uplo == ChamLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n <= m; n++) {
                    INSERT_TASK_dplssq(
                        &options,
                        VECNORMS_STEP1(m, n),
                        RESULT(0,0));
                }
            }
            /*
             *  ChamUpper
             */
            else {
                //                for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
                //                    n < A->mt; n+=A->q) {
                for(n = m; n < A->mt; n++) {
                    INSERT_TASK_dplssq(
                        &options,
                        VECNORMS_STEP1(m, n),
                        RESULT(0,0));
                }
            }
        }

        /* Compute scl * sqrt(ssq) */
        INSERT_TASK_dplssq2(
            &options,
            RESULT(0,0));

        /* Copy max norm in tiles to dispatch on every nodes */
        for(m = 0; m < A->p; m++) {
            for(n = 0; n < A->q; n++) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    RESULT(0,0), 1,
                    VECNORMS_STEP1(m, n), 1 );
            }
        }

        CHAMELEON_Desc_Flush( VECNORMS_STEP1, sequence );
        CHAMELEON_Desc_Flush( RESULT, sequence );
        RUNTIME_sequence_wait(morse, sequence);
        break;

    /*
     *  ChamMaxNorm
     */
    case ChamMaxNorm:
    default:
        /* Init workspace handle for the call to zlange but unused */
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        workm = chameleon_max( A->mt, A->p );
        workn = chameleon_max( A->nt, A->q );

        CHAMELEON_Desc_Create(&(VECNORMS_STEP1), NULL, ChamRealDouble, 1, 1, 1,
                          workm, workn, 0, 0, workm, workn, A->p, A->q);
        CHAMELEON_Desc_Create(&(RESULT), NULL, ChamRealDouble, 1, 1, 1,
                          1, 1, 0, 0, 1, 1, 1, 1);

        /* Compute local maximum to each tile */
        for(m = 0; m < A->mt; m++) {
            tempkm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);

            INSERT_TASK_zlansy(
                &options,
                ChamMaxNorm, uplo, tempkm, A->nb,
                A(m, m), ldam,
                VECNORMS_STEP1(m, m));

            /*
             *  ChamLower
             */
            if (uplo == ChamLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n < m; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    INSERT_TASK_zlange(
                        &options,
                        ChamMaxNorm, tempkm, tempkn, A->nb,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                }
            }
            /*
             *  ChamUpper
             */
            else {
                //for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
                //  n < A->mt; n+=A->q) {
                for(n = m+1; n < A->mt; n++) {
                    tempkn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    INSERT_TASK_zlange(
                        &options,
                        ChamMaxNorm, tempkm, tempkn, A->nb,
                        A(m, n), ldam,
                        VECNORMS_STEP1(m, n));
                }
            }
        }

        /* Zeroes RESULT array */
        INSERT_TASK_dlaset(
            &options,
            ChamUpperLower, 1, 1,
            0., 0.,
            RESULT(0,0), 1);

        /* Compute max norm between tiles */
        for(m = 0; m < A->mt; m++) {
            /*
             *  ChamLower
             */
            if (uplo == ChamLower) {
                //for(n = A->myrank % A->q; n < m; n+=A->q) {
                for(n = 0; n <= m; n++) {
                    INSERT_TASK_dlange_max(
                        &options,
                        VECNORMS_STEP1(m, n),
                        RESULT(0,0));
                }
            }
            /*
             *  ChamUpper
             */
            else {
                //for(n = ( part_q > part_p ?  (m/part_p)*part_p + part_q : (m/part_p)*part_p + part_q + A->q );
                //  n < A->mt; n+=A->q) {
                for(n = m; n < A->mt; n++) {
                    INSERT_TASK_dlange_max(
                        &options,
                        VECNORMS_STEP1(m, n),
                        RESULT(0,0));
                }
            }
        }

        /* Copy max norm in tiles to dispatch on every nodes */
        for(m = 0; m < A->p; m++) {
            for(n = 0; n < A->q; n++) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    RESULT(0,0), 1,
                    VECNORMS_STEP1(m, n), 1 );
            }
        }

        CHAMELEON_Desc_Flush( VECNORMS_STEP1, sequence );
        CHAMELEON_Desc_Flush( RESULT, sequence );
        RUNTIME_sequence_wait(morse, sequence);
    }

    *result = *(double *)VECNORMS_STEP1->get_blkaddr(VECNORMS_STEP1, A->myrank / A->q, A->myrank % A->q );

    CHAMELEON_Desc_Destroy( &(VECNORMS_STEP1) );
    CHAMELEON_Desc_Destroy( &(RESULT) );
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
