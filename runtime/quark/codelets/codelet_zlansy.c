/**
 *
 * @file quark/codelet_zlansy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlansy Quark codelet
 *
 * @version 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zlansy_quark(Quark *quark)
{
    CHAM_tile_t *tileNorm;
    cham_normtype_t norm;
    cham_uplo_t uplo;
    int N;
    CHAM_tile_t *tileA;
    double *work;

    quark_unpack_args_6(quark, norm, uplo, N, tileA, work, tileNorm );
    TCORE_zlansy( norm, uplo, N, tileA, work, tileNorm->mat );
}

void INSERT_TASK_zlansy(const RUNTIME_option_t *options,
                       cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANSY;
    int szeW = chameleon_max( 1, N );
    QUARK_Insert_Task(
        opt->quark, CORE_zlansy_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &norm,  VALUE,
        sizeof(int),              &uplo,  VALUE,
        sizeof(int),                     &N,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(double)*szeW,             NULL,   SCRATCH,
        sizeof(void*), RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);
}
