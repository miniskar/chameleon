/**
 *
 * @file quark/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy Quark codelet
 *
 * @version 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

static inline void CORE_zlacpy_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int M;
    int N;
    int displA;
    CHAM_tile_t *tileA;
    CHAMELEON_Complex64_t *A;
    int displB;
    CHAM_tile_t *tileB;
    CHAMELEON_Complex64_t *B;

    quark_unpack_args_7(quark, uplo, M, N, displA, tileA, displB, tileB);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

    A = tileA->mat;
    B = tileB->mat;
    CORE_zlacpy( uplo, M, N, A + displA, tileA->ld, B + displB, tileB->ld );
}

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int nb,
                          int displA, const CHAM_desc_t *A, int Am, int An,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LACPY;
    QUARK_Insert_Task(opt->quark, CORE_zlacpy_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo,   VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(int),                     &displA, VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INPUT,
        sizeof(int),                     &displB, VALUE,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             OUTPUT,
        0);
    (void)nb;
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    INSERT_TASK_zlacpyx( options, uplo, m, n, nb,
                         0, A, Am, An,
                         0, B, Bm, Bn );
}
