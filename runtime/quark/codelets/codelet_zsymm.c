/**
 *
 * @file quark/codelet_zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm Quark codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zsymm_quark(Quark *quark)
{
    cham_side_t side;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileC;

    quark_unpack_args_9(quark, side, uplo, M, N, alpha, tileA, tileB, beta, tileC);
    TCORE_zsymm(side, uplo,
        M, N,
        alpha, tileA,
        tileB,
        beta, tileC);
}

void INSERT_TASK_zsymm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                      const CHAM_desc_t *B, int Bm, int Bn,
                      CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn)
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                                    beta, C, Cm, Cn );
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int accessC = ( beta == 0. ) ? OUTPUT : INOUT;

    DAG_CORE_SYMM;
    QUARK_Insert_Task(opt->quark, CORE_zsymm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &side,    VALUE,
        sizeof(int),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha,   VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),               INPUT,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),               INPUT,
        sizeof(CHAMELEON_Complex64_t),         &beta,    VALUE,
        sizeof(void*), RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),               accessC,
        0);
}
