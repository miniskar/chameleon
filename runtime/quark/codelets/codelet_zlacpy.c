/**
 *
 * @file quark/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy Quark codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline void CORE_zlacpy_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int M;
    int N;
    int displA;
    CHAMELEON_Complex64_t *A;
    int LDA;
    int displB;
    CHAMELEON_Complex64_t *B;
    int LDB;

    quark_unpack_args_9(quark, uplo, M, N, displA, A, LDA, displB, B, LDB);
    CORE_zlacpy(uplo, M, N, A + displA, LDA, B + displB, LDB);
}

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int nb,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LACPY;
    QUARK_Insert_Task(opt->quark, CORE_zlacpy_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo,   VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(int),                     &displA, VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INPUT,
        sizeof(int),                     &lda,    VALUE,
        sizeof(int),                     &displB, VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             OUTPUT,
        sizeof(int),                     &ldb,    VALUE,
        0);
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    INSERT_TASK_zlacpyx( options, uplo, m, n, nb,
                         0, A, Am, An, lda,
                         0, B, Bm, Bn, ldb );
}
