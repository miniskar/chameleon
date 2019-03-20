/**
 *
 * @file quark/codelet_zlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum Quark codelet
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

void CORE_zlauum_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int N;
    CHAMELEON_Complex64_t *A;
    int LDA;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    CORE_zlauum(uplo, N, A, LDA);
}

void INSERT_TASK_zlauum(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(opt->quark, CORE_zlauum_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}
