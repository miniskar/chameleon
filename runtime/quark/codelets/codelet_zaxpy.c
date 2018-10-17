/**
 *
 * @file quark/codelet_zaxpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zaxpy Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zaxpy_quark(Quark *quark)
{
    int M;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int incA;
    CHAMELEON_Complex64_t *B;
    int incB;

    quark_unpack_args_6(quark, M, alpha, A, incA, B, incB);
    CORE_zaxpy(M, alpha, A, incA, B, incB);
}

void INSERT_TASK_zaxpy(const RUNTIME_option_t *options,
                      int M, CHAMELEON_Complex64_t *alpha,
                      const CHAM_desc_t *A, int Am, int An, int incA,
                      const CHAM_desc_t *B, int Bm, int Bn, int incB)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_AXPY;
    QUARK_Insert_Task(opt->quark, CORE_zaxpy_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &M,         VALUE,
        sizeof(CHAMELEON_Complex64_t),          alpha,      VALUE,
        sizeof(CHAMELEON_Complex64_t)*M,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(int),                        &incA,      VALUE,
        sizeof(CHAMELEON_Complex64_t)*M,        RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn), INOUT,
        sizeof(int),                        &incB,      VALUE,
        0);
}
