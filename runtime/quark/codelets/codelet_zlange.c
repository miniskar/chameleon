/**
 *
 * @file quark/codelet_zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlange_quark(Quark *quark)
{
    double *normA;
    cham_normtype_t norm;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int LDA;
    double *work;

    quark_unpack_args_7(quark, norm, M, N, A, LDA, work, normA);
    CORE_zlange( norm, M, N, A, LDA, work, normA);
}

void INSERT_TASK_zlange(const RUNTIME_option_t *options,
                       cham_normtype_t norm, int M, int N, int NB,
                       const CHAM_desc_t *A, int Am, int An, int LDA,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANGE;
    int szeW = chameleon_max( M, N );
    QUARK_Insert_Task(
        opt->quark, CORE_zlange_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &norm,  VALUE,
        sizeof(int),                     &M,     VALUE,
        sizeof(int),                     &N,     VALUE,
        sizeof(CHAMELEON_Complex64_t)*NB*NB, RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &LDA,   VALUE,
        sizeof(double)*szeW,             NULL,   SCRATCH,
        sizeof(double),                  RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);
}

void CORE_zlange_max_quark(Quark *quark)
{
    double *A;
    double *normA;

    quark_unpack_args_2(quark, A, normA);
    if ( A[0] > *normA )
        *normA = A[0];
}

void INSERT_TASK_zlange_max(const RUNTIME_option_t *options,
                           const CHAM_desc_t *A, int Am, int An,
                           const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANGE_MAX;
    QUARK_Insert_Task(
        opt->quark, CORE_zlange_max_quark, (Quark_Task_Flags*)opt,
        sizeof(double), RTBLKADDR(A, double, Am, An), INPUT,
        sizeof(double), RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);
}

