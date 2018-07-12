/**
 *
 * @file codelet_zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zgetrf_quark(Quark *quark)
{
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    int *IPIV;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    cham_bool_t check_info;
    int iinfo;
    int info;

    quark_unpack_args_9(quark, m, n, A, lda, IPIV, sequence, request, check_info, iinfo);
    CORE_zgetrf( m, n, A, lda, IPIV, &info );
    if ( (info != CHAMELEON_SUCCESS) && check_info ) {
        RUNTIME_sequence_flush( (CHAM_context_t*)quark, sequence, request, iinfo+info );
    }
}

void INSERT_TASK_zgetrf(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       int *IPIV,
                       cham_bool_t check_info, int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GETRF;
    QUARK_Insert_Task(opt->quark, CORE_zgetrf_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(RUNTIME_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(RUNTIME_request_t*),            &(options->request),       VALUE,
        sizeof(cham_bool_t),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}
