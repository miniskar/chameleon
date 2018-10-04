/**
 *
 * @file codelet_map.c
 *
 * @copyright 2018-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map Quark codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2018-09-24
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks.h"

void CORE_map_quark(Quark *quark)
{
    const CHAM_desc_t *desc;
    cham_uplo_t uplo;
    int m;
    int n;
    void *data;
    cham_unary_operator_t operator;
    void *op_args;

    quark_unpack_args_7( quark, desc, uplo, m, n, data, operator, op_args );
    operator( desc, uplo, m, n, data, op_args );
}

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t operator, void *op_args )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(
        opt->quark, CORE_map_quark, (Quark_Task_Flags*)opt,
        sizeof(CHAM_desc_t*),             &A,    VALUE,
        sizeof(cham_uplo_t),              &uplo, VALUE,
        sizeof(int),                      &Am,   VALUE,
        sizeof(int),                      &An,   VALUE,
        sizeof(char), RTBLKADDR(A, void, Am, An), INOUT,
        sizeof(cham_unary_operator_t),    &operator, VALUE,
        sizeof(void*),                    &op_args,  VALUE,
        0);
}
