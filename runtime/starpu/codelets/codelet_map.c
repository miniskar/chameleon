/**
 *
 * @file starpu/codelet_map.c
 *
 * @copyright 2018-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map StarPU codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2018-09-24
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

CHAMELEON_CL_CB(map, starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0, M*N);

#if !defined(CHAMELEON_SIMULATION)
static void cl_map_cpu_func(void *descr[], void *cl_arg)
{
    const CHAM_desc_t *desc;
    cham_uplo_t uplo;
    int m;
    int n;
    void *data;
    cham_unary_operator_t operator;
    void *op_args;

    data = (void *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &desc, &uplo, &m, &n, &operator, &op_args );
    operator( desc, uplo, m, n, data, op_args );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(map, 1, cl_map_cpu_func)

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t operator, void *op_args )
{

    struct starpu_codelet *codelet = &cl_map;
    void (*callback)(void*) = options->profiling ? cl_map_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &A,                      sizeof(CHAM_desc_t*),
        STARPU_VALUE,    &uplo,                   sizeof(cham_uplo_t),
        STARPU_VALUE,    &Am,                     sizeof(int),
        STARPU_VALUE,    &An,                     sizeof(int),
        STARPU_RW,        RTBLKADDR(A, void, Am, An),
        STARPU_VALUE,    &operator,                sizeof(cham_unary_operator_t),
        STARPU_VALUE,    &op_args,                 sizeof(void*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "map",
#endif
        0);
}
