/**
 *
 * @file starpu/codelet_ipiv_to_perm.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to convert pivot to permutations
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2023-08-31
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelets.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ipiv_to_perm_cpu_func( void *descr[], void *cl_arg )
{
    int m0, m, k;
    int *ipiv, *perm, *invp;

    starpu_codelet_unpack_args( cl_arg, &m0, &m, &k );

    ipiv = (int*)STARPU_VECTOR_GET_PTR(descr[0]);
    perm = (int*)STARPU_VECTOR_GET_PTR(descr[1]);
    invp = (int*)STARPU_VECTOR_GET_PTR(descr[2]);

    CORE_ipiv_to_perm( m0, m, k, ipiv, perm, invp );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
* Codelet definition
*/
static struct starpu_codelet cl_ipiv_to_perm = {
    .where        = STARPU_CPU,
#if defined(CHAMELEON_SIMULATION)
    .cpu_funcs[0] = (starpu_cpu_func_t)1,
#else
    .cpu_funcs[0] = cl_ipiv_to_perm_cpu_func,
#endif
    .nbuffers     = 3,
    .model        = NULL,
    .name         = "ipiv_to_perm"
};

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    struct starpu_codelet *codelet = &cl_ipiv_to_perm;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m0,  sizeof(int),
        STARPU_VALUE,             &m,   sizeof(int),
        STARPU_VALUE,             &k,   sizeof(int),
        STARPU_R,                 RUNTIME_ipiv_getaddr( ipivdesc, ipivk ),
        STARPU_W,                 RUNTIME_perm_getaddr( ipivdesc, ipivk ),
        STARPU_W,                 RUNTIME_invp_getaddr( ipivdesc, ipivk ),
        STARPU_PRIORITY,          options->priority,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}
