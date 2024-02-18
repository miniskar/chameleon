/**
 *
 * @file parsec/codelet_ipiv_to_perm.c
 *
 * @copyright 2023-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Parsec codelets to convert pivot to permutations
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2023-08-31
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks.h"
#include "coreblas.h"

static inline int
CORE_ipiv_to_perm_parsec( parsec_execution_stream_t *context,
                          parsec_task_t             *this_task )
{
    int m0, m, k;
    int *ipiv, *perm, *invp;

    parsec_dtd_unpack_args(
        this_task, &m0, &m, &k, &ipiv, &perm, &invp );

    CORE_ipiv_to_perm( m0, m, k, ipiv, perm, invp );
}

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ipiv_to_perm_parsec, options->priority, "ipiv_to_perm",
        sizeof(int),         &m0,           VALUE,
        sizeof(int),         &m,            VALUE,
        sizeof(int),         &k,            VALUE,
        PASSED_BY_REF, RUNTIME_ipiv_getaddr( ipivdesc, ipivk ), chameleon_parsec_get_arena_index_ipiv( ipivdesc ) | INPUT,
        PASSED_BY_REF, RUNTIME_perm_getaddr( ipivdesc, ipivk ), chameleon_parsec_get_arena_index_perm( ipivdesc ) | OUTPUT,
        PASSED_BY_REF, RUNTIME_invp_getaddr( ipivdesc, ipivk ), chameleon_parsec_get_arena_index_invp( ipivdesc ) | OUTPUT,
        PARSEC_DTD_ARG_END );
}
