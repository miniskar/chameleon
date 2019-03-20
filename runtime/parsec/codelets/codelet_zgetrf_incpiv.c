/**
 *
 * @file parsec/codelet_zgetrf_incpiv.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv PaRSEC codelet
 *
 * @version 0.9.2
 * @author Reazul Hoque
 * @date 2015-11-04
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 *  CORE_zgetrf_incpiv computes an LU factorization of a general M-by-N tile A
 *  using partial pivoting with row int erchanges.
 *
 *  The factorization has the form
 *
 *    A = P * L * U
 *
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factors L and U from the factorization
 *         A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile was int erchanged with row IPIV(i).
 *
 * @param[out] INFO
 *         See returned value.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if INFO = -k, the k-th argument had an illegal value
 * @retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 */
static inline int
CORE_zgetrf_incpiv_parsec( parsec_execution_stream_t *context,
                           parsec_task_t             *this_task )
{
    int m;
    int n;
    int ib;
    CHAMELEON_Complex64_t *A;
    int lda;
    int *IPIV;
    cham_bool_t *check_info;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &ib, &A, &lda, &IPIV, &check_info, &iinfo, &sequence, &request );

    CORE_zgetrf_incpiv( m, n, ib, A, lda, IPIV, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zgetrf_incpiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int lda,
                               const CHAM_desc_t *L, int Lm, int Ln, int ldl,
                               int *IPIV,
                               cham_bool_t check_info, int iinfo )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgetrf_incpiv_parsec, options->priority, "getrf_inc",
        sizeof(int),                 &m,                                VALUE,
        sizeof(int),                 &n,                                VALUE,
        sizeof(int),                 &ib,                               VALUE,
        PASSED_BY_REF,               RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int),                 &lda,                              VALUE,
        sizeof(int*),                &IPIV,                             VALUE,
        sizeof(int),                 &check_info,                       VALUE,
        sizeof(int),                 &iinfo,                            VALUE,
        sizeof(RUNTIME_sequence_t*), &(options->sequence),              VALUE,
        sizeof(RUNTIME_request_t*),  &(options->request),               VALUE,
        PARSEC_DTD_ARG_END );

    (void)L;
    (void)Lm;
    (void)Ln;
    (void)ldl;
    (void)nb;
}
