/**
 *
 * @file parsec/codelet_zlaset2.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset2 PaRSEC codelet
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
 *  CORE_zlaset2 - Sets the elements of the matrix A to alpha.
 *  Not LAPACK compliant! Read below.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = ChamUpper: STRICT Upper part of A is set to alpha;
 *          = ChamLower: STRICT Lower part of A is set to alpha;
 *          = ChamUpperLower: ALL elements of A are set to alpha.
 *          Not LAPACK Compliant.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set to alpha accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 */
static inline int
CORE_zlaset2_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;

    parsec_dtd_unpack_args(
        this_task, &uplo, &M, &N, &alpha, &A, &LDA );

    CORE_zlaset2( uplo, M, N, alpha, A, LDA );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlaset2(const RUNTIME_option_t *options,
                        cham_uplo_t uplo, int M, int N,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int LDA)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlaset2_parsec, options->priority, "laset2",
        sizeof(int),                &uplo,      VALUE,
        sizeof(int),                       &M,         VALUE,
        sizeof(int),                       &N,         VALUE,
        sizeof(int),                &alpha,     VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | OUTPUT | AFFINITY,
        sizeof(int),                       &LDA,       VALUE,
        PARSEC_DTD_ARG_END );
}
