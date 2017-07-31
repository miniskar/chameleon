/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 * @file codelet_zhessq.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @precisions normal z -> c
 *
 **/
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"

static int
CORE_zhessq_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    MORSE_enum *uplo;
    int *n;
    MORSE_Complex64_t *A;
    int *lda;
    double *SCALESUMSQ;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &uplo,
        UNPACK_VALUE, &n,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &lda,
        UNPACK_DATA,  &SCALESUMSQ );

    CORE_zhessq( *uplo, *n, A, *lda, &SCALESUMSQ[0], &SCALESUMSQ[1]);

    return 0;
}

void MORSE_TASK_zhessq( const MORSE_option_t *options,
                        MORSE_enum uplo, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zhessq_parsec, "hessq",
        sizeof(int),           &uplo,               VALUE,
        sizeof(int),           &n,                  VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),                    INPUT | REGION_FULL,
        sizeof(int),           &lda,                VALUE,
        PASSED_BY_REF,         RTBLKADDR( SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn ),    INOUT | REGION_FULL,
        0);
}
