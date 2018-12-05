/**
 *
 * @file openmp/codelet_ztile_zero.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztile_zero StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */

#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas.h"
/**
 *
 */
void INSERT_TASK_ztile_zero( const RUNTIME_option_t *options,
                            int X1, int X2, int Y1, int Y2,
                            const CHAM_desc_t *A, int Am, int An, int lda )
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    int x, y;
    for (x = X1; x < X2; x++)
        for (y = Y1; y < Y2; y++)
            ptrA[lda*x+y] = 0.0;
}
