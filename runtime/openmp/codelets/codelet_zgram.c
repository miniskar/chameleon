/**
 *
 * @file openmp/codelet_zgram.c
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram OpenMP codelet
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2019-04-10
 * @precisions normal z -> s d c z
 *
 */

#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din, int lddi,
                        const CHAM_desc_t *Dj, int Djm, int Djn, int lddj,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An, int lda)
{
    double *ptrDi = RTBLKADDR(Di, double, Dim, Din);
    double *ptrDj = RTBLKADDR(Dj, double, Djm, Djn);
    double *ptrD  = RTBLKADDR(D,  double, Dm, Dn);
    double *ptrA  = RTBLKADDR(A,  double, Am, An);

#pragma omp task firstprivate(uplo, m, n, mt, nt, ptrDi, lddi, ptrDj, lddj, ptrD, ptrA, lda) depend(in:ptrDi[0], ptrDj[0], ptrD[0]) depend(inout:ptrA[0])
    CORE_zgram( uplo,
                m, n, mt, nt,
                ptrDi, lddi,
                ptrDj, lddj,
                ptrD,
                ptrA, lda);
}
