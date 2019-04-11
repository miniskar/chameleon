/**
 *
 * @file starpu/runtime_codelet_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU CHAMELEON_Complex64_t codelets header
 *
 * @version 0.9.2
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2015-05-22
 * @precisions normal z -> c d s
 *
 */
#ifndef _runtime_codelet_z_h_
#define _runtime_codelet_z_h_

#include <stdio.h>
#include "runtime_codelets.h"

#include "chameleon/tasks_z.h"
#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/coreblas_z.h"
#if defined(CHAMELEON_USE_CUDA)
#include "cudablas.h"
#endif
#endif

/*
 * BLAS 1 functions
 */
CODELETS_HEADER(zaxpy)

/*
 * BLAS 3 functions
 */
CODELETS_HEADER(zgemm)
CODELETS_HEADER(zhemm)
CODELETS_HEADER(zher2k)
CODELETS_HEADER(zherk)
CODELETS_HEADER(zsymm)
CODELETS_HEADER(zsyr2k)
CODELETS_HEADER(zsyrk)
CODELETS_HEADER(ztrmm)
CODELETS_HEADER(ztrsm)

/*
 * LAPACK functions
 */
CODELETS_HEADER(zgelqt)
CODELETS_HEADER(zgeqrt)
CODELETS_HEADER(zgessm)
CODELETS_HEADER(zgessq)
CODELETS_HEADER(zgetrf)
CODELETS_HEADER(zgetrf_incpiv)
CODELETS_HEADER(zgetrf_nopiv)
CODELETS_HEADER(zherfb)
CODELETS_HEADER(zlauum)
CODELETS_HEADER(zpotrf)
CODELETS_HEADER(zssssm)
CODELETS_HEADER(zsyssq)
CODELETS_HEADER(ztrasm)
CODELETS_HEADER(ztrssq)
CODELETS_HEADER(ztrtri)
CODELETS_HEADER(ztplqt)
CODELETS_HEADER(ztpqrt)
CODELETS_HEADER(ztpmlqt)
CODELETS_HEADER(ztpmqrt)
CODELETS_HEADER(ztsmlq_hetra1)
CODELETS_HEADER(ztsmqr_hetra1)
CODELETS_HEADER(ztstrf)
CODELETS_HEADER(zunmlq)
CODELETS_HEADER(zunmqr)

/*
 * Auxiliary functions
 */
CODELETS_HEADER(zgeadd)
CODELETS_HEADER(zhe2ge)
CODELETS_HEADER(zlascal)
CODELETS_HEADER(ztradd)
CODELETS_HEADER(zlacpy)
CODELETS_HEADER(zlange)
CODELETS_HEADER(zlange_max)
CODELETS_HEADER(zlansy)
CODELETS_HEADER(zlantr)
CODELETS_HEADER(zlaset)
CODELETS_HEADER(zlaset2)
CODELETS_HEADER(zlatro)
CODELETS_HEADER(zplssq)
CODELETS_HEADER(zplssq2)

/*
 * MIXED PRECISION functions
 */
CODELETS_HEADER(zlag2c)

/*
 * DZ functions
 */
CODELETS_HEADER(dzasum)

/*
 * CPU only functions
 */
CODELETS_HEADER(zplrnt)
CODELETS_HEADER(zbuild)

#if defined(PRECISION_z) || defined(PRECISION_c)
CODELETS_HEADER(zhessq)
CODELETS_HEADER(zlanhe)
CODELETS_HEADER(zplghe)
CODELETS_HEADER(zsytrf_nopiv)
#endif
CODELETS_HEADER(zplgsy)

#endif /* _runtime_codelet_z_h_ */
