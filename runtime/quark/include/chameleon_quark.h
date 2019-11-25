/**
 *
 * @file quark/chameleon_quark.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark runtime main header
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2017-07-31
 *
 */
#ifndef _chameleon_quark_h_
#define _chameleon_quark_h_

#include "control/common.h"

#include <quark.h>
#if defined(CHAMELEON_RUNTIME_SYNC)
#define QUARk_Insert_Task QUARK_Execute_Task
#endif
#include "coreblas.h"
#include "core_blas_dag.h"

typedef struct quark_option_s {
    Quark_Task_Flags flags;
    Quark *quark;
} quark_option_t;

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (type*)RUNTIME_data_getaddr( desc, m, n ) )

#define RUNTIME_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION

#endif /* _chameleon_quark_h_ */
