/**
 *
 * @file testing_zauxiliary.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings header
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Cédric Castagnède
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#ifndef _testing_zauxiliary_h_
#define _testing_zauxiliary_h_

#include "testings.h"

/**
 *
 * Synchro for distributed computations
 *
 */
#if defined(CHAMELEON_USE_MPI)
#define START_DISTRIBUTED()  CHAMELEON_Distributed_start();
#define STOP_DISTRIBUTED()   CHAMELEON_Distributed_stop();
#else
#define START_DISTRIBUTED()  do {} while(0);
#define STOP_DISTRIBUTED()   do {} while(0);
#endif

/**
 *
 * General Macros for timing
 *
 */
#define START_TIMING( _t_ )                     \
    START_DISTRIBUTED();                        \
    (_t_) = RUNTIME_get_time();

#define STOP_TIMING( _t_ )                      \
    STOP_DISTRIBUTED();                         \
    (_t_) = RUNTIME_get_time() - (_t_);         \

#endif /* _testing_zauxiliary_h_ */
