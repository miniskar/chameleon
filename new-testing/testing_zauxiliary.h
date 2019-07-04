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
 * Macro for trace generation
 *
 */
#define START_TRACING()                         \
    RUNTIME_start_stats();                      \
    if(iparam[IPARAM_TRACE] == 2) {             \
    	RUNTIME_start_profiling();              \
    }                                           \
    if(iparam[IPARAM_BOUND]) {                  \
        CHAMELEON_Enable(CHAMELEON_BOUND);      \
    }

#define STOP_TRACING()                          \
    RUNTIME_stop_stats();                       \
    if(iparam[IPARAM_TRACE] == 2) {             \
    	RUNTIME_stop_profiling();               \
    }                                           \
    if(iparam[IPARAM_BOUND]) {                  \
        CHAMELEON_Disable(CHAMELEON_BOUND);     \
    }

/**
 *
 * Macro for DAG generation
 *
 */
#if 0
#define START_DAG()                   \
    if ( iparam[IPARAM_DAG] == 2 )    \
        CHAMELEON_Enable(CHAMELEON_DAG);

#define STOP_DAG()                    \
    if ( iparam[IPARAM_DAG] == 2 )    \
        CHAMELEON_Disable(CHAMELEON_DAG);
#else
#define START_DAG()  do {} while(0);
#define STOP_DAG()   do {} while(0);
#endif

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
/* #define START_TIMING()                          \ */
/*     START_DAG();                                \ */
/*     START_TRACING();                            \ */
/*     START_DISTRIBUTED();                        \ */
/*     t = -RUNTIME_get_time(); */

/* #define STOP_TIMING()                           \ */
/*     STOP_DISTRIBUTED();                         \ */
/*     t += RUNTIME_get_time();                    \ */
/*     STOP_TRACING();                             \ */
/*     STOP_DAG();                                 \ */
/*     if (iparam[IPARAM_PROFILE] == 2) {          \ */
/*         RUNTIME_kernelprofile_display();        \ */
/*         RUNTIME_schedprofile_display();         \ */
/*     }                                           \ */
/*     *t_ = t; */

#define START_TIMING( _t_ )                     \
    START_DISTRIBUTED();                        \
    (_t_) = RUNTIME_get_time();

#define STOP_TIMING( _t_ )                      \
    STOP_DISTRIBUTED();                         \
    (_t_) = RUNTIME_get_time() - (_t_);         \

#endif /* _testing_zauxiliary_h_ */
