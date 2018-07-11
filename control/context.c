/**
 *
 * @file context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context management routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Options
 * @brief Group routines exposed to users to handle options
 *
 */

#include <stdlib.h>
#if defined( _WIN32 ) || defined( _WIN64 )
#include "control/morsewinthread.h"
#else
#include <pthread.h>
#endif

#include "control/common.h"
#include "control/auxiliary.h"
#include "control/context.h"
#include "chameleon/runtime.h"

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas.h"
#endif

/**
 *  Global data
 */
/* master threads context lookup table */
static CHAM_context_t *morse_ctxt = NULL;

/**
 *  Create new context
 */
CHAM_context_t *morse_context_create()
{
    CHAM_context_t *morse;

    if ( morse_ctxt != NULL ) {
        morse_error("morse_context_create", "a context is already existing\n");
        return NULL;
    }

    morse = (CHAM_context_t*)malloc(sizeof(CHAM_context_t));
    if (morse == NULL) {
        morse_error("morse_context_create", "malloc() failed");
        return NULL;
    }

    /* These initializations are just in case the user
       disables autotuning and does not set nb and ib */
    morse->nb                 = 128;
    morse->ib                 = 32;
    morse->rhblock            = 4;

    morse->nworkers           = 1;
    morse->ncudas             = 0;
    morse->nthreads_per_worker= 1;

    morse->warnings_enabled     = CHAMELEON_TRUE;
    morse->autotuning_enabled   = CHAMELEON_TRUE;
    morse->parallel_enabled     = CHAMELEON_FALSE;
    morse->profiling_enabled    = CHAMELEON_FALSE;
    morse->progress_enabled     = CHAMELEON_FALSE;

    morse->householder        = ChamFlatHouseholder;
    morse->translation        = ChamOutOfPlace;


    /* Initialize scheduler */
    RUNTIME_context_create(morse);

    morse_ctxt = morse;
    return morse;
}


/**
 *  Return context for a thread
 */
CHAM_context_t *morse_context_self()
{
    return morse_ctxt;
}

/**
 *  Clean the context
 */
int morse_context_destroy(){

    RUNTIME_context_destroy(morse_ctxt);
    free(morse_ctxt);
    morse_ctxt = NULL;

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Enable - Enable CHAMELEON feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be enabled:
 *          @arg CHAMELEON_WARNINGS   printing of warning messages,
 *          @arg CHAMELEON_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg CHAMELEON_PROFILING_MODE  activate profiling of kernels
 *          @arg CHAMELEON_PROGRESS  activate progress indicator
 *          @arg CHAMELEON_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Enable(int option)
{
    CHAM_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Enable", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    switch (option)
    {
        case CHAMELEON_WARNINGS:
            morse->warnings_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_AUTOTUNING:
            morse->autotuning_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_PROFILING_MODE:
            morse->profiling_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_PROGRESS:
            morse->progress_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(1);
#else
            morse_error("CHAMELEON_Enable", "cannot enable GEMM3M (not available in cblas)");
#endif
            break;
        /* case CHAMELEON_PARALLEL: */
        /*     morse->parallel_enabled = CHAMELEON_TRUE; */
        /*     break; */
        default:
            morse_error("CHAMELEON_Enable", "illegal parameter value");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
        case CHAMELEON_BOUND:
            break;
    }

    /* Enable at the lower level if required */
    RUNTIME_enable( option );

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Disable - Disable CHAMELEON feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be disabled:
 *          @arg CHAMELEON_WARNINGS   printing of warning messages,
 *          @arg CHAMELEON_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg CHAMELEON_PROFILING_MODE  deactivate profiling of kernels
 *          @arg CHAMELEON_PROGRESS  deactivate progress indicator
 *          @arg CHAMELEON_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Disable(int option)
{
    CHAM_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Disable", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch ( option )
    {
        case CHAMELEON_WARNINGS:
            morse->warnings_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_AUTOTUNING:
            morse->autotuning_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_PROFILING_MODE:
            morse->profiling_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_PROGRESS:
            morse->progress_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(0);
#endif
            break;
        case CHAMELEON_PARALLEL_MODE:
            morse->parallel_enabled = CHAMELEON_FALSE;
            break;
        default:
            morse_error("CHAMELEON_Disable", "illegal parameter value");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    /* Disable at the lower level if required */
    RUNTIME_disable( option );

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Set - Set CHAMELEON parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg CHAMELEON_TILE_SIZE:        size matrix tile,
 *          @arg CHAMELEON_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Set(int param, int value)
{
    CHAM_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Set", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case CHAMELEON_TILE_SIZE:
            if (value <= 0) {
                morse_error("CHAMELEON_Set", "negative tile size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            morse->nb = value;
            if ( morse->autotuning_enabled ) {
                morse->autotuning_enabled = CHAMELEON_FALSE;
                morse_warning("CHAMELEON_Set", "autotuning has been automatically disable\n");
            }
            /* Limit ib to nb */
            morse->ib = chameleon_min( morse->nb, morse->ib );
            break;
        case CHAMELEON_INNER_BLOCK_SIZE:
            if (value <= 0) {
                morse_error("CHAMELEON_Set", "negative inner block size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            if (value > morse->nb) {
                morse_error("CHAMELEON_Set", "inner block larger than tile");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            /* if (morse->nb % value != 0) { */
            /*     morse_error("CHAMELEON_Set", "inner block does not divide tile"); */
            /*     return CHAMELEON_ERR_ILLEGAL_VALUE; */
            /* } */
            morse->ib = value;

            if ( morse->autotuning_enabled ) {
                morse->autotuning_enabled = CHAMELEON_FALSE;
                morse_warning("CHAMELEON_Set", "autotuning has been automatically disable\n");
            }
            break;
        case CHAMELEON_HOUSEHOLDER_MODE:
            if (value != ChamFlatHouseholder && value != ChamTreeHouseholder) {
                morse_error("CHAMELEON_Set", "illegal value of CHAMELEON_HOUSEHOLDER_MODE");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            morse->householder = value;
            break;
        case CHAMELEON_HOUSEHOLDER_SIZE:
            if (value <= 0) {
                morse_error("CHAMELEON_Set", "negative householder size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            morse->rhblock = value;
            break;
        case CHAMELEON_TRANSLATION_MODE:
            if (value != ChamInPlace && value != ChamOutOfPlace) {
                morse_error("CHAMELEON_Set", "illegal value of CHAMELEON_TRANSLATION_MODE");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            morse->translation = value;
            break;
        default:
            morse_error("CHAMELEON_Set", "unknown parameter");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Get - Get value of CHAMELEON parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg CHAMELEON_TILE_SIZE:        size matrix tile,
 *          @arg CHAMELEON_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[out] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Get(int param, int *value)
{
    CHAM_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Get", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case CHAMELEON_TILE_SIZE:
            *value = morse->nb;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_INNER_BLOCK_SIZE:
            *value = morse->ib;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_HOUSEHOLDER_MODE:
            *value = morse->householder;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_HOUSEHOLDER_SIZE:
            *value = morse->rhblock;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_TRANSLATION_MODE:
            *value = morse->translation;
            return CHAMELEON_SUCCESS;
        default:
            morse_error("CHAMELEON_Get", "unknown parameter");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    return CHAMELEON_SUCCESS;
}
