/**
 *
 * @file control.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon control routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Control
 * @brief Group routines exposed to users to control CHAMELEON state
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "control/auxiliary.h"
#include "control/common.h"
#include "chameleon/runtime.h"

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Init - Initialize CHAMELEON.
 *
 ******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use.
 *
 * @param[in] gpus
 *          Number of cuda devices to use.
 *
 ******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Init(int cores, int gpus)
{
    return CHAMELEON_InitPar(cores, gpus, -1);
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_InitPar - Initialize CHAMELEON.
 *
 ******************************************************************************
 *
 * @param[in] ncpus
 *          Number of cores to use.
 *
 * @param[in] ncudas
 *          Number of cuda devices to use.
 *
 * @param[in] nthreads_per_worker
 *          Number of threads per worker (cpu, cuda device).
 *
 ******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_InitPar(int ncpus, int ncudas, int nthreads_per_worker)
{
    CHAM_context_t *morse;

    /* Create context and insert in the context map */
    morse = morse_context_create();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_Init", "morse_context_create() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

#if defined(CHAMELEON_USE_MPI)
#  if defined(CHAMELEON_SIMULATION)
    /* Assuming that we don't initialize MPI ourself (which SMPI doesn't support anyway) */
    morse->mpi_outer_init = 1;
#  else
    {
      int flag = 0, provided = 0;
      MPI_Initialized( &flag );
      morse->mpi_outer_init = flag;
      if ( !flag ) {
          MPI_Init_thread( NULL, NULL, MPI_THREAD_MULTIPLE, &provided );
      }
    }
#  endif
#endif

    RUNTIME_init( morse, ncpus, ncudas, nthreads_per_worker );

#if defined(CHAMELEON_USE_MPI)
    morse->my_mpi_rank   = RUNTIME_comm_rank( morse );
    morse->mpi_comm_size = RUNTIME_comm_size( morse );
#endif

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Finalize - Finalize CHAMELEON.
 *
 ******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Finalize(void)
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Finalize()", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    RUNTIME_flush();
#  if !defined(CHAMELEON_SIMULATION)
    RUNTIME_barrier(morse);
#  endif
    RUNTIME_finalize( morse );

#if defined(CHAMELEON_USE_MPI)
    if (!morse->mpi_outer_init)
        MPI_Finalize();
#endif

    morse_context_destroy();
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Pause - Suspend CHAMELEON runtime to poll for new tasks.
 *
 ******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Pause(void)
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Pause()", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    RUNTIME_pause(morse);
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Resume - Symmetrical call to CHAMELEON_Pause,
 *  used to resume the workers polling for new tasks.
 *
 ******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Resume(void)
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Resume()", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    RUNTIME_resume(morse);
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Distributed_start - Prepare the distributed processes for computation
 *
 ******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Distributed_start(void)
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Finalize()", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    RUNTIME_barrier (morse);
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Distributed_stop - Clean the distributed processes after computation
 *
 ******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Distributed_stop(void)
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Finalize()", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    RUNTIME_barrier (morse);
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Comm_size - Return the size of the distributed computation
 *
 ******************************************************************************
 *
 * @retval The size of the distributed computation
 * @retval -1 if context not initialized
 *
 */
int CHAMELEON_Comm_size()
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Comm_size()", "CHAMELEON not initialized");
        return -1;
    }

    return RUNTIME_comm_size( morse );
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Comm_rank - Return the rank of the distributed computation
 *
 ******************************************************************************
 *
 * @retval The rank of the distributed computation
 * @retval -1 if context not initialized
 *
 */
int CHAMELEON_Comm_rank()
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_Comm_rank()", "CHAMELEON not initialized");
        return -1;
    }

    return RUNTIME_comm_rank( morse );
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_GetThreadNbr - Return the number of CPU workers initialized by the
 *  runtime
 *
 ******************************************************************************
 *
 * @return
 *          \retval The number of CPU workers started
 *
 */
int CHAMELEON_GetThreadNbr( )
{
    CHAM_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("CHAMELEON_GetThreadNbr()", "CHAMELEON not initialized");
        return -1;
    }

    return RUNTIME_thread_size( morse );
}
