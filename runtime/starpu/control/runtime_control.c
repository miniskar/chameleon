/**
 *
 * @file starpu/runtime_control.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU control routines
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Swartvagher
 * @author Samuel Thibault
 * @author Matthieu Kuhn
 * @date 2022-02-22
 *
 */
#include "chameleon_starpu.h"
#include <stdio.h>
#include <stdlib.h>
#if defined(STARPU_USE_FXT)
#include <starpu_fxt.h>
#endif

static int starpu_initialized = 0;

/**
 *
 */
static int chameleon_starpu_init( starpu_conf_t *conf )
{
    int hres = CHAMELEON_SUCCESS;
    int rc;

#if defined(STARPU_USE_FXT)
    starpu_fxt_autostart_profiling(0);
#endif

#if defined(CHAMELEON_USE_MPI)

    {
        int flag = 0;
#  if !defined(CHAMELEON_SIMULATION)
        MPI_Initialized( &flag );
#  endif

#  if defined(HAVE_STARPU_MPI_INIT_CONF)
        rc = starpu_mpi_init_conf(NULL, NULL, !flag, MPI_COMM_WORLD, conf);
#  else
        rc = starpu_init(conf);
        if (rc < 0) {
            return CHAMELEON_ERR_NOT_INITIALIZED;
        }
        rc = starpu_mpi_init(NULL, NULL, !flag);
#  endif
    }
#else

    rc = starpu_init(conf);

#endif
    if ( rc == -ENODEV ) {
        hres = CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Stop profiling as it seems that autostart is not sufficient */
#if defined(STARPU_USE_FXT)
    starpu_fxt_stop_profiling();
#endif

    return hres;
}

int RUNTIME_init( CHAM_context_t *chamctxt,
                  int ncpus,
                  int ncudas,
                  int nthreads_per_worker )
{
    starpu_conf_t *conf = (starpu_conf_t*)(chamctxt->schedopt);
    int hres = CHAMELEON_ERR_NOT_INITIALIZED;

    /* StarPU was already initialized by an external library */
    if (conf == NULL) {
        return 0;
    }

    if (ncpus != -1) {
        conf->ncpus = ncpus;
    }
    conf->ncuda = ncudas;
    conf->nopencl = 0;

    /* By default, use the dmdas strategy */
    if (!getenv("STARPU_SCHED")) {
        if (conf->ncuda > 0) {
            conf->sched_policy_name = "dmdas";
        }
        else {
            /**
             * Set scheduling to "ws"/"lws" if no cuda devices used because it
             * behaves better on homogneneous architectures. If the user wants
             * to use another scheduling strategy, he can set STARPU_SCHED
             * env. var. to whatever he wants
             */
#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 2))
            conf->sched_policy_name = "lws";
#else
            conf->sched_policy_name = "ws";
#endif
        }
    }

    if ((ncpus == -1)||(nthreads_per_worker == -1))
    {
        chamctxt->parallel_enabled = CHAMELEON_FALSE;

        hres = chameleon_starpu_init( conf );

        chamctxt->nworkers = ncpus;
        chamctxt->nthreads_per_worker = nthreads_per_worker;
    }
    else {
        int worker;

        chamctxt->parallel_enabled = CHAMELEON_TRUE;

        for (worker = 0; worker < ncpus; worker++)
            conf->workers_bindid[worker] = (worker+1)*nthreads_per_worker - 1;

        for (worker = 0; worker < ncpus; worker++)
            conf->workers_bindid[worker + ncudas] = worker*nthreads_per_worker;

        conf->use_explicit_workers_bindid = 1;

        hres = chameleon_starpu_init( conf );

        chamctxt->nworkers = ncpus;
        chamctxt->nthreads_per_worker = nthreads_per_worker;
    }

    if ( hres != CHAMELEON_SUCCESS ) {
        return hres;
    }

    starpu_initialized = 1;

#ifdef HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS
    starpu_malloc_on_node_set_default_flags(STARPU_MAIN_RAM, STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT
#ifdef STARPU_MALLOC_SIMULATION_FOLDED
            | STARPU_MALLOC_SIMULATION_FOLDED
#endif
            );
#endif

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    starpu_cublas_init();
#endif

    starpu_cham_tile_interface_init();

    return hres;
}

/**
 *
 */
void RUNTIME_finalize( CHAM_context_t *chamctxt )
{
    /* StarPU was already initialized by an external library or was not successfully initialized: */
    if ( (chamctxt->schedopt == NULL) || !starpu_initialized ) {
        return;
    }

    starpu_cham_tile_interface_fini();

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    starpu_cublas_shutdown();
#endif

#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_shutdown();
#endif

#if !defined(CHAMELEON_USE_MPI) || !defined(HAVE_STARPU_MPI_INIT_CONF)
    starpu_shutdown();
#endif
    return;
}

/**
 *  To suspend the processing of new tasks by workers
 */
void RUNTIME_pause( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    starpu_pause();
    return;
}

/**
 *  This is the symmetrical call to RUNTIME_pause,
 *  used to resume the workers polling for new tasks.
 */
void RUNTIME_resume( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    starpu_resume();
    return;
}

/**
 *  Busy-waiting barrier
 */
void RUNTIME_barrier( CHAM_context_t *chamctxt )
{
    (void)chamctxt;

#if defined(CHAMELEON_USE_MPI)
#  if defined(HAVE_STARPU_MPI_WAIT_FOR_ALL)
    starpu_mpi_wait_for_all(MPI_COMM_WORLD);
    starpu_mpi_barrier(MPI_COMM_WORLD);
#  else
    starpu_task_wait_for_all();
    starpu_mpi_barrier(MPI_COMM_WORLD);
#  endif
#else
    starpu_task_wait_for_all();
#endif
}

// Defined in control/auxilliary.c
extern void (*update_progress_callback)(int, int);

// no progress indicator for algorithms faster than 'PROGRESS_MINIMUM_DURATION' seconds
#define PROGRESS_MINIMUM_DURATION 10

/**
 *  Display a progress information when executing the tasks
 */
void RUNTIME_progress( CHAM_context_t *chamctxt )
{
    int tasksLeft, current, timer = 0;
    int max;

#if defined(CHAMELEON_USE_MPI)
    if ( RUNTIME_comm_rank( chamctxt ) != 0 ) {
        return;
    }
#endif

    max = starpu_task_nsubmitted();
    if ( max == 0 ) {
        return;
    }

    //  update_progress_callback(0, max);
    while ((tasksLeft = starpu_task_nsubmitted()) > 0) {
        current = max - tasksLeft;
        if (timer > PROGRESS_MINIMUM_DURATION) {
            update_progress_callback(current, max);
        }
        sleep(1);
        timer++;
    }
    if (timer > PROGRESS_MINIMUM_DURATION) {
        update_progress_callback(max, max);
    }

    (void)chamctxt;
    return;
}

/**
 * Thread rank.
 */
int RUNTIME_thread_rank( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return starpu_worker_get_id();
}

/**
 * Number of threads.
 */
int RUNTIME_thread_size( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return starpu_worker_get_count_by_type( STARPU_CPU_WORKER );
}

/**
 *  The process rank
 */
int RUNTIME_comm_rank( CHAM_context_t *chamctxt )
{
    int rank = 0;

#if defined(CHAMELEON_USE_MPI)
#  if defined(HAVE_STARPU_MPI_COMM_RANK)
    starpu_mpi_comm_rank( MPI_COMM_WORLD, &rank );
#  else
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#  endif
#endif

    (void)chamctxt;
    return rank;
}

/**
 *  This returns the size of the distributed computation
 */
int RUNTIME_comm_size( CHAM_context_t *chamctxt )
{
    int size;
#if defined(CHAMELEON_USE_MPI)
#  if defined(HAVE_STARPU_MPI_COMM_RANK)
    starpu_mpi_comm_size( MPI_COMM_WORLD, &size );
#  else
    MPI_Comm_size( MPI_COMM_WORLD, &size );
#  endif
#else
    size = 1;
#endif

    (void)chamctxt;
    return size;
}

void RUNTIME_set_minmax_submitted_tasks( int min, int max ){
#if defined(HAVE_STARPU_SET_LIMIT_SUBMITTED_TASKS)
    starpu_set_limit_min_submitted_tasks( min );
    starpu_set_limit_max_submitted_tasks( max );
#else
    fprintf( stderr,
             "RUNTIME_set_minmax_submitted_tasks: StarPU version does not support dynamic limit setting.\n"
             "Please use setting through environment variables:\n"
             "    export STARPU_LIMIT_MIN_SUBMITTED_TASKS=%d\n"
             "    export STARPU_LIMIT_MAX_SUBMITTED_TASKS=%d\n", min, max );
#endif
}
