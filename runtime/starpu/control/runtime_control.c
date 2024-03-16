/**
 *
 * @file starpu/runtime_control.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU control routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Swartvagher
 * @author Samuel Thibault
 * @author Matthieu Kuhn
 * @author Loris Lucido
 * @author Terry Cojean
 * @date 2024-03-16
 *
 */
#include "chameleon_starpu.h"
#include <stdio.h>
#include <stdlib.h>
#if defined(STARPU_USE_FXT)
#include <starpu_fxt.h>
#endif

static int starpu_initialized = 0;

#if defined(STARPU_HAVE_HWLOC) && defined(HAVE_STARPU_PARALLEL_WORKER)
void chameleon_starpu_parallel_worker_init( starpu_sched_opt_t *sched_opt )
{
    char *env_pw_level = chameleon_getenv( "CHAMELEON_PARALLEL_WORKER_LEVEL" );

    if (env_pw_level != NULL) {
        struct starpu_parallel_worker_config *pw_config = NULL;
        hwloc_obj_type_t pw_level;
        int  pw_level_number = 1;
        char level[256];

        int argc  = strchr( env_pw_level, ':') == NULL ? 1 : 2;
        int match = sscanf( env_pw_level, "%255[^:]:%d", level, &pw_level_number );

#if !defined(CHAMELEON_KERNELS_MT)
        chameleon_warning("chameleon_starpu_parallel_worker_init()", "CHAMELEON has been compiled with multi-threaded kernels disabled (-DCHAMELEON_KERNELS_MT=OFF). This won't break the execution, but you may not obtain the performance gain expected. It is recommended to recompile with -DCHAMELEON_KERNELS_MT=ON.\n");
#endif

        if ( (match != argc) ||
             ((match == 2) && (pw_level_number < 0) ) )
        {
            fprintf( stderr, "error CHAMELEON_PARALLEL_WORKER_LEVEL \"%s\"  does not match the format level[:number] where number > 0.\n", env_pw_level );
            exit(1);
        }

        if ( hwloc_type_sscanf( level, &pw_level, NULL, 0 ) == -1 )
        {
            fprintf( stderr, "error CHAMELEON_PARALLEL_WORKER_LEVEL \"%s\"  does not match an hwloc level.\n", level );
            exit(1);
        }

        pw_config = starpu_parallel_worker_init( pw_level,
                                                 STARPU_PARALLEL_WORKER_NB, pw_level_number,
                                                 STARPU_PARALLEL_WORKER_TYPE, STARPU_PARALLEL_WORKER_GNU_OPENMP_MKL,
                                                 0 );

        if ( pw_config == NULL )
        {
            fprintf( stderr, "error CHAMELEON_PARALLEL_WORKER_LEVEL : cannot create a parallel worker at %s level.\n", level );
            exit(1);
        }

        if ( chameleon_env_on_off( "CHAMELEON_PARALLEL_WORKER_SHOW", CHAMELEON_FALSE ) == CHAMELEON_TRUE ) {
            starpu_parallel_worker_print( pw_config );
        }

        sched_opt->pw_config = pw_config;
    }

    chameleon_cleanenv( env_pw_level );
}

void chameleon_starpu_parallel_worker_fini( starpu_sched_opt_t *sched_opt )
{
    if ( sched_opt->pw_config != NULL ) {
        starpu_parallel_worker_shutdown( sched_opt->pw_config );
        sched_opt->pw_config = NULL;
    }
}
#else
#define chameleon_starpu_parallel_worker_init(sched_opt) do { (void) sched_opt; } while(0)
#define chameleon_starpu_parallel_worker_fini(sched_opt) do { (void) sched_opt; } while(0)
#endif

/**
 *
 */
static int chameleon_starpu_init( MPI_Comm comm, struct starpu_conf *conf )
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
        rc = starpu_mpi_init_conf(NULL, NULL, !flag, comm, conf);
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
    starpu_sched_opt_t *sched_opt = (starpu_sched_opt_t*)(chamctxt->schedopt);
    struct starpu_conf *conf = &sched_opt->starpu_conf;
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
        hres = chameleon_starpu_init( chamctxt->comm, conf );

        chamctxt->nworkers = ncpus;
        chamctxt->nthreads_per_worker = nthreads_per_worker;
    }
    else {
        int worker;

        for (worker = 0; worker < ncpus; worker++)
            conf->workers_bindid[worker] = (worker+1)*nthreads_per_worker - 1;

        for (worker = 0; worker < ncpus; worker++)
            conf->workers_bindid[worker + ncudas] = worker*nthreads_per_worker;

        conf->use_explicit_workers_bindid = 1;

        hres = chameleon_starpu_init( chamctxt->comm, conf );

        chamctxt->nworkers = ncpus;
        chamctxt->nthreads_per_worker = nthreads_per_worker;
    }

    if ( hres != CHAMELEON_SUCCESS ) {
        return hres;
    }

    starpu_initialized = 1;

#ifdef HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS
    starpu_malloc_on_node_set_default_flags( STARPU_MAIN_RAM,
                                             STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT
#ifdef STARPU_MALLOC_SIMULATION_FOLDED
                                             | STARPU_MALLOC_SIMULATION_FOLDED
#endif
                                             );
#endif

#if !defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA)
    starpu_cublas_init();
#elif defined(CHAMELEON_USE_HIP)
    starpu_hipblas_init();
#endif
#endif

    starpu_cham_tile_interface_init();
    cppi_interface_init();

    chameleon_starpu_parallel_worker_init( sched_opt );
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

    starpu_sched_opt_t *sched_opt = (starpu_sched_opt_t*)(chamctxt->schedopt);
    chameleon_starpu_parallel_worker_fini( sched_opt );

    starpu_cham_tile_interface_fini();

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    starpu_cublas_shutdown();
#endif
#if defined(CHAMELEON_USE_HIP) && !defined(CHAMELEON_SIMULATION)
    starpu_hipblas_shutdown();
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
    starpu_mpi_wait_for_all( chamctxt->comm );
    starpu_mpi_barrier( chamctxt->comm );
#  else
    starpu_task_wait_for_all();
    starpu_mpi_barrier( chamctxt->comm );
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
    starpu_mpi_comm_rank( chamctxt->comm, &rank );
#  else
    MPI_Comm_rank( chamctxt->comm, &rank );
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
    starpu_mpi_comm_size( chamctxt->comm, &size );
#  else
    MPI_Comm_size( chamctxt->comm, &size );
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
