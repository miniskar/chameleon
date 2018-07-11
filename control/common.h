/**
 *
 * @file common.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon common header file
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
/**
 *  CHAMELEON facilities of interest to both CHAMELEON core developer
 *  and also of interest to CHAMELEON community contributor.
 */
#ifndef _CHAMELEON_COMMON_H_
#define _CHAMELEON_COMMON_H_


#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif

/**
 * Implementation headers
 */
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas.h>
#include <cublas_v2.h>
#else
#include <cublas.h>
#endif
#endif

#if defined(CHAMELEON_USE_OPENCL) && !defined(CHAMELEON_SIMULATION)
#include <OpenCL/cl.h>
#endif

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/**
 *  Line to avoid conflict with other linear algebra libraries, because, we
 *  don't know why but lapacke provide a wrong interface of lapack in fortran
 */
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif

/**
 *  Chameleon header files
 */
#include "chameleon.h"

#include "control/global.h"
#include "control/auxiliary.h"
#include "control/context.h"
#include "control/descriptor.h"
#include "control/async.h"

/**
 *  Global shortcuts
 */
#define CHAMELEON_RANK        morse_rank(morse)
#define CHAMELEON_SIZE        morse->world_size
#define CHAMELEON_GRPSIZE     morse->group_size
#define CHAMELEON_NB          morse->nb
#define CHAMELEON_IB          morse->ib
#define CHAMELEON_SCHEDULING  morse->scheduling
#define CHAMELEON_RHBLK       morse->rhblock
#define CHAMELEON_TRANSLATION morse->translation
#define CHAMELEON_PARALLEL    morse->parallel_enabled
#define CHAMELEON_PROFILING   morse->profiling_enabled
#if defined(CHAMELEON_USE_MPI)
#define CHAMELEON_MPI_RANK    morse->my_mpi_rank
#define CHAMELEON_MPI_SIZE    morse->mpi_comm_size
#endif

/**
 *  IPT internal define
 */
#define ChamIPT_NoDep   0
#define ChamIPT_Panel   1
#define ChamIPT_All     2

/**
 *  Global array of LAPACK constants
 */
extern char *morse_lapack_constants[];
#define morse_lapack_const(morse_const) morse_lapack_constants[morse_const][0]

#ifdef __cplusplus
extern "C" {
#endif

#include "control/compute_s.h"
#include "control/compute_d.h"
#include "control/compute_c.h"
#include "control/compute_z.h"

/*
void morse_pdlag2s(CHAM_context_t *morse);
void morse_pzlag2c(CHAM_context_t *morse);
void morse_pslag2d(CHAM_context_t *morse);
void morse_pclag2z(CHAM_context_t *morse);
*/

#ifdef __cplusplus
}
#endif

#endif
