/**
 *
 * @file chameleon.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon main header
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2018-02-08
 *
 */
#ifndef _chameleon_h_
#define _chameleon_h_

/* ****************************************************************************
 * CHAMELEON types and constants
 */
#include "chameleon/config.h"
#include "chameleon/constants.h"
#include "chameleon/types.h"
#include "chameleon/struct.h"

/* ****************************************************************************
 * CHAMELEON runtime common API
 */
#include "chameleon/runtime.h"


/* ****************************************************************************
 * CHAMELEON Simulation mode
 */
#include "chameleon/simulate.h"

/* ****************************************************************************
 * Include LibHQR for hierarchical trees QR/LQ factorizations
 */
#include "libhqr.h"

/* ****************************************************************************
 * CHAMELEON Tasks
 */
#include "chameleon/tasks.h"

#include "chameleon/chameleon_z.h"
#include "chameleon/chameleon_c.h"
#include "chameleon/chameleon_d.h"
#include "chameleon/chameleon_s.h"
#include "chameleon/chameleon_zc.h"
#include "chameleon/chameleon_ds.h"

BEGIN_C_DECLS

/* ****************************************************************************
 * CHAMELEON functionnalities
 */
int CHAMELEON_map_Tile( cham_uplo_t           uplo,
                        CHAM_desc_t          *A,
                        cham_unary_operator_t op_fct,
                        void                 *op_args );
int CHAMELEON_map_Tile_Async( cham_uplo_t           uplo,
                              CHAM_desc_t          *A,
                              cham_unary_operator_t op_fct,
                              void                 *op_args,
                              RUNTIME_sequence_t   *sequence,
                              RUNTIME_request_t    *request );

/* ****************************************************************************
 * CHAMELEON Functions
 */

/* Auxiliary */
int CHAMELEON_Version           (int *ver_major, int *ver_minor, int *ver_micro);
int CHAMELEON_My_Mpi_Rank       (void);
int __chameleon_init            (int nworkers, int ncudas);
int __chameleon_initpar         (int nworkers, int ncudas, int nthreads_per_worker);
int __chameleon_finalize        (void);
int CHAMELEON_Pause             (void);
int CHAMELEON_Resume            (void);
int CHAMELEON_Distributed_start (void);
int CHAMELEON_Distributed_stop  (void);
int CHAMELEON_Comm_size         (void);
int CHAMELEON_Comm_rank         (void);
int CHAMELEON_Lapack_to_Tile    (void *Af77, int LDA, CHAM_desc_t *A);
int CHAMELEON_Tile_to_Lapack    (CHAM_desc_t *A, void *Af77, int LDA);
int CHAMELEON_Distributed_start (void);
int CHAMELEON_Distributed_stop  (void);
int CHAMELEON_Distributed_size  (int *size);
int CHAMELEON_Distributed_rank  (int *rank);
int CHAMELEON_GetThreadNbr      (void);

/* Descriptor */
int CHAMELEON_Element_Size(int type);

int CHAMELEON_Desc_Create_User( CHAM_desc_t **desc, void *mat, cham_flttype_t dtyp, int mb, int nb, int bsiz,
                                int lm, int ln, int i, int j, int m, int n, int p, int q,
                                blkaddr_fct_t get_blkaddr, blkldd_fct_t get_blkldd, blkrankof_fct_t get_rankof );

int CHAMELEON_Desc_Create( CHAM_desc_t **desc, void *mat, cham_flttype_t dtyp,
                           int mb, int nb, int bsiz, int lm, int ln,
                           int i, int j, int m, int n, int p, int q );

int CHAMELEON_Desc_Create_OOC_User( CHAM_desc_t **desc, cham_flttype_t dtyp,
                                    int mb, int nb, int bsiz, int lm, int ln,
                                    int i, int j, int m, int n, int p, int q,
                                    blkrankof_fct_t get_rankof );
int CHAMELEON_Desc_Create_OOC( CHAM_desc_t **desc, cham_flttype_t dtyp,
                               int mb, int nb, int bsiz, int lm, int ln,
                               int i, int j, int m, int n, int p, int q );

CHAM_desc_t *CHAMELEON_Desc_Copy( CHAM_desc_t *descin, void *mat );
CHAM_desc_t *CHAMELEON_Desc_CopyOnZero( CHAM_desc_t *descin, void *mat );
CHAM_desc_t *CHAMELEON_Desc_SubMatrix( CHAM_desc_t *descA, int i, int j, int m, int n );

int CHAMELEON_Desc_Destroy( CHAM_desc_t **desc );
int CHAMELEON_Desc_Acquire( CHAM_desc_t  *desc );
int CHAMELEON_Desc_Release( CHAM_desc_t  *desc );
int CHAMELEON_Desc_Flush  ( CHAM_desc_t  *desc, RUNTIME_sequence_t *sequence );

/* Workspaces */
int CHAMELEON_Dealloc_Workspace (CHAM_desc_t **desc);

/* Options */
int  CHAMELEON_Enable  (int option);
int  CHAMELEON_Disable (int option);
int  CHAMELEON_Set     (int param, int  value);
int  CHAMELEON_Get     (int param, int *value);
int  CHAMELEON_Set_update_progress_callback(void (*p)(int, int)) ;
void CHAMELEON_user_tag_size(int, int);

/* Sequences */
int CHAMELEON_Sequence_Create  (RUNTIME_sequence_t **sequence);
int CHAMELEON_Sequence_Destroy (RUNTIME_sequence_t *sequence);
int CHAMELEON_Sequence_Wait    (RUNTIME_sequence_t *sequence);

/**
 *
 * @ingroup Control
 *
 * @brief Initialize CHAMELEON.
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
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
#if defined(CHAMELEON_SCHED_OPENMP)
#define CHAMELEON_Init(nworkers, ncudas)\
    __chameleon_init(nworkers, ncudas);\
    _Pragma("omp parallel")\
    _Pragma("omp master")\
    {
#define CHAMELEON_InitPar(nworkers, ncudas, nthreads_per_worker)\
    __chameleon_initpar(nworkers, ncudas, nthreads_per_worker);\
    _Pragma("omp parallel")\
    _Pragma("omp master")\
    {
#define CHAMELEON_Finalize()\
    }\
    __chameleon_finalize();
#else
#define CHAMELEON_Init(nworkers, ncudas)\
  __chameleon_init(nworkers, ncudas);
#define CHAMELEON_InitPar(nworkers, ncudas, nthreads_per_worker)\
  __chameleon_initpar(nworkers, ncudas, nthreads_per_worker);
#define CHAMELEON_Finalize()\
  __chameleon_finalize();
#endif

END_C_DECLS

#endif /* _chameleon_h_ */
