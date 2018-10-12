/**
 *
 * @file chameleon.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon main header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2012-09-15
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

/* ****************************************************************************
 * CHAMELEON functionnalities
 */
int CHAMELEON_map_Tile( cham_uplo_t           uplo,
                        CHAM_desc_t          *A,
                        cham_unary_operator_t operator,
                        void                 *op_args );
int CHAMELEON_map_Tile_Async( cham_uplo_t           uplo,
                              CHAM_desc_t          *A,
                              cham_unary_operator_t operator,
                              void                 *op_args,
                              RUNTIME_sequence_t   *sequence,
                              RUNTIME_request_t    *request );

#include "chameleon/chameleon_z.h"
#include "chameleon/chameleon_c.h"
#include "chameleon/chameleon_d.h"
#include "chameleon/chameleon_s.h"
#include "chameleon/chameleon_zc.h"
#include "chameleon/chameleon_ds.h"

/* ****************************************************************************
 * CHAMELEON Functions
 */
BEGIN_C_DECLS

/* Auxiliary */
int CHAMELEON_Version           (int *ver_major, int *ver_minor, int *ver_micro);
int CHAMELEON_My_Mpi_Rank       (void);
int CHAMELEON_Init              (int nworkers, int ncudas);
int CHAMELEON_InitPar           (int nworkers, int ncudas, int nthreads_per_worker);
int CHAMELEON_Finalize          (void);
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
int CHAMELEON_Desc_Create  (CHAM_desc_t **desc, void *mat, cham_flttype_t dtyp,
                        int mb, int nb, int bsiz, int lm, int ln,
                        int i, int j, int m, int n, int p, int q);
int CHAMELEON_Desc_Create_User(CHAM_desc_t **desc, void *mat, cham_flttype_t dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n, int p, int q,
                           void* (*get_blkaddr)( const CHAM_desc_t*, int, int ),
                           int (*get_blkldd)( const CHAM_desc_t*, int ),
                           int (*get_rankof)( const CHAM_desc_t*, int, int ));
int CHAMELEON_Desc_Create_OOC(CHAM_desc_t **desc, cham_flttype_t dtyp,
                          int mb, int nb, int bsiz, int lm, int ln,
                          int i, int j, int m, int n, int p, int q);
int CHAMELEON_Desc_Create_OOC_User(CHAM_desc_t **desc, cham_flttype_t dtyp,
                               int mb, int nb, int bsiz, int lm, int ln,
                               int i, int j, int m, int n, int p, int q,
                               int (*get_rankof)( const CHAM_desc_t*, int, int ));
int CHAMELEON_Desc_Destroy (CHAM_desc_t **desc);
int CHAMELEON_Desc_Acquire (CHAM_desc_t  *desc);
int CHAMELEON_Desc_Release (CHAM_desc_t  *desc);
int CHAMELEON_Desc_Flush   (CHAM_desc_t  *desc, RUNTIME_sequence_t *sequence);
void CHAMELEON_user_tag_size(int, int) ;

/* Workspaces */
int CHAMELEON_Dealloc_Workspace (CHAM_desc_t **desc);

/* Options */
int CHAMELEON_Enable  (int option);
int CHAMELEON_Disable (int option);
int CHAMELEON_Set     (int param, int  value);
int CHAMELEON_Get     (int param, int *value);
int CHAMELEON_Set_update_progress_callback(void (*p)(int, int)) ;

/* Sequences */
int CHAMELEON_Sequence_Create  (RUNTIME_sequence_t **sequence);
int CHAMELEON_Sequence_Destroy (RUNTIME_sequence_t *sequence);
int CHAMELEON_Sequence_Wait    (RUNTIME_sequence_t *sequence);

END_C_DECLS

#endif /* _chameleon_h_ */
