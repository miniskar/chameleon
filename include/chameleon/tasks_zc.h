/**
 *
 * @file chameleon_tasks_zc.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon complex mixed precision elementary tasks header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#ifndef _chameleon_tasks_zc_h_
#define _chameleon_tasks_zc_h_

/**
 *  Declarations of QUARK wrappers (called by CHAMELEON) - alphabetical order
 */
void INSERT_TASK_clag2z( const RUNTIME_option_t *options,
                        int m, int n, int nb,
                        CHAM_desc_t *A, int Am, int An, int lda,
                        CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_zlag2c( const RUNTIME_option_t *options,
                        int m, int n, int nb,
                        CHAM_desc_t *A, int Am, int An, int lda,
                        CHAM_desc_t *B, int Bm, int Bn, int ldb );

#endif
