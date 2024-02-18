/**
 *
 * @file coreblas_zctile.h
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800 ), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CPU kernel CHAM_tile_t interface
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2022-02-22
 * @precisions mixed zc -> ds
 *
 */
#ifndef _coreblas_zctile_h_
#define _coreblas_zctile_h_

void TCORE_clag2z( int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B );
void TCORE_zlag2c( int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B, int *info );

#endif /* _coreblas_zctile_h_ */
