/**
 *
 * @file descriptor_helpers.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Set of functions to help the user to declare matrix descriptors (allocation, mapping... )
 *
 * @version 1.2.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @date 2020-03-03
 *
 * @addtogroup chameleon_descriptors
 * @{
 *   @brief Set of predefined functions to produce standard matrix mappings
 *
 *   This module provides the set of functions to produce standard mapping of
 *   the matrix tiles among processors, as well as associated function to: get a
 *   tile address within a specific storage, get a tile leading dimension, ...
 *
 */
#ifndef _chameleon_descriptor_helpers_h_
#define _chameleon_descriptor_helpers_h_

#include <chameleon/struct.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @name Mapping functions
 * @{
 */
int chameleon_getrankof_2d      ( const CHAM_desc_t *A, int m, int n );
int chameleon_getrankof_2d_diag ( const CHAM_desc_t *A, int m, int n );

/**
 * @}
 * @name Block address functions
 * @{
 */
void* chameleon_getaddr_cm  ( const CHAM_desc_t *A, int m, int n );
void* chameleon_getaddr_ccrb( const CHAM_desc_t *A, int m, int n );
void* chameleon_getaddr_null( const CHAM_desc_t *A, int m, int n );
void* chameleon_getaddr_diag( const CHAM_desc_t *A, int m, int n );

/**
 * @}
 * @name Leading dimension functions
 * @{
 */
int chameleon_getblkldd_cm  ( const CHAM_desc_t *A, int m );
int chameleon_getblkldd_ccrb( const CHAM_desc_t *A, int m );
/**
 * @}
 */

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_descriptor_helpers_h_ */
