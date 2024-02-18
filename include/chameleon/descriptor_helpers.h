/**
 *
 * @file descriptor_helpers.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Set of functions to help the user to declare matrix descriptors (allocation, mapping... )
 *
 * @version 1.3.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @author Lionel Eyraud-Dubois
 * @date 2023-07-05
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
int chameleon_getrankof_2d     ( const CHAM_desc_t *A, int m, int n );
int chameleon_getrankof_2d_diag( const CHAM_desc_t *A, int m, int n );

typedef struct custom_dist_s{
    int *blocks_dist;         // Matrix of size dist_m times dist_n with values from 1 to number of process MPI
    int dist_m, dist_n;       // The matrix has dist_m rows of dist_n elements
    const char* dist_file;    // Name of the file that contains the distribution
} custom_dist_t;

int chameleon_getrankof_custom_init   ( custom_dist_t **dist, const char *filename );
int chameleon_getrankof_custom_destroy( custom_dist_t **dist );
int chameleon_getrankof_custom        ( const CHAM_desc_t *A, int m, int n );

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
