/**
 *
 * @file cblas_wrapper.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cblas header wrapper
 *
 * @version 1.2.0
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2022-02-22
 *
 */
#ifndef _cblas_wrapper_h_
#define _cblas_wrapper_h_

/**
 *  CBLAS requires for scalar arguments to be passed
 *        by address rather than by value
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif
#include "coreblas/cblas.h"

/**
 * CBlas enum
 */
#define CBLAS_ORDER     enum CBLAS_ORDER
#define CBLAS_TRANSPOSE enum CBLAS_TRANSPOSE
#define CBLAS_UPLO      enum CBLAS_UPLO
#define CBLAS_DIAG      enum CBLAS_DIAG
#define CBLAS_SIDE      enum CBLAS_SIDE

#endif /* _cblas_wrapper_h_ */
