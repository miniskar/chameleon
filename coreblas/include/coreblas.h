/**
 *
 * @file coreblas.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU kernels main header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2020-03-03
 *
 */
#ifndef _coreblas_h_
#define _coreblas_h_

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/**
 *  CBLAS requires for scalar arguments to be passed
 *        by address rather than by value
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif
#include "coreblas/cblas.h"

/**
 * CHAMELEON types and constants
 */
#include "chameleon/config.h"
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "chameleon/constants.h"

/**
 * CORE BLAS headers
 */
BEGIN_C_DECLS

#include "coreblas/coreblas_z.h"
#include "coreblas/coreblas_d.h"
#include "coreblas/coreblas_c.h"
#include "coreblas/coreblas_s.h"
#include "coreblas/coreblas_zc.h"
#include "coreblas/coreblas_ds.h"

END_C_DECLS

/**
 * Coreblas Error
 */
#define coreblas_error(k, str) do {                                     \
        fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str) ;  \
        assert(0);                                                      \
    } while(0)

/**
 * CBlas enum
 */
#define CBLAS_TRANSPOSE enum CBLAS_TRANSPOSE
#define CBLAS_UPLO      enum CBLAS_UPLO
#define CBLAS_DIAG      enum CBLAS_DIAG
#define CBLAS_SIDE      enum CBLAS_SIDE

/**
 *  LAPACK Constants
 */
BEGIN_C_DECLS

extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

void set_coreblas_gemm3m_enabled( int v );
int  get_coreblas_gemm3m_enabled( void );


#if defined(CHAMELEON_KERNELS_TRACE)

void __coreblas_kernel_trace( const char *func, ... );
#define coreblas_kernel_trace( ... )                                    \
    do {                                                                \
        __coreblas_kernel_trace( __chameleon_func__, ##__VA_ARGS__, NULL );  \
    } while(0)

#else

#define coreblas_kernel_trace( ... ) do {} while(0)

#endif

END_C_DECLS

#endif /* _coreblas_h_ */
