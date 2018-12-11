/**
 *
 * @file openmp/chameleon_openmp.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP runtime main header
 *
 * @version 1.0.0
 * @author Philippe Virouleau
 * @date 2018-06-21
 *
 */
#ifndef _CHAMELEON_OPENMP_H_
#define _CHAMELEON_OPENMP_H_

#include "coreblas.h"

#include "control/common.h"
#include <omp.h>

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (type*)RUNTIME_data_getaddr( desc, m, n ) )


#endif /* _CHAMELEON_OPENMP_H_ */
