/**
 *
 * @file auxiliary.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#ifndef _CHAMELEON_AUXILIARY_H_
#define _CHAMELEON_AUXILIARY_H_

#include "chameleon/struct.h"
#include "chameleon/tasks.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
void morse_warning      (const char *func_name, const char* msg_text);
void morse_error        (const char *func_name, const char* msg_text);
void morse_fatal_error  (const char *func_name, const char* msg_text);
int  morse_rank         (CHAM_context_t *morse);
int  morse_tune         (cham_tasktype_t func, int M, int N, int NRHS);

/**
 *  API routines
 */
int  CHAMELEON_Version      (int *ver_major, int *ver_minor, int *ver_micro);
int  CHAMELEON_Element_Size (int type);
int  CHAMELEON_My_Mpi_Rank  (void);

#ifdef __cplusplus
}
#endif

#endif
