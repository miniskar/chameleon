/**
 *
 * @file async.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon asynchronous management header
 *
 * @version 0.9.2
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @date 2014-11-16
 *
 */
#ifndef _chameleon_async_h_
#define _chameleon_async_h_

#include "chameleon/struct.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
int chameleon_request_fail     (RUNTIME_sequence_t *sequence, RUNTIME_request_t *request, int error);
int chameleon_sequence_create  (CHAM_context_t *CHAMELEON, RUNTIME_sequence_t **sequence);
int chameleon_sequence_destroy (CHAM_context_t *CHAMELEON, RUNTIME_sequence_t *sequence);
int chameleon_sequence_wait    (CHAM_context_t *CHAMELEON, RUNTIME_sequence_t *sequence);

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_async_h_ */
