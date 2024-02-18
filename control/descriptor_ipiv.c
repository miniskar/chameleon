/**
 *
 * @file descriptor_ipiv.c
 *
 * @copyright 2022-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2023-08-22
 *
 ***
 *
 * @defgroup Descriptor
 * @brief Group descriptor routines exposed to users to manipulate IPIV data structures
 *
 */
#define _GNU_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "control/common.h"
#include "control/descriptor.h"
#include "chameleon/runtime.h"

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to create tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to initialize.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated ipiv descriptor must be generated.
 *
 * @param[in] data
 *          The pointer to the original vector where to store the pivot values.
 *
 ******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on success, CHAMELEON_ERR_NOT_INITIALIZED otherwise.
 *
 */
int chameleon_ipiv_init( CHAM_ipiv_t *ipiv, const CHAM_desc_t *desc, void *data )
{
    CHAM_context_t *chamctxt;
    int rc = CHAMELEON_SUCCESS;

    memset( ipiv, 0, sizeof(CHAM_ipiv_t) );

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    ipiv->desc = desc;
    ipiv->data = data;
    ipiv->i    = 0;
    ipiv->m    = chameleon_min( desc->m, desc->n );
    ipiv->mb   = desc->mb;
    ipiv->mt   = chameleon_ceil( ipiv->m, ipiv->mb );

    /* Create runtime specific structure like registering data */
    RUNTIME_ipiv_create( ipiv );

    return rc;
}

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to destroy a tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to destroy.
 *
 */
void chameleon_ipiv_destroy( CHAM_ipiv_t *ipiv )
{
    RUNTIME_ipiv_destroy( ipiv );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Create a tiled ipiv descriptor associated to a given matrix.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to initialize.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated ipiv descriptor must be generated.
 *
 * @param[in] data
 *          The pointer to the original vector where to store the pivot values.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS on successful exit
 * @retval CHAMELEON_ERR_NOT_INITIALIZED if failed to initialize the descriptor.
 * @retval CHAMELEON_ERR_OUT_OF_RESOURCES if failed to allocated some ressources.
 *
 */
int CHAMELEON_Ipiv_Create( CHAM_ipiv_t **ipivptr, const CHAM_desc_t *desc, void *data )
{
    CHAM_context_t *chamctxt;
    CHAM_ipiv_t *ipiv;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Allocate memory and initialize the ipivriptor */
    ipiv = (CHAM_ipiv_t*)malloc(sizeof(CHAM_ipiv_t));
    if (ipiv == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Create", "malloc() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    chameleon_ipiv_init( ipiv, desc, data );

    *ipivptr = ipiv;
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Destroys an ipiv tile descriptor.
 *
 ******************************************************************************
 *
 * @param[in] ipivptr
 *          The Ipiv tile descriptor to destroy.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Destroy(CHAM_ipiv_t **ipivptr)
{
    CHAM_context_t *chamctxt;
    CHAM_ipiv_t *ipiv;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Destroy", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ((ipivptr == NULL) || (*ipivptr == NULL)) {
        chameleon_error("CHAMELEON_Ipiv_Destroy", "attempting to destroy a NULL descriptor");
        return CHAMELEON_ERR_UNALLOCATED;
    }

    ipiv = *ipivptr;
    chameleon_ipiv_destroy( ipiv );
    free(ipiv);
    *ipivptr = NULL;
    return CHAMELEON_SUCCESS;
}

 /**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Flushes the data in the sequence when they won't be reused. This calls
 * cleans up the distributed communication caches, and transfer the data back to
 * the CPU.
 *
 ******************************************************************************
 *
 * @param[in] ipiv
 *          ipiv vector descriptor.
 *
 * @param[in] sequence
 *          The seqeunce in which to submit the calls to flush the data.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Flush( const CHAM_ipiv_t        *ipiv,
                          const RUNTIME_sequence_t *sequence )
{
    RUNTIME_ipiv_flush( ipiv, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Gathers an IPIV tile descriptor in a single vector on the given root node.
 *
 ******************************************************************************
 *
 * @param[in] ipivdesc
 *          the ipiv vector descriptor to gather.
 *
 * @param[in] ipiv
 *          The ipiv vector where to store the result. Allocated vector of size
 *          ipivdesc->m on root, not referenced on other nodes.
 *
 * @param[in] root
 *          root node on which to gather the data.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Gather( CHAM_ipiv_t *ipivdesc, int *ipiv, int root )
{
    RUNTIME_ipiv_gather( ipivdesc, ipiv, root );
    return CHAMELEON_SUCCESS;
}
