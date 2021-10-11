/**
 *
 * @file descriptor_rec.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Gwenole Lucas
 * @date 2020-03-03
 *
 */
#define _GNU_SOURCE 1
#include "control/common.h"
#include <stdlib.h>
#include <stdio.h>
#include "chameleon/runtime.h"

static int
chameleon_recdesc_create( const char *name, CHAM_desc_t **descptr, void *mat, cham_flttype_t dtyp,
                          int *mb, int *nb,
                          int lm, int ln, int m, int n, int p, int q,
                          blkaddr_fct_t get_blkaddr, blkldd_fct_t get_blkldd, blkrankof_fct_t get_rankof )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t    *desc;
    int rc;

    *descptr = NULL;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_error( "CHAMELEON_Recursive_Desc_Create", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( chamctxt->scheduler != RUNTIME_SCHED_STARPU ) {
        chameleon_error( "CHAMELEON_Recursive_Desc_Create", "CHAMELEON Recursive descriptors only available with StaRPU" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Let's make sure we have at least one couple (mb, nb) defined */
    assert( (mb[0] > 0) && (nb[0] > 0) );

    /* Create the current layer descriptor */
    desc = (CHAM_desc_t*)malloc(sizeof(CHAM_desc_t));
    rc = chameleon_desc_init_internal( desc, name, mat, dtyp, mb[0], nb[0],
                                       lm, ln, m, n, p, q,
                                       get_blkaddr, get_blkldd, get_rankof );
    *descptr = desc;

    if ( rc != CHAMELEON_SUCCESS ) {
        return rc;
    }

    /* Move to the next tile size to recurse */
    mb++;
    nb++;
    if ( (mb[0] <= 0) || (nb[0] <= 0) ) {
        return CHAMELEON_SUCCESS;
    }

    for ( n=0; n<desc->nt; n++ ) {
        for ( m=0; m<desc->mt; m++ ) {
            CHAM_desc_t *tiledesc;
            CHAM_tile_t *tile;
            int tempmm, tempnn;
            char *subname;

            tile = desc->get_blktile( desc, m, n );
            tempmm = m == desc->mt-1 ? desc->m - m * desc->mb : desc->mb;
            tempnn = n == desc->nt-1 ? desc->n - n * desc->nb : desc->nb;

            chameleon_asprintf( &subname, "%s[%d,%d]", name, m, n );

            chameleon_recdesc_create( subname, &tiledesc, tile->mat,
                                      desc->dtyp, mb, nb,
                                      tile->ld, tempnn, /* Abuse as ln is not used */
                                      tempmm, tempnn,
                                      1, 1,             /* can recurse only on local data */
                                      chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL);

            tile->format = CHAMELEON_TILE_DESC;
            tile->mat = tiledesc;

            if ( rc != CHAMELEON_SUCCESS ) {
                return rc;
            }
        }
    }

    return CHAMELEON_SUCCESS;
}

int
CHAMELEON_Recursive_Desc_Create( CHAM_desc_t **descptr, void *mat, cham_flttype_t dtyp,
                                 int *mb, int *nb, int lm, int ln, int m, int n, int p, int q,
                                 blkaddr_fct_t get_blkaddr, blkldd_fct_t get_blkldd, blkrankof_fct_t get_rankof )
{
    /*
     * The first layer must be allocated, otherwise we will give unitialized
     * pointers to the lower layers
     */
    assert( (mat != CHAMELEON_MAT_ALLOC_TILE) &&
            (mat != CHAMELEON_MAT_OOC) );

    return chameleon_recdesc_create( "A", descptr, mat, dtyp,
                                     mb, nb, lm, ln, m, n, p, q,
                                     get_blkaddr, get_blkldd, get_rankof );
}
