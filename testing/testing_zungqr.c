/**
 *
 * @file testing_zungqr.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqr testing
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zungqr( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descT, *descQ;

    /* Reads arguments */
    int    nb    = run_arg_get_int( args, "nb", 320 );
    int    ib    = run_arg_get_int( args, "ib", 48 );
    int    P     = parameters_getvalue_int( "P" );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    K     = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    RH    = run_arg_get_int( args, "qra", 0 );
    int    seedA = run_arg_get_int( args, "seedA", random() );
    int    Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zungqr( M, N, K );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( K > chameleon_min( M, N ) ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameters for ungqr (K > min(M,N))\n" );
        }
        return -1;
    }

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, M, K, P, Q );
    CHAMELEON_Desc_Create(
        &descQ, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, K, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    hres = CHAMELEON_zgeqrf_Tile( descA, descT );

    /* Calculates the solution */
    START_TIMING( t );
    CHAMELEON_zungqr_Tile( descA, descT, descQ );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        hres += check_zortho( args, descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descQ );

    run_id++;
    return hres;
}

testing_t   test_zungqr;
const char *zungqr_params[] = { "nb", "ib", "m", "n", "k", "lda", "qra", "seedA", NULL };
const char *zungqr_output[] = { NULL };
const char *zungqr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zungqr_init( void ) __attribute__( ( constructor ) );
void
testing_zungqr_init( void )
{
    test_zungqr.name        = "zungqr";
    test_zungqr.helper      = "Q generation (QR)";
    test_zungqr.params      = zungqr_params;
    test_zungqr.output      = zungqr_output;
    test_zungqr.outchk      = zungqr_outchk;
    test_zungqr.fptr        = testing_zungqr;
    test_zungqr.next        = NULL;

    testing_register( &test_zungqr );
}
