/**
 *
 * @file testing_zunmqr.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zunmqr_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_int( args, "nb", 320 );
    int          ib     = run_arg_get_int( args, "ib", 48 );
    int          P      = parameters_getvalue_int( "P" );
    cham_side_t  side   = run_arg_get_side( args, "side", ChamLeft );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int          LDA    = run_arg_get_int( args, "LDA", ( side == ChamLeft ) ? M : N );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    int          RH     = run_arg_get_int( args, "qra", 4 );
    int          seedA  = run_arg_get_int( args, "seedA", random() );
    int          seedC  = run_arg_get_int( args, "seedC", random() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    int          Am;
    CHAM_desc_t *descA, *descT, *descC;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Calculates the dimensions according to the transposition and the side */
    Am = ( side == ChamLeft ) ? M : N;

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, Am, K, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( Am, K, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Computes the factorization */
    hres = CHAMELEON_zgeqrf_Tile( descA, descT );
    assert( hres == 0 );

    /* Computes unmqr */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zunmqr_Tile_Async( side, trans, descA, descT, descC,
                                             test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zunmqr_Tile( side, trans, descA, descT, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zunmqr( side, M, N, K ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descC0 = CHAMELEON_Desc_Copy( descC, NULL );
        CHAM_desc_t *descQ;

        CHAMELEON_zplrnt_Tile( descC0, seedC );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, Am, Am, 0, 0, Am, Am, P, Q );
        CHAMELEON_zungqr_Tile( descA, descT, descQ );

        hres += check_zqc( args, side, trans, descC0, descQ, descC );

        CHAMELEON_Desc_Destroy( &descC0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

testing_t   test_zunmqr;
const char *zunmqr_params[] = { "mtxfmt", "nb",  "ib",  "side", "trans", "m",     "n",
                                "k",      "lda", "ldc", "qra",  "seedA", "seedC", NULL };
const char *zunmqr_output[] = { NULL };
const char *zunmqr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunmqr_init( void ) __attribute__( ( constructor ) );
void
testing_zunmqr_init( void )
{
    test_zunmqr.name   = "zunmqr";
    test_zunmqr.helper = "Q application (QR)";
    test_zunmqr.params = zunmqr_params;
    test_zunmqr.output = zunmqr_output;
    test_zunmqr.outchk = zunmqr_outchk;
    test_zunmqr.fptr_desc = testing_zunmqr_desc;
    test_zunmqr.fptr_std  = NULL;
    test_zunmqr.next   = NULL;

    testing_register( &test_zunmqr );
}
