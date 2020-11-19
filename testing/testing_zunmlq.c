/**
 *
 * @file testing_zunmlq.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq testing
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zunmlq( run_arg_list_t *args, int check )
{
    int          An;
    int          hres;
    CHAM_desc_t *descA, *descT, *descC;

    /* Reads arguments */
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_int( args, "nb", 320 );
    int          ib     = run_arg_get_int( args, "ib", 48 );
    int          P      = parameters_getvalue_int( "P" );
    cham_side_t  side   = run_arg_get_uplo( args, "side", ChamLeft );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", N );
    int          LDA    = run_arg_get_int( args, "LDA", K );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    int          RH     = run_arg_get_int( args, "qra", 4 );
    int          seedA  = run_arg_get_int( args, "seedA", random() );
    int          seedC  = run_arg_get_int( args, "seedC", random() );
    int          Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zunmlq( side, M, N, K );

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
    An = ( side == ChamLeft ) ? M : N;

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, K, An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, An, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Computes the factorization */
    hres = CHAMELEON_zgelqf_Tile( descA, descT );
    assert( hres == 0 );

    /* Computes unmlq */
    START_TIMING( t );
    hres += CHAMELEON_zunmlq_Tile( side, trans, descA, descT, descC );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descC0 = CHAMELEON_Desc_Copy( descC, NULL );
        CHAM_desc_t *descQ;

        CHAMELEON_zplrnt_Tile( descC0, seedC );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, An, An, 0, 0, An, An, P, Q );
        CHAMELEON_zunglq_Tile( descA, descT, descQ );

        hres += check_zqc( args, side, trans, descC0, descQ, descC );

        CHAMELEON_Desc_Destroy( &descC0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

testing_t   test_zunmlq;
const char *zunmlq_params[] = { "mtxfmt", "nb", "ib",  "side", "trans", "m",     "n", "k",
                                "lda", "ldc", "qra",   "seedA", "seedC", NULL };
const char *zunmlq_output[] = { NULL };
const char *zunmlq_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunmlq_init( void ) __attribute__( ( constructor ) );
void
testing_zunmlq_init( void )
{
    test_zunmlq.name        = "zunmlq";
    test_zunmlq.helper      = "Q application (LQ)";
    test_zunmlq.params      = zunmlq_params;
    test_zunmlq.output      = zunmlq_output;
    test_zunmlq.outchk      = zunmlq_outchk;
    test_zunmlq.fptr        = testing_zunmlq;
    test_zunmlq.next        = NULL;

    testing_register( &test_zunmlq );
}
