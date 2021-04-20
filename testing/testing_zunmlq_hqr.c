/**
 *
 * @file testing_zunmlq_hqr.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq_hqr testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
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
testing_zunmlq_hqr( run_arg_list_t *args, int check )
{
    int          An;
    int          hres;
    CHAM_desc_t *descA, *descTS, *descTT, *descC;

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
    int          qr_a   = run_arg_get_int( args, "qra", -1 );
    int          qr_p   = run_arg_get_int( args, "qrp", -1 );
    int          llvl   = run_arg_get_int( args, "llvl", -1 );
    int          hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int          domino = run_arg_get_int( args, "domino", -1 );
    int          seedA  = run_arg_get_int( args, "seedA", random() );
    int          seedC  = run_arg_get_int( args, "seedC", random() );
    int          Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zunmlq( side, M, N, K );

    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Calculates the dimensions according to the transposition and the side */
    An = ( side == ChamLeft ) ? M : N;

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, K, An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, An, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, An, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_LQ, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Computes the factorization */
    hres = CHAMELEON_zgelqf_param_Tile( &qrtree, descA, descTS, descTT );
    assert( hres == 0 );

    /* Computes unmlq_hqr */
    START_TIMING( t );
    hres += CHAMELEON_zunmlq_param_Tile( &qrtree, side, trans, descA, descTS, descTT, descC );
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
        CHAMELEON_zunglq_param_Tile( &qrtree, descA, descTS, descTT, descQ );

        hres += check_zqc( args, side, trans, descC0, descQ, descC );

        CHAMELEON_Desc_Destroy( &descC0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    CHAMELEON_Desc_Destroy( &descC );
    libhqr_finalize( &qrtree );

    return hres;
}

testing_t   test_zunmlq_hqr;
const char *zunmlq_hqr_params[] = { "mtxfmt", "nb",  "ib",     "side",  "trans", "m",   "n",
                                    "k",    "lda",    "ldc",   "qra",   "qrp", "llvl",
                                    "hlvl", "domino", "seedA", "seedC", NULL };
const char *zunmlq_hqr_output[] = { NULL };
const char *zunmlq_hqr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunmlq_hqr_init( void ) __attribute__( ( constructor ) );
void
testing_zunmlq_hqr_init( void )
{
    test_zunmlq_hqr.name   = "zunmlq_hqr";
    test_zunmlq_hqr.helper = "Q application with hierarchical reduction trees (LQ)";
    test_zunmlq_hqr.params = zunmlq_hqr_params;
    test_zunmlq_hqr.output = zunmlq_hqr_output;
    test_zunmlq_hqr.outchk = zunmlq_hqr_outchk;
    test_zunmlq_hqr.fptr = testing_zunmlq_hqr;
    test_zunmlq_hqr.next = NULL;

    testing_register( &test_zunmlq_hqr );
}
