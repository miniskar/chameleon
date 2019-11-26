/**
 *
 * @file testing_zunglq_hqr.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq_hqr testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-09-09
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zunglq_hqr( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descTS, *descTT, *descQ;

    /* Reads arguments */
    int    nb     = run_arg_get_int( args, "nb", 320 );
    int    ib     = run_arg_get_int( args, "ib", 48 );
    int    P      = parameters_getvalue_int( "P" );
    int    N      = run_arg_get_int( args, "N", 1000 );
    int    M      = run_arg_get_int( args, "M", N );
    int    K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int    LDA    = run_arg_get_int( args, "LDA", M );
    int    qr_a   = run_arg_get_int( args, "qra", -1 );
    int    qr_p   = run_arg_get_int( args, "qrp", -1 );
    int    llvl   = run_arg_get_int( args, "llvl", -1 );
    int    hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int    domino = run_arg_get_int( args, "domino", -1 );
    int    seedA  = run_arg_get_int( args, "seedA", random() );
    int    Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zunglq( M, N, K );

    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( K > chameleon_min( M, N ) ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameters for unglq_hqr (K > min(M,N))\n" );
        }
        return -1;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, K, N, P, Q );
    CHAMELEON_Desc_Create(
        &descQ, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, N, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, N, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_LQ, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    hres = CHAMELEON_zgelqf_param_Tile( &qrtree, descA, descTS, descTT );

    /* Calculates the solution */
    START_TIMING( t );
    CHAMELEON_zunglq_param_Tile( &qrtree, descA, descTS, descTT, descQ );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        hres += check_zortho( args, descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    CHAMELEON_Desc_Destroy( &descQ );
    libhqr_finalize( &qrtree );

    run_id++;
    return hres;
}

testing_t   test_zunglq_hqr;
const char *zunglq_hqr_params[] = { "nb",  "ib",   "m",    "n",      "k",     "lda", "qra",
                                    "qrp", "llvl", "hlvl", "domino", "seedA", NULL };
const char *zunglq_hqr_output[] = { NULL };
const char *zunglq_hqr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunglq_hqr_init( void ) __attribute__( ( constructor ) );
void
testing_zunglq_hqr_init( void )
{
    test_zunglq_hqr.name        = "zunglq_hqr";
    test_zunglq_hqr.helper      = "zunglq_hqr";
    test_zunglq_hqr.params      = zunglq_hqr_params;
    test_zunglq_hqr.output      = zunglq_hqr_output;
    test_zunglq_hqr.outchk      = zunglq_hqr_outchk;
    test_zunglq_hqr.params_list = "nb;ib;P;m;n;k;lda;qra;qrp;llvl;hlvl;domino;seedA";
    test_zunglq_hqr.fptr        = testing_zunglq_hqr;
    test_zunglq_hqr.next        = NULL;

    testing_register( &test_zunglq_hqr );
}
