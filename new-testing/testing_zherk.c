/**
 *
 * @file testing_zherk.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-08-09
 * @precisions normal z -> z c
 *
 */
#include <chameleon.h>
#include "flops.h"
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zherk( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          Am, An;
    int          hres = 0;
    CHAM_desc_t *descA, *descC, *descCinit;

    /* Reads arguments */
    int          nb    = run_arg_get_int( args, "nb", 320 );
    int          P     = parameters_getvalue_int( "P" );
    cham_trans_t trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t  uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int          N     = run_arg_get_int( args, "N", 1000 );
    int          K     = run_arg_get_int( args, "K", N );
    int          LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int          LDC   = run_arg_get_int( args, "LDC", N );
    double       alpha = testing_dalea();
    double       beta  = testing_dalea();
    double       bump  = testing_dalea();
    int          seedA = run_arg_get_int( args, "seedA", random() );
    int          seedC = run_arg_get_int( args, "seedC", random() );
    int          Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zherk( K, N );

    alpha = run_arg_get_double( args, "alpha", alpha );
    beta  = run_arg_get_double( args, "beta", beta );
    bump  = run_arg_get_double( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplghe_Tile( bump, uplo, descC, seedC );

    /* Calculates the product */
    START_TIMING( t );
    hres = CHAMELEON_zherk_Tile( uplo, trans, alpha, descA, beta, descC );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplghe_Tile( bump, uplo, descCinit, seedC );

        hres +=
            check_zsyrk( args, ChamHermitian, uplo, trans, alpha, descA, NULL, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descC );

    run_id++;
    return hres;
}

testing_t   test_zherk;
const char *zherk_params[] = { "nb",    "trans", "uplo",  "n",     "k",    "lda", "ldc",
                               "alpha", "beta",  "seedA", "seedC", "bump", NULL };
const char *zherk_output[] = { NULL };
const char *zherk_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zherk_init( void ) __attribute__( ( constructor ) );
void
testing_zherk_init( void )
{
    test_zherk.name        = "zherk";
    test_zherk.helper      = "zherk";
    test_zherk.params      = zherk_params;
    test_zherk.output      = zherk_output;
    test_zherk.outchk      = zherk_outchk;
    test_zherk.params_list = "nb;P;trans;uplo;n;k;lda;ldc;alpha;beta;seedA;seedC;bump";
    test_zherk.fptr        = testing_zherk;
    test_zherk.next        = NULL;

    testing_register( &test_zherk );
}
