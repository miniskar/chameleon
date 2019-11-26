/**
 *
 * @file testing_zher2k.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-08-09
 * @precisions normal z -> z c
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zher2k( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          Am, An;
    int          hres = 0;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;

    /* Read arguments */
    int          nb    = run_arg_get_int( args, "nb", 320 );
    int          P     = parameters_getvalue_int( "P" );
    cham_trans_t trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t  uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int          N     = run_arg_get_int( args, "N", 1000 );
    int          K     = run_arg_get_int( args, "K", N );
    int          LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int          LDB   = run_arg_get_int( args, "LDB", ( ( trans == ChamNoTrans ) ? N : K ) );
    int          LDC   = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    double                beta  = testing_dalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    double                bump  = testing_dalea();
    bump                        = run_arg_get_double( args, "bump", bump );
    int    Q                    = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zher2k( K, N );

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_double( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descB, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplghe_Tile( bump, uplo, descC, seedC );

    /* Calculate the product */
    START_TIMING( t );
    hres = CHAMELEON_zher2k_Tile( uplo, trans, alpha, descA, descB, beta, descC );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", gflops );

    /* Check the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplghe_Tile( bump, uplo, descCinit, seedC );

        hres +=
            check_zsyrk( args, ChamHermitian, uplo, trans, alpha, descA, descB, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descC );

    run_id++;
    return hres;
}

testing_t   test_zher2k;
const char *zher2k_params[] = { "nb",    "trans", "uplo",  "n",     "k",     "lda",  "ldb", "ldc",
                                "alpha", "beta",  "seedA", "seedB", "seedC", "bump", NULL };
const char *zher2k_output[] = { NULL };
const char *zher2k_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zher2k_init( void ) __attribute__( ( constructor ) );
void
testing_zher2k_init( void )
{
    test_zher2k.name        = "zher2k";
    test_zher2k.helper      = "zher2k";
    test_zher2k.params      = zher2k_params;
    test_zher2k.output      = zher2k_output;
    test_zher2k.outchk      = zher2k_outchk;
    test_zher2k.params_list = "nb;P;trans;uplo;n;k;lda;ldb;ldc;alpha;beta;seedA;seedB;seedC;bump";
    test_zher2k.fptr        = testing_zher2k;
    test_zher2k.next        = NULL;

    testing_register( &test_zher2k );
}
