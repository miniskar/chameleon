/**
 *
 * @file testing_zhemm.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-08-08
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zhemm( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          Am;
    int          hres = 0;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;

    /* Reads arguments */
    int                   nb    = run_arg_get_int( args, "nb", 320 );
    int                   P     = parameters_getvalue_int( "P" );
    cham_side_t           side  = run_arg_get_uplo( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( side == ChamLeft ) ? M : N ) );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    int                   LDC   = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    double                bump  = testing_dalea();
    bump                        = run_arg_get_double( args, "bump", bump );
    int    Q                    = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zhemm( side, M, N );

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the side */
    if ( side == ChamLeft ) {
        Am = M;
    }
    else {
        Am = N;
    }

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, Am, 0, 0, Am, Am, P, Q );
    CHAMELEON_Desc_Create(
        &descB, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descC, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( bump, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Calculates the product */
    START_TIMING( t );
    hres = CHAMELEON_zhemm_Tile( side, uplo, alpha, descA, descB, beta, descC );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descCinit, seedC );

        hres +=
            check_zsymm( args, ChamHermitian, side, uplo, alpha, descA, descB, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descC );

    run_id++;
    return hres;
}

testing_t   test_zhemm;
const char *zhemm_params[] = { "nb",    "side", "uplo",  "m",     "n",     "lda",  "ldb", "ldc",
                               "alpha", "beta", "seedA", "seedB", "seedC", "bump", NULL };
const char *zhemm_output[] = { NULL };
const char *zhemm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zhemm_init( void ) __attribute__( ( constructor ) );
void
testing_zhemm_init( void )
{
    test_zhemm.name        = "zhemm";
    test_zhemm.helper      = "zhemm";
    test_zhemm.params      = zhemm_params;
    test_zhemm.output      = zhemm_output;
    test_zhemm.outchk      = zhemm_outchk;
    test_zhemm.params_list = "nb;P;side;uplo;m;n;lda;ldb;ldc;alpha;beta;seedA;seedB;seedC;bump";
    test_zhemm.fptr        = testing_zhemm;
    test_zhemm.next        = NULL;

    testing_register( &test_zhemm );
}
