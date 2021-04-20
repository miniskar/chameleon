/**
 *
 * @file testing_zsymm.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsymm testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zsymm( run_arg_list_t *args, int check )
{
    int          Am;
    int          hres = 0;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;

    /* Reads arguments */
    intptr_t              mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int                   nb     = run_arg_get_int( args, "nb", 320 );
    int                   P      = parameters_getvalue_int( "P" );
    cham_side_t           side   = run_arg_get_uplo( args, "side", ChamLeft );
    cham_uplo_t           uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N      = run_arg_get_int( args, "N", 1000 );
    int                   M      = run_arg_get_int( args, "M", N );
    int                   LDA    = run_arg_get_int( args, "LDA", ( ( side == ChamLeft ) ? M : N ) );
    int                   LDB    = run_arg_get_int( args, "LDB", M );
    int                   LDC    = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha  = testing_zalea();
    CHAMELEON_Complex64_t beta   = testing_zalea();
    int                   seedA  = run_arg_get_int( args, "seedA", random() );
    int                   seedB  = run_arg_get_int( args, "seedB", random() );
    int                   seedC  = run_arg_get_int( args, "seedC", random() );
    double                bump   = testing_dalea();
    bump                         = run_arg_get_double( args, "bump", bump );
    int                   Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zsymm( side, M, N );

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
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, Am, 0, 0, Am, Am, P, Q );
    CHAMELEON_Desc_Create(
        &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy_Tile( bump, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Calculates the product */
    START_TIMING( t );
    hres = CHAMELEON_zsymm_Tile( side, uplo, alpha, descA, descB, beta, descC );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descCinit, seedC );

        hres +=
            check_zsymm( args, ChamSymmetric, side, uplo, alpha, descA, descB, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

testing_t   test_zsymm;
const char *zsymm_params[] = { "mtxfmt", "nb",   "side", "uplo",  "m",     "n",     "lda",  "ldb", "ldc",
                               "alpha", "beta", "seedA", "seedB", "seedC", "bump", NULL };
const char *zsymm_output[] = { NULL };
const char *zsymm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsymm_init( void ) __attribute__( ( constructor ) );
void
testing_zsymm_init( void )
{
    test_zsymm.name        = "zsymm";
    test_zsymm.helper      = "Symmetric matrix-matrix multiply";
    test_zsymm.params      = zsymm_params;
    test_zsymm.output      = zsymm_output;
    test_zsymm.outchk      = zsymm_outchk;
    test_zsymm.fptr        = testing_zsymm;
    test_zsymm.next        = NULL;

    testing_register( &test_zsymm );
}
