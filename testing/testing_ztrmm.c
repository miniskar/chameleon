/**
 *
 * @file testing_ztrmm.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm testing
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
testing_ztrmm( run_arg_list_t *args, int check )
{
    int          Bm, Bn;
    int          hres = 0;
    CHAM_desc_t *descA, *descB, *descBinit;

    /* Reads arguments */
    intptr_t              mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int                   nb     = run_arg_get_int( args, "nb", 320 );
    int                   P      = parameters_getvalue_int( "P" );
    cham_trans_t          trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t           side   = run_arg_get_uplo( args, "side", ChamLeft );
    cham_uplo_t           uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t           diag   = run_arg_get_diag( args, "diag", ChamNonUnit );
    int                   N      = run_arg_get_int( args, "N", 1000 );
    int                   K      = run_arg_get_int( args, "K", N );
    int                   LDA    = run_arg_get_int( args, "LDA", N );
    int                   LDB    = run_arg_get_int( args, "LDB", N );
    CHAMELEON_Complex64_t alpha  = testing_zalea();
    int                   seedA  = run_arg_get_int( args, "seedA", random() );
    int                   seedB  = run_arg_get_int( args, "seedB", random() );
    int                   Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_ztrmm( side, N, K );

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the side */
    if ( side == ChamLeft ) {
        Bm = N;
        Bn = K;
    }
    else {
        Bm = K;
        Bn = N;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create(
        &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, Bn, 0, 0, Bm, Bn, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Calculates the product */
    START_TIMING( t );
    hres = CHAMELEON_ztrmm_Tile( side, uplo, trans, diag, alpha, descA, descB );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descBinit, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, Bn, 0, 0, Bm, Bn, P, Q );
        CHAMELEON_zplrnt_Tile( descBinit, seedB );

        hres += check_ztrmm( args, CHECK_TRMM, side, uplo, trans, diag, alpha, descA, descB, descBinit );

        CHAMELEON_Desc_Destroy( &descBinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );

    return hres;
}

testing_t   test_ztrmm;
const char *ztrmm_params[] = { "mtxfmt", "nb", "trans", "side",  "uplo",  "diag",  "n", "k",
                               "lda", "ldb",   "alpha", "seedA", "seedB", NULL };
const char *ztrmm_output[] = { NULL };
const char *ztrmm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrmm_init( void ) __attribute__( ( constructor ) );
void
testing_ztrmm_init( void )
{
    test_ztrmm.name        = "ztrmm";
    test_ztrmm.helper      = "Triangular matrix-matrix multiply";
    test_ztrmm.params      = ztrmm_params;
    test_ztrmm.output      = ztrmm_output;
    test_ztrmm.outchk      = ztrmm_outchk;
    test_ztrmm.fptr        = testing_ztrmm;
    test_ztrmm.next        = NULL;

    testing_register( &test_ztrmm );
}
