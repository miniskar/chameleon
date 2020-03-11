/**
 *
 * @file testing_ztrsm.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm testing
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
testing_ztrsm( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres = 0;
    CHAM_desc_t *descA, *descB, *descBinit;

    /* Reads arguments */
    int                   nb    = run_arg_get_int( args, "nb", 320 );
    int                   P     = parameters_getvalue_int( "P" );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t           side  = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t           diag  = run_arg_get_diag( args, "diag", ChamNonUnit );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   Ak    = ( side == ChamLeft ) ? M : N;
    int                   LDA   = run_arg_get_int( args, "LDA", Ak );
    int                   LDB   = run_arg_get_int( args, "LDB", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_ztrsm( side, M, N );

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, Ak, 0, 0, Ak, Ak, P, Q );
    CHAMELEON_Desc_Create(
        &descB, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    /* We bump a little bit the diagonal to make it stable */
    CHAMELEON_zplgsy_Tile( 2., uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Calculates the product */
    START_TIMING( t );
    hres = CHAMELEON_ztrsm_Tile( side, uplo, trans, diag, alpha, descA, descB );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descBinit, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descBinit, seedB );

        hres += check_ztrmm( args, CHECK_TRSM, side, uplo, trans, diag, alpha, descA, descB, descBinit );

        CHAMELEON_Desc_Destroy( &descBinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );

    run_id++;
    return hres;
}

testing_t   test_ztrsm;
const char *ztrsm_params[] = { "nb",  "side",  "uplo",  "trans", "diag",  "m", "n",
                               "lda", "ldb",   "alpha", "seedA", "seedB", NULL };
const char *ztrsm_output[] = { NULL };
const char *ztrsm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrsm_init( void ) __attribute__( ( constructor ) );
void
testing_ztrsm_init( void )
{
    test_ztrsm.name        = "ztrsm";
    test_ztrsm.helper      = "Triangular matrix solve";
    test_ztrsm.params      = ztrsm_params;
    test_ztrsm.output      = ztrsm_output;
    test_ztrsm.outchk      = ztrsm_outchk;
    test_ztrsm.params_list = "nb;P;trans;side;uplo;diag;n;k;lda;ldb;alpha;seedA;seedB";
    test_ztrsm.fptr        = testing_ztrsm;
    test_ztrsm.next        = NULL;

    testing_register( &test_ztrsm );
}
