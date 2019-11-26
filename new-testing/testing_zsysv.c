/**
 *
 * @file testing_zsysv.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsysv testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-08-12
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

static cham_fixdbl_t
flops_zsysv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotrs( N, NRHS );
    return flops;
}

int
testing_zsysv( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descX;

    /* Reads arguments */
    int         nb    = run_arg_get_int( args, "nb", 320 );
    int         P     = parameters_getvalue_int( "P" );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         LDB   = run_arg_get_int( args, "LDB", N );
    int         seedA = run_arg_get_int( args, "seedA", random() );
    int         seedB = run_arg_get_int( args, "seedB", random() );
    int         Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zsysv( N, NRHS );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy_Tile( (double)N, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    START_TIMING( t );
    hres = CHAMELEON_zsysv_Tile( uplo, descA, descX );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );

        hres += check_zxxtrf( args, ChamSymmetric, uplo, descA0, descA );

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, NULL );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );
        hres += check_zsolve( args, ChamSymmetric, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    run_id++;
    return hres;
}

testing_t   test_zsysv;
const char *zsysv_params[] = { "nb", "uplo", "n", "nrhs", "lda", "ldb", "seedA", "seedB", NULL };
const char *zsysv_output[] = { NULL };
const char *zsysv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsysv_init( void ) __attribute__( ( constructor ) );
void
testing_zsysv_init( void )
{
    test_zsysv.name        = "zsysv";
    test_zsysv.helper      = "zsysv";
    test_zsysv.params      = zsysv_params;
    test_zsysv.output      = zsysv_output;
    test_zsysv.outchk      = zsysv_outchk;
    test_zsysv.params_list = "nb;P;uplo;n;nrhs;lda;ldb;seedA;seedB";
    test_zsysv.fptr        = testing_zsysv;
    test_zsysv.next        = NULL;

    testing_register( &test_zsysv );
}
