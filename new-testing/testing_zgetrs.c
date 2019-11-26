/**
 *
 * @file testing_zgetrs.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrs testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-09-09
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zgetrs( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres;
    CHAM_desc_t *descA, *descX;

    /* Reads arguments */
    int    nb    = run_arg_get_int( args, "nb", 320 );
    int    P     = parameters_getvalue_int( "P" );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int    LDA   = run_arg_get_int( args, "LDA", N );
    int    LDB   = run_arg_get_int( args, "LDB", N );
    int    seedA = run_arg_get_int( args, "seedA", random() );
    int    seedB = run_arg_get_int( args, "seedB", random() );
    int    Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgetrs( N, NRHS );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    hres = CHAMELEON_zgetrf_nopiv_Tile( descA );
    assert( hres == 0 );

    /* Calculates the solution */
    START_TIMING( t );
    hres += CHAMELEON_zgetrs_nopiv_Tile( descA, descX );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, NULL );

        CHAMELEON_zplrnt_Tile( descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres += check_zsolve( args, ChamGeneral, ChamNoTrans, ChamUpperLower, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    run_id++;
    return hres;
}

testing_t   test_zgetrs;
const char *zgetrs_params[] = { "nb", "n", "nrhs", "lda", "ldb", "seedA", "seedB", NULL };
const char *zgetrs_output[] = { NULL };
const char *zgetrs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgetrs_init( void ) __attribute__( ( constructor ) );
void
testing_zgetrs_init( void )
{
    test_zgetrs.name        = "zgetrs";
    test_zgetrs.helper      = "zgetrs";
    test_zgetrs.params      = zgetrs_params;
    test_zgetrs.output      = zgetrs_output;
    test_zgetrs.outchk      = zgetrs_outchk;
    test_zgetrs.params_list = "nb;P;n;nrhs;lda;ldb;seedA;seedB";
    test_zgetrs.fptr        = testing_zgetrs;
    test_zgetrs.next        = NULL;

    testing_register( &test_zgetrs );
}
