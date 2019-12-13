/**
 *
 * @file testing_zgetrf.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf testing
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
testing_zgetrf( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA;

    /* Reads arguments */
    int    nb    = run_arg_get_int( args, "nb", 320 );
    int    P     = parameters_getvalue_int( "P" );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    seedA = run_arg_get_int( args, "seedA", random() );
    int    Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgetrf( M, N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the solution */
    START_TIMING( t );
    hres = CHAMELEON_zgetrf_nopiv_Tile( descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        hres += check_zxxtrf( args, ChamGeneral, ChamUpperLower, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    run_id++;
    return hres;
}

testing_t   test_zgetrf;
const char *zgetrf_params[] = { "nb", "m", "n", "lda", "seedA", NULL };
const char *zgetrf_output[] = { NULL };
const char *zgetrf_outchk[] = { "||A||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgetrf_init( void ) __attribute__( ( constructor ) );
void
testing_zgetrf_init( void )
{
    test_zgetrf.name        = "zgetrf";
    test_zgetrf.helper      = "zgetrf";
    test_zgetrf.params      = zgetrf_params;
    test_zgetrf.output      = zgetrf_output;
    test_zgetrf.outchk      = zgetrf_outchk;
    test_zgetrf.params_list = "nb;P;m;n;lda;seedA";
    test_zgetrf.fptr        = testing_zgetrf;
    test_zgetrf.next        = NULL;

    testing_register( &test_zgetrf );
}
