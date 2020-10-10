/**
 *
 * @file testing_zplrnk.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnk testing
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @date 2020-03-05
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zplrnk( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descC;

    /* Reads arguments */
    int             nb    = run_arg_get_int( args, "nb", 320 );
    int             P     = parameters_getvalue_int( "P" );
    int             N     = run_arg_get_int( args, "N", 1000 );
    int             M     = run_arg_get_int( args, "M", N );
    int             K     = run_arg_get_int( args, "K", N );
    int             LDC   = run_arg_get_int( args, "LDC", M );
    int             seedA = run_arg_get_int( args, "seedA", random() );
    int             seedB = run_arg_get_int( args, "seedB", random() );
    int             Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    /* We consider the gemm cost used in this operation as the cost */
    cham_fixdbl_t flops = flops_zgemm( M, N, K );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descC, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );

    /* Calculates the random rank-k matrix */
    START_TIMING( t );
    hres = CHAMELEON_zplrnk_Tile( K, descC, seedA, seedB );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        hres = check_zrankk( args, K, descC );
    }

    CHAMELEON_Desc_Destroy( &descC );

    run_id++;
    return hres;
}

testing_t   test_zplrnk;
const char *zplrnk_params[] = { "nb", "m", "n", "k", "ldc", "seedA", "seedB", NULL };
const char *zplrnk_output[] = { NULL };
const char *zplrnk_outchk[] = { "||A||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zplrnk_init( void ) __attribute__( ( constructor ) );
void
testing_zplrnk_init( void )
{
    test_zplrnk.name        = "zplrnk";
    test_zplrnk.helper      = "General rank-k matrix generation";
    test_zplrnk.params      = zplrnk_params;
    test_zplrnk.output      = zplrnk_output;
    test_zplrnk.outchk      = zplrnk_outchk;
    test_zplrnk.fptr        = testing_zplrnk;
    test_zplrnk.next        = NULL;

    testing_register( &test_zplrnk );
}
