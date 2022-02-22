/**
 *
 * @file testing_zplrnk.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnk testing
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zplrnk_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int async = parameters_getvalue_int( "async" );
    int nb    = run_arg_get_int( args, "nb", 320 );
    int P     = parameters_getvalue_int( "P" );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int K     = run_arg_get_int( args, "K", N );
    int LDC   = run_arg_get_int( args, "LDC", M );
    int seedA = run_arg_get_int( args, "seedA", random() );
    int seedB = run_arg_get_int( args, "seedB", random() );
    int Q     = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descC;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descC, NULL, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );

    /* Calculates the random rank-k matrix */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zplrnk_Tile_Async( K, descC, seedA, seedB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zplrnk_Tile( K, descC, seedA, seedB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm( M, N, K ) );

    /* Checks the solution */
    if ( check ) {
        hres = check_zrankk( args, K, descC );
    }

    CHAMELEON_Desc_Destroy( &descC );

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
    test_zplrnk.name   = "zplrnk";
    test_zplrnk.helper = "General rank-k matrix generation";
    test_zplrnk.params = zplrnk_params;
    test_zplrnk.output = zplrnk_output;
    test_zplrnk.outchk = zplrnk_outchk;
    test_zplrnk.fptr_desc = testing_zplrnk_desc;
    test_zplrnk.fptr_std  = NULL;
    test_zplrnk.next   = NULL;

    testing_register( &test_zplrnk );
}
