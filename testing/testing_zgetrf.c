/**
 *
 * @file testing_zgetrf.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zgetrf( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_int( args, "nb", 320 );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      seedA  = run_arg_get_int( args, "seedA", random() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgetrf_nopiv_Tile_Async( descA, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgetrf_nopiv_Tile( descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrf( M, N ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        hres += check_zxxtrf( args, ChamGeneral, ChamUpperLower, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zgetrf;
const char *zgetrf_params[] = { "mtxfmt", "nb", "m", "n", "lda", "seedA", NULL };
const char *zgetrf_output[] = { NULL };
const char *zgetrf_outchk[] = { "||A||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgetrf_init( void ) __attribute__( ( constructor ) );
void
testing_zgetrf_init( void )
{
    test_zgetrf.name   = "zgetrf";
    test_zgetrf.helper = "General factorization (LU without pivoting)";
    test_zgetrf.params = zgetrf_params;
    test_zgetrf.output = zgetrf_output;
    test_zgetrf.outchk = zgetrf_outchk;
    test_zgetrf.fptr   = testing_zgetrf;
    test_zgetrf.next   = NULL;

    testing_register( &test_zgetrf );
}
