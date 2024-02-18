/**
 *
 * @file testing_zgetrf_nopiv.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv testing
 *
 * @version 1.3.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2023-07-05
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zgetrf_nopiv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int    async = parameters_getvalue_int( "async" );
    int    nb    = run_arg_get_nb( args );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    double bump  = run_arg_get_double( args, "bump", (double)N );

    /* Descriptors */
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplgtr_Tile( 0,    ChamUpper, descA, seedA   );
    CHAMELEON_zplgtr_Tile( bump, ChamLower, descA, seedA+1 );

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
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplgtr_Tile( 0,    ChamUpper, descA0, seedA   );
        CHAMELEON_zplgtr_Tile( bump, ChamLower, descA0, seedA+1 );

        hres += check_zxxtrf( args, ChamGeneral, ChamUpperLower, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

    return hres;
}

int
testing_zgetrf_nopiv_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int    nb    = run_arg_get_nb( args );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    double bump  = run_arg_get_double( args, "bump", (double)N );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplgtr( 0,    ChamUpper, M, N, A, LDA, seedA   );
    CHAMELEON_zplgtr( bump, ChamLower, M, N, A, LDA, seedA+1 );

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zgetrf_nopiv( M, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrf( M, N ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAMELEON_Complex64_t *A0 = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplgtr( 0,    ChamUpper, M, N, A0, LDA, seedA   );
        CHAMELEON_zplgtr( bump, ChamLower, M, N, A0, LDA, seedA+1 );

        hres += check_zxxtrf_std( args, ChamGeneral, ChamUpperLower, M, N, A0, A, LDA );

        free( A0 );
    }

    free( A );

    return hres;
}

testing_t   test_zgetrf_nopiv;
const char *zgetrf_nopiv_params[] = { "mtxfmt", "nb", "m", "n", "lda", "seedA", "bump", NULL };
const char *zgetrf_nopiv_output[] = { NULL };
const char *zgetrf_nopiv_outchk[] = { "||A||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgetrf_nopiv_init( void ) __attribute__( ( constructor ) );
void
testing_zgetrf_nopiv_init( void )
{
    test_zgetrf_nopiv.name   = "zgetrf_nopiv";
    test_zgetrf_nopiv.helper = "General factorization (LU without pivoting)";
    test_zgetrf_nopiv.params = zgetrf_nopiv_params;
    test_zgetrf_nopiv.output = zgetrf_nopiv_output;
    test_zgetrf_nopiv.outchk = zgetrf_nopiv_outchk;
    test_zgetrf_nopiv.fptr_desc = testing_zgetrf_nopiv_desc;
    test_zgetrf_nopiv.fptr_std  = testing_zgetrf_nopiv_std;
    test_zgetrf_nopiv.next   = NULL;

    testing_register( &test_zgetrf_nopiv );
}
