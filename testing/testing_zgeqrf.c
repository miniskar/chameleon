/**
 *
 * @file testing_zgeqrf.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf testing
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
testing_zgeqrf_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_nb( args );
    int      ib     = run_arg_get_ib( args );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      RH     = run_arg_get_int( args, "qra", 4 );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descT;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgeqrf_Tile_Async( descA, descT, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgeqrf_Tile( descA, descT );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgeqrf( M, N ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descQ;
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, M, M, 0, 0, M, M, P, Q );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zungqr_Tile( descA, descT, descQ );

        hres += check_zgeqrf( args, descA0, descA, descQ );
        hres += check_zortho( args, descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );

    return hres;
}

int
testing_zgeqrf_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int nb    = run_arg_get_nb( args );
    int ib    = run_arg_get_ib( args );
    int P     = parameters_getvalue_int( "P" );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int LDA   = run_arg_get_int( args, "LDA", M );
    int RH    = run_arg_get_int( args, "qra", 4 );
    int seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int Q     = parameters_compute_q( P );
    int K     = chameleon_min( M, N );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;
    CHAM_desc_t           *descT;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zgeqrf( M, N, A, LDA, descT );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgeqrf( M, N ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAMELEON_Complex64_t *Qlap = malloc( M*M*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_Complex64_t *A0   = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

        CHAMELEON_zplrnt( M, N, A0, LDA, seedA );

        CHAMELEON_zungqr( M, M, K, A, LDA, descT, Qlap, M );

        hres += check_zgeqrf_std( args, M, M, K, A0, A, LDA, Qlap, M );
        hres += check_zortho_std( args, M, M, Qlap, M );

        free( A0 );
        free( Qlap );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descT );

    return hres;
}

testing_t   test_zgeqrf;
const char *zgeqrf_params[] = { "mtxfmt", "nb", "ib", "m", "n", "lda", "qra", "seedA", NULL };
const char *zgeqrf_output[] = { NULL };
const char *zgeqrf_outchk[] = { "||A||", "||I-QQ'||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgeqrf_init( void ) __attribute__( ( constructor ) );
void
testing_zgeqrf_init( void )
{
    test_zgeqrf.name   = "zgeqrf";
    test_zgeqrf.helper = "General QR factorization";
    test_zgeqrf.params = zgeqrf_params;
    test_zgeqrf.output = zgeqrf_output;
    test_zgeqrf.outchk = zgeqrf_outchk;
    test_zgeqrf.fptr_desc = testing_zgeqrf_desc;
    test_zgeqrf.fptr_std  = testing_zgeqrf_std;
    test_zgeqrf.next   = NULL;

    testing_register( &test_zgeqrf );
}
