/**
 *
 * @file testing_zpoinv.c
 *
 * @copyright 2019-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpoinv testing
 *
 * @version 1.3.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2023-07-05
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if defined(CHAMELEON_TESTINGS_VENDOR)
#include <coreblas/lapacke.h>
#include <coreblas.h>
#endif

static cham_fixdbl_t
flops_zpoinv( int N )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotri( N );
    return flops;
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zpoinv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );

    /* Initialise the matrix with the random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the inversed matrix */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zpoinv_Tile_Async( uplo, descA, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zpoinv_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpoinv( N ) );

    /* Check the inverse */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_ztrtri( args, ChamHermitian, uplo, ChamNonUnit, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

    return hres;
}
#endif

int
testing_zpoinv_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

    /* Initialise the matrix with the random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );

    /* Calculates the inversed matrix */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres =  LAPACKE_zpotrf( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
    hres += LAPACKE_zpotri( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpoinv( N ) );
#else
    testing_start( &test_data );
    hres = CHAMELEON_zpoinv( uplo, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpoinv( N ) );

    /* Check the inverse */
    if ( check ) {
        CHAMELEON_Complex64_t *A0 = malloc( LDA * N * sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_ztrtri_std( args, ChamHermitian, uplo, ChamNonUnit, N, A0, A, LDA );

        free( A0 );
    }
#endif

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zpoinv;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zpoinv_params[] = { "uplo", "n", "lda", "seedA", NULL };
#else
const char *zpoinv_params[] = { "mtxfmt", "nb", "uplo", "n", "lda", "seedA", NULL };
#endif
const char *zpoinv_output[] = { NULL };
const char *zpoinv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpoinv_init( void ) __attribute__( ( constructor ) );
void
testing_zpoinv_init( void )
{
    test_zpoinv.name   = "zpoinv";
    test_zpoinv.helper = "Hermitian positive definite matrix inversion";
    test_zpoinv.params = zpoinv_params;
    test_zpoinv.output = zpoinv_output;
    test_zpoinv.outchk = zpoinv_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zpoinv.fptr_desc = NULL;
#else
    test_zpoinv.fptr_desc = testing_zpoinv_desc;
#endif
    test_zpoinv.fptr_std  = testing_zpoinv_std;
    test_zpoinv.next   = NULL;

    testing_register( &test_zpoinv );
}
