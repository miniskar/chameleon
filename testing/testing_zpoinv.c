/**
 *
 * @file testing_zpoinv.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpoinv testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-02
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zpoinv( int N )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotri( N );
    return flops;
}

int
testing_zpoinv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         async  = parameters_getvalue_int( "async" );
    intptr_t    mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int         nb     = run_arg_get_int( args, "nb", 320 );
    int         P      = parameters_getvalue_int( "P" );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         LDA    = run_arg_get_int( args, "LDA", N );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

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
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_ztrtri( args, ChamHermitian, uplo, ChamNonUnit, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zpoinv;
const char *zpoinv_params[] = { "mtxfmt", "nb", "uplo", "n", "lda", "seedA", NULL };
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
    test_zpoinv.fptr_desc = testing_zpoinv_desc;
    test_zpoinv.fptr_std  = NULL;
    test_zpoinv.next   = NULL;

    testing_register( &test_zpoinv );
}
