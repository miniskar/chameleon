/**
 *
 * @file testing_zlanhe.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlanhe testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zlanhe( cham_normtype_t ntype, int N )
{
    cham_fixdbl_t flops   = 0.;
    double        coefabs = 1.;
#if defined( PRECISION_z ) || defined( PRECISION_c )
    coefabs = 3.;
#endif

    switch ( ntype ) {
        case ChamMaxNorm:
            flops = coefabs * ( N * ( N + 1 ) ) / 2.;
            break;
        case ChamOneNorm:
        case ChamInfNorm:
            flops = coefabs * ( N * ( N + 1 ) ) / 2. + N * ( N - 1 );
            break;
        case ChamFrobeniusNorm:
            flops = ( coefabs + 1. ) * ( N * ( N + 1 ) ) / 2.;
            break;
        default:;
    }
    return flops;
}

int
testing_zlanhe( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int             async     = parameters_getvalue_int( "async" );
    intptr_t        mtxfmt    = parameters_getvalue_int( "mtxfmt" );
    int             nb        = run_arg_get_int( args, "nb", 320 );
    int             P         = parameters_getvalue_int( "P" );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t     uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             LDA       = run_arg_get_int( args, "LDA", N );
    int             seedA     = run_arg_get_int( args, "seedA", random() );
    double          bump      = testing_dalea();
    int             Q         = parameters_compute_q( P );

    /* Descriptors */
    double       norm;
    CHAM_desc_t *descA;

    bump = run_arg_get_double( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( bump, uplo, descA, seedA );

    /* Calculates the norm */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlanhe_Tile_Async( norm_type, uplo, descA, &norm,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        norm = CHAMELEON_zlanhe_Tile( norm_type, uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlanhe( norm_type, N ) );

    /* Checks the solution */
    if ( check ) {
        hres = check_znorm( args, ChamHermitian, norm_type, uplo, ChamNonUnit, norm, descA );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zlanhe;
const char *zlanhe_params[] = { "mtxfmt", "nb", "norm", "uplo", "n", "lda", "seedA", "bump", NULL };
const char *zlanhe_output[] = { NULL };
const char *zlanhe_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlanhe_init( void ) __attribute__( ( constructor ) );
void
testing_zlanhe_init( void )
{
    test_zlanhe.name   = "zlanhe";
    test_zlanhe.helper = "Hermitian matrix norm";
    test_zlanhe.params = zlanhe_params;
    test_zlanhe.output = zlanhe_output;
    test_zlanhe.outchk = zlanhe_outchk;
    test_zlanhe.fptr   = testing_zlanhe;
    test_zlanhe.next   = NULL;

    testing_register( &test_zlanhe );
}
