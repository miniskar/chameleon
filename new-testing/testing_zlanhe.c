/**
 *
 * @file testing_zlanhe.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlanhe testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-07-17
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

static cham_fixdbl_t
flops_zlanhe( cham_normtype_t ntype, cham_uplo_t uplo, int N )
{
    cham_fixdbl_t flops   = 0.;
    double coefabs = 1.;
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
    static int   run_id = 0;
    int          hres   = 0;
    double       norm;
    CHAM_desc_t *descA;

    /* Reads arguments */
    int             nb        = run_arg_get_int( args, "nb", 320 );
    int             P         = parameters_getvalue_int( "P" );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t     uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             LDA       = run_arg_get_int( args, "LDA", N );
    int             seedA     = run_arg_get_int( args, "seedA", random() );
    double          bump      = testing_dalea();
    int             Q         = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zlanhe( norm_type, uplo, N );

    bump = run_arg_get_double( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( bump, uplo, descA, seedA );

    /* Calculates the norm */
    START_TIMING( t );
    norm = CHAMELEON_zlanhe_Tile( norm_type, uplo, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        hres = check_znorm( args, ChamHermitian, norm_type, uplo, ChamNonUnit, norm, descA );
    }

    CHAMELEON_Desc_Destroy( &descA );

    run_id++;
    return hres;
}

testing_t   test_zlanhe;
const char *zlanhe_params[] = { "nb", "norm", "uplo", "n", "lda", "seedA", "bump", NULL };
const char *zlanhe_output[] = { NULL };
const char *zlanhe_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlanhe_init( void ) __attribute__( ( constructor ) );
void
testing_zlanhe_init( void )
{
    test_zlanhe.name        = "zlanhe";
    test_zlanhe.helper      = "zlanhe";
    test_zlanhe.params      = zlanhe_params;
    test_zlanhe.output      = zlanhe_output;
    test_zlanhe.outchk      = zlanhe_outchk;
    test_zlanhe.params_list = "nb;P;norm;uplo;n;lda;seedA;bump";
    test_zlanhe.fptr        = testing_zlanhe;
    test_zlanhe.next        = NULL;

    testing_register( &test_zlanhe );
}
