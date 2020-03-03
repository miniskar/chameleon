/**
 *
 * @file testing_zlange.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange testing
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include "flops.h"

static cham_fixdbl_t
flops_zlange( cham_normtype_t ntype, int M, int N )
{
    cham_fixdbl_t flops   = 0.;
    double coefabs = 1.;
#if defined( PRECISION_z ) || defined( PRECISION_c )
    coefabs = 3.;
#endif

    switch ( ntype ) {
        case ChamMaxNorm:
            flops = coefabs * M * N;
            break;
        case ChamOneNorm:
            flops = coefabs * M * N + M * ( N - 1 );
            break;
        case ChamInfNorm:
            flops = coefabs * M * N + N * ( M - 1 );
            break;
        case ChamFrobeniusNorm:
            flops = ( coefabs + 1. ) * M * N;
            break;
        default:;
    }
    return flops;
}

int
testing_zlange( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    double       norm;
    CHAM_desc_t *descA;

    /* Reads arguments */
    int             nb        = run_arg_get_int( args, "nb", 320 );
    int             P         = parameters_getvalue_int( "P" );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             M         = run_arg_get_int( args, "M", N );
    int             LDA       = run_arg_get_int( args, "LDA", M );
    int             seedA     = run_arg_get_int( args, "seedA", random() );
    int             Q         = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zlange( norm_type, M, N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the norm */
    START_TIMING( t );
    norm = CHAMELEON_zlange_Tile( norm_type, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( norm >= 0. ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        hres = check_znorm( args, ChamGeneral, norm_type, ChamUpperLower, ChamNonUnit, norm, descA );
    }

    CHAMELEON_Desc_Destroy( &descA );

    run_id++;
    return hres;
}

testing_t   test_zlange;
const char *zlange_params[] = { "nb", "norm", "m", "n", "lda", "seedA", NULL };
const char *zlange_output[] = { NULL };
const char *zlange_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlange_init( void ) __attribute__( ( constructor ) );
void
testing_zlange_init( void )
{
    test_zlange.name        = "zlange";
    test_zlange.helper      = "General matrix norm";
    test_zlange.params      = zlange_params;
    test_zlange.output      = zlange_output;
    test_zlange.outchk      = zlange_outchk;
    test_zlange.params_list = "nb;P;norm;m;n;lda;seedA";
    test_zlange.fptr        = testing_zlange;
    test_zlange.next        = NULL;

    testing_register( &test_zlange );
}
