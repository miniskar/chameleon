/**
 *
 * @file testing_zlantr.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlantr testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zlantr( cham_normtype_t ntype, cham_uplo_t uplo, int M, int N )
{
    cham_fixdbl_t flops = 0.;
    double coefabs = 1.;
#if defined( PRECISION_z ) || defined( PRECISION_c )
    coefabs = 3.;
#endif

    switch ( uplo ) {
        case ChamUpper:
            if ( N > M ) {
                flops = ( M * ( M + 1 ) / 2 ) + M * ( N - M );
            }
            else {
                flops = N * ( N + 1 ) / 2;
            }
            break;
        case ChamLower:
            if ( M > N ) {
                flops = ( N * ( N + 1 ) / 2 ) + N * ( M - N );
            }
            else {
                flops = M * ( M + 1 ) / 2;
            }
            break;
        case ChamUpperLower:
        default:
            flops = M * N;
    }
    flops *= coefabs;

    switch ( ntype ) {
        case ChamOneNorm:
            flops += N;
            break;
        case ChamInfNorm:
            flops += M;
            break;
        case ChamMaxNorm:
        case ChamFrobeniusNorm:
        default:;
    }
    return flops;
}

int
testing_zlantr( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    double       norm;
    CHAM_desc_t *descA;

    /* Reads arguments */
    intptr_t        mtxfmt    = parameters_getvalue_int( "mtxfmt" );
    int             nb        = run_arg_get_int( args, "nb", 320 );
    int             P         = parameters_getvalue_int( "P" );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t     uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t     diag      = run_arg_get_diag( args, "diag", ChamNonUnit );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             M         = run_arg_get_int( args, "M", N );
    int             LDA       = run_arg_get_int( args, "LDA", M );
    int             seedA     = run_arg_get_int( args, "seedA", random() );
    int             Q         = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zlantr( norm_type, uplo, M, N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the norm */
    START_TIMING( t );
    norm = CHAMELEON_zlantr_Tile( norm_type, uplo, diag, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( norm >= 0. ) ? gflops : -1. );

    /* Checks the solution */
    if ( check ) {
        hres = check_znorm( args, ChamTriangular, norm_type, uplo, diag, norm, descA );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zlantr;
const char *zlantr_params[] = { "mtxfmt", "nb","norm", "uplo", "diag", "m", "n", "lda", "seedA", NULL };
const char *zlantr_output[] = { NULL };
const char *zlantr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlantr_init( void ) __attribute__( ( constructor ) );
void
testing_zlantr_init( void )
{
    test_zlantr.name        = "zlantr";
    test_zlantr.helper      = "Triangular matrix norm";
    test_zlantr.params      = zlantr_params;
    test_zlantr.output      = zlantr_output;
    test_zlantr.outchk      = zlantr_outchk;
    test_zlantr.fptr        = testing_zlantr;
    test_zlantr.next        = NULL;

    testing_register( &test_zlantr );
}
