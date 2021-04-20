/**
 *
 * @file testing_zlacpy.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy testing
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
flops_zlacpy( cham_uplo_t uplo, int M, int N )
{
    cham_fixdbl_t flops;

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
    flops *= sizeof( CHAMELEON_Complex64_t );

    return flops;
}

int
testing_zlacpy( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA, *descB;

    /* Reads arguments */
    intptr_t    mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int         nb     = run_arg_get_int( args, "nb", 320 );
    int         P      = parameters_getvalue_int( "P" );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         M      = run_arg_get_int( args, "M", N );
    int         LDA    = run_arg_get_int( args, "LDA", M );
    int         LDB    = run_arg_get_int( args, "LDB", M );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zlacpy( uplo, M, N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates two different matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );

    /* Fills each matrix with different random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    /* We use seedA + 1, just to create a variation in B */
    CHAMELEON_zplrnt_Tile( descB, seedA + 1 );

    /* Makes a copy of descA to descB */
    START_TIMING( t );
    hres = CHAMELEON_zlacpy_Tile( uplo, descA, descB );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks their differences */
    if ( check ) {
        hres += check_zmatrices( args, uplo, descA, descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );

    return hres;
}

testing_t   test_zlacpy;
const char *zlacpy_params[] = { "mtxfmt", "nb","uplo", "m", "n", "lda", "ldb", "seedA", NULL };
const char *zlacpy_output[] = { NULL };
const char *zlacpy_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlacpy_init( void ) __attribute__( ( constructor ) );
void
testing_zlacpy_init( void )
{
    test_zlacpy.name        = "zlacpy";
    test_zlacpy.helper      = "General matrix copy";
    test_zlacpy.params      = zlacpy_params;
    test_zlacpy.output      = zlacpy_output;
    test_zlacpy.outchk      = zlacpy_outchk;
    test_zlacpy.fptr        = testing_zlacpy;
    test_zlacpy.next        = NULL;

    testing_register( &test_zlacpy );
}
