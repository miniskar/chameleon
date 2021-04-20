/**
 *
 * @file testing_zlauum.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum testing
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
flops_zlauum( int N )
{
    cham_fixdbl_t flops = flops_zpotri( N ) - flops_ztrtri( N );
    return flops;
}

int
testing_zlauum( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA;

    /* Reads arguments */
    intptr_t    mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int         nb     = run_arg_get_int( args, "nb", 320 );
    int         P      = parameters_getvalue_int( "P" );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         LDA    = run_arg_get_int( args, "LDA", N );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zlauum( N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

    /* Initialises the matrices with the same values */
    CHAMELEON_zplghe_Tile( 0., uplo, descA, seedA );

    /* Calculates the matrix product */
    START_TIMING( t );
    hres = CHAMELEON_zlauum_Tile( uplo, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplghe_Tile( 0., uplo, descA0, seedA );

        hres += check_zlauum( args, uplo, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zlauum;
const char *zlauum_params[] = { "mtxfmt", "nb","uplo", "n", "lda", "seedA", NULL };
const char *zlauum_output[] = { NULL };
const char *zlauum_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlauum_init( void ) __attribute__( ( constructor ) );
void
testing_zlauum_init( void )
{
    test_zlauum.name        = "zlauum";
    test_zlauum.helper      = "Trianguilar in-place matrix-matrix computation for Cholesky inversion";
    test_zlauum.params      = zlauum_params;
    test_zlauum.output      = zlauum_output;
    test_zlauum.outchk      = zlauum_outchk;
    test_zlauum.fptr        = testing_zlauum;
    test_zlauum.next        = NULL;

    testing_register( &test_zlauum );
}
