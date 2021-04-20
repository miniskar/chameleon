/**
 *
 * @file testing_zgelqf.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zgelqf( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA, *descT;

    /* Reads arguments */
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_int( args, "nb", 320 );
    int      ib     = run_arg_get_int( args, "ib", 48 );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      RH     = run_arg_get_int( args, "qra", 4 );
    int      seedA  = run_arg_get_int( args, "seedA", random() );
    int      Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgelqf( M, N );

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
    START_TIMING( t );
    hres = CHAMELEON_zgelqf_Tile( descA, descT );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descQ;
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, N, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zunglq_Tile( descA, descT, descQ );

        hres += check_zgelqf( args, descA0, descA, descQ );
        hres += check_zortho( args, descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );

    return hres;
}

testing_t   test_zgelqf;
const char *zgelqf_params[] = { "mtxfmt", "nb","ib", "m", "n", "lda", "qra", "seedA", NULL };
const char *zgelqf_output[] = { NULL };
const char *zgelqf_outchk[] = { "||A||", "||I-QQ'||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgelqf_init( void ) __attribute__( ( constructor ) );
void
testing_zgelqf_init( void )
{
    test_zgelqf.name        = "zgelqf";
    test_zgelqf.helper      = "General LQ factorization";
    test_zgelqf.params      = zgelqf_params;
    test_zgelqf.output      = zgelqf_output;
    test_zgelqf.outchk      = zgelqf_outchk;
    test_zgelqf.fptr        = testing_zgelqf;
    test_zgelqf.next        = NULL;

    testing_register( &test_zgelqf );
}
