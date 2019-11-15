/**
 *
 * @file testing_zunglq.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-09-09
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zunglq( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descT, *descQ;

    /* Reads arguments */
    int    nb    = run_arg_get_int( args, "nb", 320 );
    int    ib    = run_arg_get_int( args, "ib", 48 );
    int    P     = parameters_getvalue_int( "P" );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    K     = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    RH    = run_arg_get_int( args, "qra", 0 );
    int    seedA = run_arg_get_int( args, "seedA", random() );
    int    Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zunglq( M, N, K );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( K > chameleon_min( M, N ) ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameters for unglq (K > min(M,N))\n" );
        }
        return -1;
    }

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, K, N, P, Q );
    CHAMELEON_Desc_Create(
        &descQ, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    hres = CHAMELEON_zgelqf_Tile( descA, descT );

    /* Calculates the solution */
    START_TIMING( t );
    CHAMELEON_zunglq_Tile( descA, descT, descQ );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", gflops );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        hres += check_zortho( args, descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descQ );

    run_id++;
    return hres;
}

testing_t   test_zunglq;
const char *zunglq_params[] = { "nb", "ib", "m", "n", "k", "lda", "qra", "seedA", NULL };
const char *zunglq_output[] = { NULL };
const char *zunglq_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunglq_init( void ) __attribute__( ( constructor ) );
void
testing_zunglq_init( void )
{
    test_zunglq.name        = "zunglq";
    test_zunglq.helper      = "zunglq";
    test_zunglq.params      = zunglq_params;
    test_zunglq.output      = zunglq_output;
    test_zunglq.outchk      = zunglq_outchk;
    test_zunglq.params_list = "nb;ib;P;m;n;k;lda;rh;seedA";
    test_zunglq.fptr        = testing_zunglq;
    test_zunglq.next        = NULL;

    testing_register( &test_zunglq );
}
