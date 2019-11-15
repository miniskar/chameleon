/**
 *
 * @file testing_zgeqrf.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf testing
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
testing_zgeqrf( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descT;

    /* Reads arguments */
    int    nb    = run_arg_get_int( args, "nb", 320 );
    int    ib    = run_arg_get_int( args, "ib", 48 );
    int    P     = parameters_getvalue_int( "P" );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    RH    = run_arg_get_int( args, "qra", 4 );
    int    seedA = run_arg_get_int( args, "seedA", random() );
    int    Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgeqrf( M, N );

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
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the solution */
    START_TIMING( t );
    hres = CHAMELEON_zgeqrf_Tile( descA, descT );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", gflops );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descQ;
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );

        CHAMELEON_Desc_Create(
            &descQ, NULL, ChamComplexDouble, nb, nb, nb * nb, M, M, 0, 0, M, M, P, Q );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zungqr_Tile( descA, descT, descQ );

        hres += check_zgeqrf( descA0, descA, descQ );
        hres += check_zortho( descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );

    run_id++;
    return hres;
}

testing_t   test_zgeqrf;
const char *zgeqrf_params[] = { "nb", "ib", "m", "n", "lda", "qra", "seedA", NULL };
const char *zgeqrf_output[] = { NULL };
const char *zgeqrf_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgeqrf_init( void ) __attribute__( ( constructor ) );
void
testing_zgeqrf_init( void )
{
    test_zgeqrf.name        = "zgeqrf";
    test_zgeqrf.helper      = "zgeqrf";
    test_zgeqrf.params      = zgeqrf_params;
    test_zgeqrf.output      = zgeqrf_output;
    test_zgeqrf.outchk      = zgeqrf_outchk;
    test_zgeqrf.params_list = "nb;ib;P;m;n;lda;rh;seedA";
    test_zgeqrf.fptr        = testing_zgeqrf;
    test_zgeqrf.next        = NULL;

    testing_register( &test_zgeqrf );
}
