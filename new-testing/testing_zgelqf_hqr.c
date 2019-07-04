/**
 *
 * @file testing_zgelqf_hqr.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf_param testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-09-10
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zgelqf_hqr( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descTS, *descTT;

    /* Reads arguments */
    int    nb     = run_arg_get_int( args, "nb", 320 );
    int    ib     = run_arg_get_int( args, "ib", 48 );
    int    P      = parameters_getvalue_int( "P" );
    int    N      = run_arg_get_int( args, "N", 1000 );
    int    M      = run_arg_get_int( args, "M", N );
    int    LDA    = run_arg_get_int( args, "LDA", M );
    int    qr_a   = run_arg_get_int( args, "qra", -1 );
    int    qr_p   = run_arg_get_int( args, "qrp", -1 );
    int    llvl   = run_arg_get_int( args, "llvl", -1 );
    int    hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int    domino = run_arg_get_int( args, "domino", -1 );
    int    seedA  = run_arg_get_int( args, "seedA", random() );
    int    Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgelqf( M, N );

    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_LQ, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the solution */
    START_TIMING( t );
    hres = CHAMELEON_zgelqf_param_Tile( &qrtree, descA, descTS, descTT );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", gflops );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descQ;
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );

        CHAMELEON_Desc_Create(
            &descQ, NULL, ChamComplexDouble, nb, nb, nb * nb, N, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zunglq_param_Tile( &qrtree, descA, descTS, descTT, descQ );

        hres += check_zgelqf( descA0, descA, descQ );
        hres += check_zortho( descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    libhqr_finalize( &qrtree );

    run_id++;
    return hres;
}

testing_t   test_zgelqf_hqr;
const char *zgelqf_hqr_params[] = { "nb",  "ib",   "m",    "n",      "lda",   "qra",
                                    "qrp", "llvl", "hlvl", "domino", "seedA", NULL };
const char *zgelqf_hqr_output[] = { NULL };
const char *zgelqf_hqr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgelqf_hqr_init( void ) __attribute__( ( constructor ) );
void
testing_zgelqf_hqr_init( void )
{
    test_zgelqf_hqr.name        = "zgelqf_hqr";
    test_zgelqf_hqr.helper      = "zgelqf_hqr";
    test_zgelqf_hqr.params      = zgelqf_hqr_params;
    test_zgelqf_hqr.output      = zgelqf_hqr_output;
    test_zgelqf_hqr.outchk      = zgelqf_hqr_outchk;
    test_zgelqf_hqr.params_list = "nb;ib;P;m;n;lda;qra;qrp;llvl;hlvl;domino;seedA";
    test_zgelqf_hqr.fptr        = testing_zgelqf_hqr;
    test_zgelqf_hqr.next        = NULL;

    testing_register( &test_zgelqf_hqr );
}
