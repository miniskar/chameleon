/**
 *
 * @file testing_zunmqr_hqr.c
 *
 * @copyright 2019-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr_hqr testing
 *
 * @version 1.3.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2023-07-05
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#include <assert.h>

int
testing_zunmqr_hqr_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_nb( args );
    int          ib     = run_arg_get_ib( args );
    int          P      = parameters_getvalue_int( "P" );
    cham_side_t  side   = run_arg_get_side( args, "side", ChamLeft );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int          LDA    = run_arg_get_int( args, "LDA", ( side == ChamLeft ) ? M : N );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    int          qr_a   = run_arg_get_int( args, "qra", -1 );
    int          qr_p   = run_arg_get_int( args, "qrp", -1 );
    int          llvl   = run_arg_get_int( args, "llvl", -1 );
    int          hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int          domino = run_arg_get_int( args, "domino", -1 );
    int          seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedC  = run_arg_get_int( args, "seedC", testing_ialea() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    int             Am;
    CHAM_desc_t    *descA, *descTS, *descTT, *descC;
    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Calculates the dimensions according to the transposition and the side */
    Am = ( side == ChamLeft ) ? M : N;

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, Am, K, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( Am, K, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( Am, K, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Computes the factorization */
    hres = CHAMELEON_zgeqrf_param_Tile( &qrtree, descA, descTS, descTT );
    assert( hres == 0 );

    /* Computes unmqr_hqr */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zunmqr_param_Tile_Async( &qrtree, side, trans, descA, descTS, descTT, descC,
                                                   test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descTS, test_data.sequence );
        CHAMELEON_Desc_Flush( descTT, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zunmqr_param_Tile( &qrtree, side, trans, descA, descTS, descTT, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zunmqr( side, M, N, K ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descC0 = CHAMELEON_Desc_Copy( descC, CHAMELEON_MAT_ALLOC_TILE );
        CHAM_desc_t *descQ;

        CHAMELEON_zplrnt_Tile( descC0, seedC );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, Am, Am, 0, 0, Am, Am, P, Q );
        CHAMELEON_zungqr_param_Tile( &qrtree, descA, descTS, descTT, descQ );

        hres += check_zqc( args, side, trans, descC0, descQ, descC );

        CHAMELEON_Desc_Destroy( &descC0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    CHAMELEON_Desc_Destroy( &descC );
    libhqr_finalize( &qrtree );

    return hres;
}

int
testing_zunmqr_hqr_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          nb     = run_arg_get_nb( args );
    int          ib     = run_arg_get_ib( args );
    int          P      = parameters_getvalue_int( "P" );
    cham_side_t  side   = run_arg_get_side( args, "side", ChamLeft );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int          LDA    = run_arg_get_int( args, "LDA", ( side == ChamLeft ) ? M : N );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    int          qr_a   = run_arg_get_int( args, "qra", -1 );
    int          qr_p   = run_arg_get_int( args, "qrp", -1 );
    int          llvl   = run_arg_get_int( args, "llvl", -1 );
    int          hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int          domino = run_arg_get_int( args, "domino", -1 );
    int          seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedC  = run_arg_get_int( args, "seedC", testing_ialea() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    int                    Am;
    CHAMELEON_Complex64_t *A, *C;
    CHAM_desc_t           *descTS, *descTT;
    libhqr_tree_t          qrtree;
    libhqr_matrix_t        matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Calculates the dimensions according to the transposition and the side */
    Am = ( side == ChamLeft ) ? M : N;

    /* Creates the matrices */
    A = malloc( LDA*K*sizeof(CHAMELEON_Complex64_t) );
    C = malloc( LDC*N*sizeof(CHAMELEON_Complex64_t) );
    CHAMELEON_Alloc_Workspace_zgels( Am, K, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( Am, K, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( Am, K, A, LDA, seedA );
    CHAMELEON_zplrnt( M,  N, C, LDC, seedC );

    /* Computes the factorization */
    hres = CHAMELEON_zgeqrf_param( &qrtree, Am, K, A, LDA, descTS, descTT );

    /* Computes unmqr_hqr */
    testing_start( &test_data );
    hres += CHAMELEON_zunmqr_param( &qrtree, side, trans, M, N, K, A, LDA, descTS, descTT, C, LDC );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zunmqr( side, M, N, K ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAMELEON_Complex64_t *C0   = malloc( LDC*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_Complex64_t *Qlap = malloc( Am*Am*sizeof(CHAMELEON_Complex64_t) );

        CHAMELEON_zplrnt( M, N, C0, LDC, seedC );
        CHAMELEON_zungqr_param( &qrtree, Am, Am, K, A, LDA, descTS, descTT, Qlap, Am );

        hres += check_zqc_std( args, side, trans, M, N, C0, C, LDC, Qlap, Am );

        free( C0 );
        free( Qlap );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    free( C );
    libhqr_finalize( &qrtree );

    return hres;
}

testing_t   test_zunmqr_hqr;
const char *zunmqr_hqr_params[] = { "mtxfmt", "nb",   "ib",     "side",  "trans", "m",
                                    "n",      "k",    "lda",    "ldc",   "qra",   "qrp",
                                    "llvl",   "hlvl", "domino", "seedA", "seedC", NULL };
const char *zunmqr_hqr_output[] = { NULL };
const char *zunmqr_hqr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunmqr_hqr_init( void ) __attribute__( ( constructor ) );
void
testing_zunmqr_hqr_init( void )
{
    test_zunmqr_hqr.name   = "zunmqr_hqr";
    test_zunmqr_hqr.helper = "Q application with hierarchical reduction trees (QR)";
    test_zunmqr_hqr.params = zunmqr_hqr_params;
    test_zunmqr_hqr.output = zunmqr_hqr_output;
    test_zunmqr_hqr.outchk = zunmqr_hqr_outchk;
    test_zunmqr_hqr.fptr_desc = testing_zunmqr_hqr_desc;
    test_zunmqr_hqr.fptr_std  = testing_zunmqr_hqr_std;
    test_zunmqr_hqr.next   = NULL;

    testing_register( &test_zunmqr_hqr );
}
