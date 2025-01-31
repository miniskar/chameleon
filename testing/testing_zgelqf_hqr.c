/**
 *
 * @file testing_zgelqf_hqr.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf_param testing
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

int
testing_zgelqf_hqr_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_nb( args );
    int      ib     = run_arg_get_ib( args );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      qr_a   = run_arg_get_int( args, "qra", -1 );
    int      qr_p   = run_arg_get_int( args, "qrp", -1 );
    int      llvl   = run_arg_get_int( args, "llvl", -1 );
    int      hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int      domino = run_arg_get_int( args, "domino", -1 );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t    *descA, *descTS, *descTT;
    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
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
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgelqf_param_Tile_Async( &qrtree, descA, descTS, descTT,
                                                  test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descTS, test_data.sequence );
        CHAMELEON_Desc_Flush( descTT, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgelqf_param_Tile( &qrtree, descA, descTS, descTT );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgelqf( M, N ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descQ;
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, N, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zunglq_param_Tile( &qrtree, descA, descTS, descTT, descQ );

        hres += check_zgelqf( args, descA0, descA, descQ );
        hres += check_zortho( args, descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    libhqr_finalize( &qrtree );

    return hres;
}

int
testing_zgelqf_hqr_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int nb     = run_arg_get_nb( args );
    int ib     = run_arg_get_ib( args );
    int P      = parameters_getvalue_int( "P" );
    int N      = run_arg_get_int( args, "N", 1000 );
    int M      = run_arg_get_int( args, "M", N );
    int LDA    = run_arg_get_int( args, "LDA", M );
    int qr_a   = run_arg_get_int( args, "qra", -1 );
    int qr_p   = run_arg_get_int( args, "qrp", -1 );
    int llvl   = run_arg_get_int( args, "llvl", -1 );
    int hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int domino = run_arg_get_int( args, "domino", -1 );
    int seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int Q      = parameters_compute_q( P );
    int K      = chameleon_min( M, N );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;
    CHAM_desc_t           *descTS, *descTT;
    libhqr_tree_t          qrtree;
    libhqr_matrix_t        matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_LQ, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zgelqf_param( &qrtree, M, N, A, LDA, descTS, descTT );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgelqf( M, N ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAMELEON_Complex64_t *Qlap = malloc( N*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_Complex64_t *A0   = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

        CHAMELEON_zplrnt( M, N, A0, LDA, seedA );

        CHAMELEON_zunglq_param( &qrtree, N, N, K, A, LDA, descTS, descTT, Qlap, N );

        hres += check_zgelqf_std( args, N, N, K, A0, A, LDA, Qlap, N );
        hres += check_zortho_std( args, N, N, Qlap, N );

        free( A0 );
        free( Qlap );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    libhqr_finalize( &qrtree );

    return hres;
}

testing_t   test_zgelqf_hqr;
const char *zgelqf_hqr_params[] = { "mtxfmt", "nb",   "ib",   "m",      "n",     "lda", "qra",
                                    "qrp",    "llvl", "hlvl", "domino", "seedA", NULL };
const char *zgelqf_hqr_output[] = { NULL };
const char *zgelqf_hqr_outchk[] = { "||A||", "||I-QQ'||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgelqf_hqr_init( void ) __attribute__( ( constructor ) );
void
testing_zgelqf_hqr_init( void )
{
    test_zgelqf_hqr.name   = "zgelqf_hqr";
    test_zgelqf_hqr.helper = "General LQ factorization with hierachical reduction trees";
    test_zgelqf_hqr.params = zgelqf_hqr_params;
    test_zgelqf_hqr.output = zgelqf_hqr_output;
    test_zgelqf_hqr.outchk = zgelqf_hqr_outchk;
    test_zgelqf_hqr.fptr_desc = testing_zgelqf_hqr_desc;
    test_zgelqf_hqr.fptr_std  = testing_zgelqf_hqr_std;
    test_zgelqf_hqr.next   = NULL;

    testing_register( &test_zgelqf_hqr );
}
