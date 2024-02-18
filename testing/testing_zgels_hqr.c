/**
 *
 * @file testing_zgels_hqr.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgels_hqr testing
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @author Alycia Lisito
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zgels_hqr( cham_trans_t trans, int M, int N, int NRHS )
{
    cham_fixdbl_t flops = 0.;
    (void)trans;
    (void)M;
    (void)N;
    (void)NRHS;
    return flops;
}

int
testing_zgels_hqr_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_nb( args );
    int          ib     = run_arg_get_ib( args );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          maxMN  = chameleon_max( M, N );
    int          NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int          LDA    = run_arg_get_int( args, "LDA", M );
    int          LDB    = run_arg_get_int( args, "LDB", maxMN );
    int          qr_a   = run_arg_get_int( args, "qra", -1 );
    int          qr_p   = run_arg_get_int( args, "qrp", -1 );
    int          llvl   = run_arg_get_int( args, "llvl", -1 );
    int          hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int          domino = run_arg_get_int( args, "domino", -1 );
    int          seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedB  = run_arg_get_int( args, "seedB", testing_ialea() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t    *descA, *descX, *descTS, *descTT;
    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, maxMN, NRHS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, ( M >= N ) ? LIBHQR_QR : LIBHQR_LQ, &matrix,
                     llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Computes the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgels_param_Tile_Async( &qrtree, trans, descA, descTS, descTT, descX,
                                                 test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descTS, test_data.sequence );
        CHAMELEON_Desc_Flush( descTT, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgels_param_Tile( &qrtree, trans, descA, descTS, descTT, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgels_hqr( trans, M, N, NRHS ) );

    if ( check ) {
        CHAM_desc_t *descA0, *descB;
        CHAM_desc_t *subX, *subB;

        CHAMELEON_Desc_Create(
            &descA0, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
        CHAMELEON_Desc_Create(
            &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, maxMN, NRHS, P, Q );

        CHAMELEON_zplrnt_Tile( descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        if ( trans == ChamNoTrans ) {
            subX = CHAMELEON_Desc_SubMatrix( descX, 0, 0, N, NRHS );
            subB = CHAMELEON_Desc_SubMatrix( descB, 0, 0, M, NRHS );
        }
        else {
            subX = CHAMELEON_Desc_SubMatrix( descX, 0, 0, M, NRHS );
            subB = CHAMELEON_Desc_SubMatrix( descB, 0, 0, N, NRHS );
        }

        /* Check the factorization and the residual */
        hres = check_zgels( args, trans, descA0, subX, subB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );

        free( subB );
        free( subX );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    CHAMELEON_Desc_Destroy( &descX );
    libhqr_finalize( &qrtree );

    return hres;
}

int
testing_zgels_hqr_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          nb     = run_arg_get_nb( args );
    int          ib     = run_arg_get_ib( args );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          maxMN  = chameleon_max( M, N );
    int          NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int          LDA    = run_arg_get_int( args, "LDA", M );
    int          LDB    = run_arg_get_int( args, "LDB", maxMN );
    int          qr_a   = run_arg_get_int( args, "qra", -1 );
    int          qr_p   = run_arg_get_int( args, "qrp", -1 );
    int          llvl   = run_arg_get_int( args, "llvl", -1 );
    int          hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int          domino = run_arg_get_int( args, "domino", -1 );
    int          seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedB  = run_arg_get_int( args, "seedB", testing_ialea() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *X;
    CHAM_desc_t           *descTS, *descTT;
    libhqr_tree_t          qrtree;
    libhqr_matrix_t        matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    A = malloc( LDA*N*   sizeof(CHAMELEON_Complex64_t) );
    X = malloc( LDB*NRHS*sizeof(CHAMELEON_Complex64_t) );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, ( M >= N ) ? LIBHQR_QR : LIBHQR_LQ, &matrix,
                     llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M,     N,    A, LDA, seedA );
    CHAMELEON_zplrnt( maxMN, NRHS, X, LDB, seedB );

    /* Computes the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zgels_param( &qrtree, trans, M, N, NRHS, A, LDA, descTS, descTT, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgels_hqr( trans, M, N, NRHS ) );

    if ( check ) {
        CHAMELEON_Complex64_t *A0, *B;

        A0 = malloc( LDA*N*   sizeof(CHAMELEON_Complex64_t) );
        B  = malloc( LDB*NRHS*sizeof(CHAMELEON_Complex64_t) );

        CHAMELEON_zplrnt( M,     N,    A0, LDA, seedA );
        CHAMELEON_zplrnt( maxMN, NRHS, B,  LDB, seedB );

        /* Check the factorization and the residual */
        hres = check_zgels_std( args, trans, M, N, NRHS, A0, LDA, X, LDB, B, LDB );

        free( A0 );
        free( B );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    free( X );
    libhqr_finalize( &qrtree );

    return hres;
}

testing_t   test_zgels_hqr;
const char *zgels_hqr_params[] = { "mtxfmt", "nb",   "ib",     "trans", "m",     "n",
                                   "k",      "lda",  "ldb",    "qra",   "qra",   "qrp",
                                   "llvl",   "hlvl", "domino", "seedA", "seedB", NULL };
const char *zgels_hqr_output[] = { NULL };
const char *zgels_hqr_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgels_hqr_init( void ) __attribute__( ( constructor ) );
void
testing_zgels_hqr_init( void )
{
    test_zgels_hqr.name = "zgels_hqr";
    test_zgels_hqr.helper =
        "Linear least squares with general matrix using hierarchical reduction trees";
    test_zgels_hqr.params = zgels_hqr_params;
    test_zgels_hqr.output = zgels_hqr_output;
    test_zgels_hqr.outchk = zgels_hqr_outchk;
    test_zgels_hqr.fptr_desc = testing_zgels_hqr_desc;
    test_zgels_hqr.fptr_std  = testing_zgels_hqr_std;
    test_zgels_hqr.next   = NULL;

    testing_register( &test_zgels_hqr );
}
