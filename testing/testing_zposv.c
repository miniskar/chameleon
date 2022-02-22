/**
 *
 * @file testing_zposv.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zposv testing
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zposv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotrs( N, NRHS );
    return flops;
}

int
testing_zposv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         async  = parameters_getvalue_int( "async" );
    intptr_t    mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int         nb     = run_arg_get_int( args, "nb", 320 );
    int         P      = parameters_getvalue_int( "P" );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int         LDA    = run_arg_get_int( args, "LDA", N );
    int         LDB    = run_arg_get_int( args, "LDB", N );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         seedB  = run_arg_get_int( args, "seedB", random() );
    int         Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descX;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zposv_Tile_Async( uplo, descA, descX,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zposv_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zposv( N, NRHS ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_zxxtrf( args, ChamHermitian, uplo, descA0, descA );

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, NULL );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );
        hres += check_zsolve( args, ChamHermitian, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    return hres;
}

int
testing_zposv_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         nb    = run_arg_get_int( args, "nb", 320 );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         LDB   = run_arg_get_int( args, "LDB", N );
    int         seedA = run_arg_get_int( args, "seedA", random() );
    int         seedB = run_arg_get_int( args, "seedB", random() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *X;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( LDA*N   *sizeof(CHAMELEON_Complex64_t) );
    X = malloc( LDB*NRHS*sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zposv( uplo, N, NRHS, A, LDA, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zposv( N, NRHS ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAMELEON_Complex64_t *A0 = malloc( LDA*N*   sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_Complex64_t *B  = malloc( LDB*NRHS*sizeof(CHAMELEON_Complex64_t) );

        /* Check the factorization */
        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_zxxtrf_std( args, ChamHermitian, uplo, N, N, A0, A, LDA );

        /* Check the solve */
        CHAMELEON_zplrnt( N, NRHS, B, LDB, seedB );

        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );
        hres += check_zsolve_std( args, ChamHermitian, ChamNoTrans, uplo, N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }

    free( A );
    free( X );

    return hres;
}

testing_t   test_zposv;
const char *zposv_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                               "lda",    "ldb", "seedA", "seedB", NULL };
const char *zposv_output[] = { NULL };
const char *zposv_outchk[] = { "||A||", "||A-fact(A)||", "||X||", "||B||", "||Ax-b||", "||Ax-b||/N/eps/(||A||||x||+||b||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zposv_init( void ) __attribute__( ( constructor ) );
void
testing_zposv_init( void )
{
    test_zposv.name   = "zposv";
    test_zposv.helper = "Hermitian positive definite linear system solve (Cholesky)";
    test_zposv.params = zposv_params;
    test_zposv.output = zposv_output;
    test_zposv.outchk = zposv_outchk;
    test_zposv.fptr_desc = testing_zposv_desc;
    test_zposv.fptr_std  = testing_zposv_std;
    test_zposv.next   = NULL;

    testing_register( &test_zposv );
}
