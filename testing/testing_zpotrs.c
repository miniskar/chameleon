/**
 *
 * @file testing_zpotrs.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrs testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-07
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zpotrs_desc( run_arg_list_t *args, int check )
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

    hres = CHAMELEON_zpotrf_Tile( uplo, descA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zpotrs_Tile_Async( uplo, descA, descX,
                                             test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zpotrs_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrs( N, NRHS ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, NULL );

        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres += check_zsolve( args, ChamHermitian, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    return hres;
}

int
testing_zpotrs_std( run_arg_list_t *args, int check )
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
    A = malloc( LDA*N*   sizeof(CHAMELEON_Complex64_t) );
    X = malloc( LDB*NRHS*sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

    hres = CHAMELEON_zpotrf( uplo, N, A, LDA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    hres += CHAMELEON_zpotrs( uplo, N, NRHS, A, LDA, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrs( N, NRHS ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAMELEON_Complex64_t *A0 = malloc( LDA*N*   sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_Complex64_t *B  = malloc( LDB*NRHS*sizeof(CHAMELEON_Complex64_t) );

        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );
        CHAMELEON_zplrnt( N, NRHS, B, LDB, seedB );

        hres += check_zsolve_std( args, ChamHermitian, ChamNoTrans, uplo, N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }

    free( A );
    free( X );

    return hres;
}

testing_t   test_zpotrs;
const char *zpotrs_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                                "lda",    "ldb", "seedA", "seedB", NULL };
const char *zpotrs_output[] = { NULL };
const char *zpotrs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotrs_init( void ) __attribute__( ( constructor ) );
void
testing_zpotrs_init( void )
{
    test_zpotrs.name   = "zpotrs";
    test_zpotrs.helper = "Hermitian positive definite solve (Cholesky)";
    test_zpotrs.params = zpotrs_params;
    test_zpotrs.output = zpotrs_output;
    test_zpotrs.outchk = zpotrs_outchk;
    test_zpotrs.fptr_desc = testing_zpotrs_desc;
    test_zpotrs.fptr_std  = testing_zpotrs_std;
    test_zpotrs.next   = NULL;

    testing_register( &test_zpotrs );
}
