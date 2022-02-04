/**
 *
 * @file testing_zsytrs.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrs testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zsytrs_desc( run_arg_list_t *args, int check )
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
    CHAMELEON_zplgsy_Tile( (double)N, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    hres = CHAMELEON_zsytrf_Tile( uplo, descA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zsytrs_Tile_Async( uplo, descA, descX,
                                             test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zsytrs_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, 0 /*flops_zsytrs( N, NRHS )*/ );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, NULL );

        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres += check_zsolve( args, ChamSymmetric, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    return hres;
}

testing_t   test_zsytrs;
const char *zsytrs_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                                "lda",    "ldb", "seedA", "seedB", NULL };
const char *zsytrs_output[] = { NULL };
const char *zsytrs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsytrs_init( void ) __attribute__( ( constructor ) );
void
testing_zsytrs_init( void )
{
    test_zsytrs.name   = "zsytrs";
    test_zsytrs.helper = "Symmetric triangular solve";
    test_zsytrs.params = zsytrs_params;
    test_zsytrs.output = zsytrs_output;
    test_zsytrs.outchk = zsytrs_outchk;
    test_zsytrs.fptr_desc = testing_zsytrs_desc;
    test_zsytrs.fptr_std  = NULL;
    test_zsytrs.next   = NULL;

    testing_register( &test_zsytrs );
}
