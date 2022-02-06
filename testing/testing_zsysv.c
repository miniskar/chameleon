/**
 *
 * @file testing_zsysv.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsysv testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zsysv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotrs( N, NRHS );
    return flops;
}

int
testing_zsysv_desc( run_arg_list_t *args, int check )
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

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zsysv_Tile_Async( uplo, descA, descX,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zsysv_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zsysv( N, NRHS ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );

        hres += check_zxxtrf( args, ChamSymmetric, uplo, descA0, descA );

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, NULL );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );
        hres += check_zsolve( args, ChamSymmetric, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    return hres;
}

testing_t   test_zsysv;
const char *zsysv_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                               "lda",    "ldb", "seedA", "seedB", NULL };
const char *zsysv_output[] = { NULL };
const char *zsysv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsysv_init( void ) __attribute__( ( constructor ) );
void
testing_zsysv_init( void )
{
    test_zsysv.name   = "zsysv";
    test_zsysv.helper = "Symmetrix linear system solve";
    test_zsysv.params = zsysv_params;
    test_zsysv.output = zsysv_output;
    test_zsysv.outchk = zsysv_outchk;
    test_zsysv.fptr_desc = testing_zsysv_desc;
    test_zsysv.fptr_std  = NULL;
    test_zsysv.next   = NULL;

    testing_register( &test_zsysv );
}
