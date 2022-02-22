/**
 *
 * @file testing_zgesv.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesv testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zgesv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zgetrf( N, N ) + flops_zgetrs( N, NRHS );
    return flops;
}

int
testing_zgesv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_int( args, "nb", 320 );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int      LDA    = run_arg_get_int( args, "LDA", N );
    int      LDB    = run_arg_get_int( args, "LDB", N );
    int      seedA  = run_arg_get_int( args, "seedA", random() );
    int      seedB  = run_arg_get_int( args, "seedB", random() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descX;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgesv_nopiv_Tile_Async( descA, descX,
                                                 test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgesv_nopiv_Tile( descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgesv( N, NRHS ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        hres += check_zxxtrf( args, ChamGeneral, ChamUpperLower, descA0, descA );

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, NULL );
        CHAMELEON_zplrnt_Tile( descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres +=
            check_zsolve( args, ChamGeneral, ChamNoTrans, ChamUpperLower, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    return hres;
}

testing_t   test_zgesv;
const char *zgesv_params[] = { "mtxfmt", "nb", "n", "nrhs", "lda", "ldb", "seedA", "seedB", NULL };
const char *zgesv_output[] = { NULL };
const char *zgesv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgesv_init( void ) __attribute__( ( constructor ) );
void
testing_zgesv_init( void )
{
    test_zgesv.name   = "zgesv";
    test_zgesv.helper = "General linear system solve (LU without pivoting)";
    test_zgesv.params = zgesv_params;
    test_zgesv.output = zgesv_output;
    test_zgesv.outchk = zgesv_outchk;
    test_zgesv.fptr_desc = testing_zgesv_desc;
    test_zgesv.fptr_std  = NULL;
    test_zgesv.next   = NULL;

    testing_register( &test_zgesv );
}
