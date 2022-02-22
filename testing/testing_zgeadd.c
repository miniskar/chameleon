/**
 *
 * @file testing_zgeadd.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd testing
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
flops_zgeadd( int M, int N )
{
    cham_fixdbl_t flops = 0.;
#if defined( PRECISION_z ) || defined( PRECISION_c )
    /* 2 multiplications and 1 addition per element */
    flops = ( 2. * 6. + 2. ) * M * N;
#else
    flops = ( 2. + 1. ) * M * N;
#endif

    return flops;
}

int
testing_zgeadd_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_int( args, "nb", 320 );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? M : N ) );
    int          LDB    = run_arg_get_int( args, "LDB", M );
    int          seedA  = run_arg_get_int( args, "seedA", random() );
    int          seedB  = run_arg_get_int( args, "seedB", random() );
    int          Q      = parameters_compute_q( P );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descB;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    if ( trans != ChamNoTrans ) {
        Am = N;
        An = M;
    }
    else {
        Am = M;
        An = N;
    }

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Compute the sum */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgeadd_Tile_Async( trans, alpha, descA, beta, descB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgeadd_Tile( trans, alpha, descA, beta, descB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgeadd( M, N ) );

    /* Check the solution */
    if ( check ) {
        CHAM_desc_t *descB0 = CHAMELEON_Desc_Copy( descB, NULL );
        CHAMELEON_zplrnt_Tile( descB0, seedB );

        hres += check_zsum( args, ChamUpperLower, trans, alpha, descA, beta, descB0, descB );

        CHAMELEON_Desc_Destroy( &descB0 );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );

    return hres;
}

testing_t   test_zgeadd;
const char *zgeadd_params[] = { "mtxfmt", "nb",    "trans", "m",     "n",     "lda",
                                "ldb",    "alpha", "beta",  "seedA", "seedB", NULL };
const char *zgeadd_output[] = { NULL };
const char *zgeadd_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgeadd_init( void ) __attribute__( ( constructor ) );
void
testing_zgeadd_init( void )
{
    test_zgeadd.name   = "zgeadd";
    test_zgeadd.helper = "General matrix-matrix addition";
    test_zgeadd.params = zgeadd_params;
    test_zgeadd.output = zgeadd_output;
    test_zgeadd.outchk = zgeadd_outchk;
    test_zgeadd.fptr_desc = testing_zgeadd_desc;
    test_zgeadd.fptr_std  = NULL;
    test_zgeadd.next   = NULL;

    testing_register( &test_zgeadd );
}
