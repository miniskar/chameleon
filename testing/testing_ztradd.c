/**
 *
 * @file testing_ztradd.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd testing
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
flops_ztradd( cham_uplo_t uplo, int M, int N )
{
    cham_fixdbl_t flops = 0.;
    int           minMN = chameleon_min( M, N );
    switch ( uplo ) {
        case ChamUpper:
            flops = ( minMN * ( minMN + 1 ) / 2 ) + M * chameleon_max( 0, N - M );
            break;
        case ChamLower:
            flops = ( minMN * ( minMN + 1 ) / 2 ) + N * chameleon_max( 0, M - N );
            break;
        case ChamUpperLower:
        default:
            flops = M * N;
    }

#if defined( PRECISION_z ) || defined( PRECISION_c )
    /* 2 multiplications and 1 addition per element */
    flops *= ( 2. * 6. + 2. );
#else
    flops *= ( 2. + 1. );
#endif

    return flops;
}

int
testing_ztradd_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_int( args, "nb", 320 );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t  uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? M : N ) );
    int          LDB    = run_arg_get_int( args, "LDB", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   Q     = parameters_compute_q( P );

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

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    switch ( uplo ) {
        case ChamUpper:
        case ChamLower:
            if ( trans == ChamNoTrans ) {
                CHAMELEON_zplgsy_Tile( 0., uplo, descA, seedA );
            }
            else {
                cham_uplo_t uplo_inv = ( uplo == ChamUpper ) ? ChamLower : ChamUpper;
                CHAMELEON_zplgsy_Tile( 0., uplo_inv, descA, seedA );
            }
            CHAMELEON_zplgsy_Tile( 0., uplo, descB, seedB );
            break;
        case ChamUpperLower:
        default:
            CHAMELEON_zplrnt_Tile( descA, seedA );
            CHAMELEON_zplrnt_Tile( descB, seedB );
            break;
    }

    /* Calculates the sum */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_ztradd_Tile_Async( uplo, trans, alpha, descA, beta, descB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_ztradd_Tile( uplo, trans, alpha, descA, beta, descB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztradd( uplo, M, N ) );

    /* Checks the solution */
    if ( check ) {
        CHAM_desc_t *descB0 = CHAMELEON_Desc_Copy( descB, NULL );

        if ( uplo == ChamUpperLower ) {
            CHAMELEON_zplrnt_Tile( descB0, seedB );
        }
        else {
            CHAMELEON_zplgsy_Tile( 0., uplo, descB0, seedB );
        }
        hres += check_zsum( args, uplo, trans, alpha, descA, beta, descB0, descB );

        CHAMELEON_Desc_Destroy( &descB0 );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );

    return hres;
}

testing_t   test_ztradd;
const char *ztradd_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "m",     "n", "lda",
                                "ldb",    "alpha", "beta",  "seedA", "seedB", NULL };
const char *ztradd_output[] = { NULL };
const char *ztradd_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztradd_init( void ) __attribute__( ( constructor ) );
void
testing_ztradd_init( void )
{
    test_ztradd.name   = "ztradd";
    test_ztradd.helper = "Triangular matrix-matrix addition";
    test_ztradd.params = ztradd_params;
    test_ztradd.output = ztradd_output;
    test_ztradd.outchk = ztradd_outchk;
    test_ztradd.fptr_desc = testing_ztradd_desc;
    test_ztradd.fptr_std  = NULL;
    test_ztradd.next   = NULL;

    testing_register( &test_ztradd );
}
