/**
 *
 * @file testing_zlascal.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2014-07-13
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

static cham_fixdbl_t
flops_zlascal( cham_uplo_t uplo, int M, int N )
{
    cham_fixdbl_t flops = 0.;
    int    minMN = chameleon_min( M, N );
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
    /* 1 multiplications per element */
    flops *= 6.;
#endif

    return flops;
}

int
testing_zlascal( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descAinit;

    /* Reads arguments */
    int                   nb    = run_arg_get_int( args, "nb", 320 );
    int                   P     = parameters_getvalue_int( "P" );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", M );
    CHAMELEON_Complex64_t alpha = run_arg_get_Complex64( args, "alpha", 1. );
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zlascal( uplo, M, N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Scales the matrix */
    START_TIMING( t );
    hres = CHAMELEON_zlascal_Tile( uplo, alpha, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", gflops );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descAinit, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descAinit, seedA );

        hres += check_zscale( uplo, alpha, descAinit, descA );

        CHAMELEON_Desc_Destroy( &descAinit );
    }

    CHAMELEON_Desc_Destroy( &descA );

    run_id++;
    return hres;
}

testing_t   test_zlascal;
const char *zlascal_params[] = { "nb", "uplo", "m", "n", "lda", "alpha", "seedA", NULL };
const char *zlascal_output[] = { NULL };
const char *zlascal_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlascal_init( void ) __attribute__( ( constructor ) );
void
testing_zlascal_init( void )
{
    test_zlascal.name        = "zlascal";
    test_zlascal.helper      = "zlascal";
    test_zlascal.params      = zlascal_params;
    test_zlascal.output      = zlascal_output;
    test_zlascal.outchk      = zlascal_outchk;
    test_zlascal.params_list = "nb;P;uplo;m;n;lda;alpha;seedA";
    test_zlascal.fptr        = testing_zlascal;
    test_zlascal.next        = NULL;

    testing_register( &test_zlascal );
}
