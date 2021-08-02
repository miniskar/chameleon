/**
 *
 * @file testing_zcesca.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca testing
 *
 * @version 1.1.0
 * @author Florent Pruvost
 * @date 2021-05-10
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zcesca( int M, int N )
{
    cham_fixdbl_t flops = 0.;
#if defined( PRECISION_z ) || defined( PRECISION_c )
    /*  2 multiplications and 5 addition per element */
    flops = ( 2. * 6. + 10. ) * M * N;
#else
    flops = ( 2. + 5. ) * M * N;
#endif

    return flops;
}

int
testing_zcesca( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA;

    /* Read arguments */
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_int( args, "nb", 320 );
    int          P      = parameters_getvalue_int( "P" );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          LDA    = run_arg_get_int( args, "LDA", M );
    int          seedA  = run_arg_get_int( args, "seedA", random() );
    int          Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zcesca( M, N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Compute the centered-scaled matrix transformation */
    START_TIMING( t );
    hres = CHAMELEON_zcesca_Tile( 1, 1, ChamColumnwise, descA, NULL, NULL );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zcesca;
const char *zcesca_params[] = { "mtxfmt", "nb", "trans", "m", "n", "lda", "seedA", NULL };
const char *zcesca_output[] = { NULL };
const char *zcesca_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zcesca_init( void ) __attribute__( ( constructor ) );
void
testing_zcesca_init( void )
{
    test_zcesca.name        = "zcesca";
    test_zcesca.helper      = "General centered-scaled matrix transformation";
    test_zcesca.params      = zcesca_params;
    test_zcesca.output      = zcesca_output;
    test_zcesca.outchk      = zcesca_outchk;
    test_zcesca.fptr        = testing_zcesca;
    test_zcesca.next        = NULL;

    testing_register( &test_zcesca );
}
