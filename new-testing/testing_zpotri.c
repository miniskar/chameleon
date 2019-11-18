/**
 *
 * @file testing_zpotri.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotri testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-08-13
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <assert.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"

int
testing_zpotri( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres;
    CHAM_desc_t *descA;

    /* Reads arguments */
    int         nb    = run_arg_get_int( args, "nb", 320 );
    int         P     = parameters_getvalue_int( "P" );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", random() );
    int         Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zpotri( N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

    /* Initialise the matrix with the random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the inversed matrix */
    START_TIMING( t );
    hres = CHAMELEON_zpotrf_Tile( uplo, descA );
    assert( hres == 0 );
    hres = CHAMELEON_zpotri_Tile( uplo, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", gflops );

    /* Check the inverse */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_ztrtri( args, ChamHermitian, uplo, ChamNonUnit, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    run_id++;
    return hres;
}

testing_t   test_zpotri;
const char *zpotri_params[] = { "nb", "uplo", "n", "lda", "seedA", NULL };
const char *zpotri_output[] = { NULL };
const char *zpotri_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotri_init( void ) __attribute__( ( constructor ) );
void
testing_zpotri_init( void )
{
    test_zpotri.name        = "zpotri";
    test_zpotri.helper      = "zpotri";
    test_zpotri.params      = zpotri_params;
    test_zpotri.output      = zpotri_output;
    test_zpotri.outchk      = zpotri_outchk;
    test_zpotri.params_list = "nb;P;uplo;n;lda;seedA";
    test_zpotri.fptr        = testing_zpotri;
    test_zpotri.next        = NULL;

    testing_register( &test_zpotri );
}
