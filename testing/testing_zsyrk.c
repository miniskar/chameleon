/**
 *
 * @file testing_zsyrk.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-04
 * @precisions normal z -> z c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zsyrk_desc( run_arg_list_t *args, int check )
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
    int          K      = run_arg_get_int( args, "K", N );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int          LDC    = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    CHAMELEON_Complex64_t bump  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descC, *descCinit;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );
    bump  = run_arg_get_complex64( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplgsy_Tile( bump, uplo, descC, seedC );

    /* Calculates the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zsyrk_Tile_Async( uplo, trans, alpha, descA, beta, descC,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zsyrk_Tile( uplo, trans, alpha, descA, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zsyrk( K, N ) );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplgsy_Tile( bump, uplo, descCinit, seedC );

        hres += check_zsyrk( args, ChamSymmetric, uplo, trans, alpha, descA, NULL,
                             beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

int
testing_zsyrk_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          nb    = run_arg_get_int( args, "nb", 320 );
    cham_trans_t trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t  uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int          N     = run_arg_get_int( args, "N", 1000 );
    int          K     = run_arg_get_int( args, "K", N );
    int          LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int          LDC   = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    CHAMELEON_Complex64_t bump  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );

    /* Descriptors */
    int                    Am, An;
    CHAMELEON_Complex64_t *A, *C, *Cinit;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );
    bump  = run_arg_get_complex64( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Creates the matrices */
    A = malloc( LDA*An*sizeof(CHAMELEON_Complex64_t) );
    C = malloc( LDC*N *sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
    CHAMELEON_zplgsy( bump, uplo, N, C, LDC, seedC );

    /* Calculates the product */
    testing_start( &test_data );
    hres = CHAMELEON_zsyrk( uplo, trans, N, K, alpha, A, LDA, beta, C, LDC );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zsyrk( K, N ) );

    /* Checks the solution */
    if ( check ) {
        Cinit = malloc( LDC*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplgsy( bump, uplo, N, Cinit, LDC, seedC );

        hres += check_zsyrk_std( args, ChamSymmetric, uplo, trans, N, K, alpha, A, LDA, NULL, LDA, beta, Cinit, C, LDC );

        free( Cinit );
    }

    free( A );
    free( C );

    return hres;
}

testing_t   test_zsyrk;
const char *zsyrk_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "n",     "k",    "lda",
                               "ldc",    "alpha", "beta",  "seedA", "seedC", "bump", NULL };
const char *zsyrk_output[] = { NULL };
const char *zsyrk_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsyrk_init( void ) __attribute__( ( constructor ) );
void
testing_zsyrk_init( void )
{
    test_zsyrk.name   = "zsyrk";
    test_zsyrk.helper = "Symmetrix matrix-matrix rank k update";
    test_zsyrk.params = zsyrk_params;
    test_zsyrk.output = zsyrk_output;
    test_zsyrk.outchk = zsyrk_outchk;
    test_zsyrk.fptr_desc = testing_zsyrk_desc;
    test_zsyrk.fptr_std  = testing_zsyrk_std;
    test_zsyrk.next   = NULL;

    testing_register( &test_zsyrk );
}
