/**
 *
 * @file testing_zherk.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-04
 * @precisions normal z -> z c
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zherk_desc( run_arg_list_t *args, int check )
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
    double       alpha  = testing_dalea();
    double       beta   = testing_dalea();
    double       bump   = testing_dalea();
    int          seedA  = run_arg_get_int( args, "seedA", random() );
    int          seedC  = run_arg_get_int( args, "seedC", random() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descC, *descCinit;

    alpha = run_arg_get_double( args, "alpha", alpha );
    beta  = run_arg_get_double( args, "beta", beta );
    bump  = run_arg_get_double( args, "bump", bump );

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
    CHAMELEON_zplghe_Tile( bump, uplo, descC, seedC );

    /* Calculates the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zherk_Tile_Async( uplo, trans, alpha, descA, beta, descC,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zherk_Tile( uplo, trans, alpha, descA, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zherk( K, N ) );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplghe_Tile( bump, uplo, descCinit, seedC );

        hres += check_zsyrk( args, ChamHermitian, uplo, trans, alpha, descA, NULL,
                             beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

int
testing_zherk_std( run_arg_list_t *args, int check )
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
    double       alpha = testing_dalea();
    double       beta  = testing_dalea();
    double       bump  = testing_dalea();
    int          seedA = run_arg_get_int( args, "seedA", random() );
    int          seedC = run_arg_get_int( args, "seedC", random() );


    /* Descriptors */
    int                    Am, An;
    CHAMELEON_Complex64_t *A, *C, *Cinit;

    alpha = run_arg_get_double( args, "alpha", alpha );
    beta  = run_arg_get_double( args, "beta", beta );
    bump  = run_arg_get_double( args, "bump", bump );

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
    C = malloc( LDC*N*sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
    CHAMELEON_zplghe( bump, uplo, N, C, LDC, seedC );

    /* Calculates the product */
    testing_start( &test_data );
    hres = CHAMELEON_zherk( uplo, trans, N, K, alpha, A, LDA, beta, C, LDC );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zherk( K, N ) );

    /* Checks the solution */
    if ( check ) {
        Cinit = malloc( LDC*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplghe( bump, uplo, N, Cinit, LDC, seedC );

        // hres += check_zsyrk( args, ChamHermitian, uplo, trans, alpha, descA, NULL,
        //                      beta, descCinit, descC );

        free( Cinit );
    }

    free( A );
    free( C );

    return hres;
}

testing_t   test_zherk;
const char *zherk_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "n",     "k",    "lda",
                               "ldc",    "alpha", "beta",  "seedA", "seedC", "bump", NULL };
const char *zherk_output[] = { NULL };
const char *zherk_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zherk_init( void ) __attribute__( ( constructor ) );
void
testing_zherk_init( void )
{
    test_zherk.name   = "zherk";
    test_zherk.helper = "Hermitian matrix-matrix rank k update";
    test_zherk.params = zherk_params;
    test_zherk.output = zherk_output;
    test_zherk.outchk = zherk_outchk;
    test_zherk.fptr_desc = testing_zherk_desc;
    test_zherk.fptr_std  = NULL;
    test_zherk.next   = NULL;

    testing_register( &test_zherk );
}
