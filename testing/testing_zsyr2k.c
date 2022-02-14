/**
 *
 * @file testing_zsyr2k.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyr2k testing
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
testing_zsyr2k_desc( run_arg_list_t *args, int check )
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
    int          LDB    = run_arg_get_int( args, "LDB", ( ( trans == ChamNoTrans ) ? N : K ) );
    int          LDC    = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    double                bump  = testing_dalea();
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;

    bump  = run_arg_get_double( args, "bump", bump );
    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, An, 0, 0, Am, An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplgsy_Tile( bump, uplo, descC, seedC );

    /* Calculate the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zsyr2k_Tile_Async( uplo, trans, alpha, descA, descB, beta, descC,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zsyr2k_Tile( uplo, trans, alpha, descA, descB, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zher2k( K, N ) );

    /* Check the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplgsy_Tile( bump, uplo, descCinit, seedC );

        hres += check_zsyrk( args, ChamSymmetric, uplo, trans, alpha, descA, descB,
                             beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

int
testing_zsyr2k_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int                   nb    = run_arg_get_int( args, "nb", 320 );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   K     = run_arg_get_int( args, "K", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDB   = run_arg_get_int( args, "LDB", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDC   = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    double                bump  = testing_dalea();

    /* Descriptors */
    int                    Am, An;
    CHAMELEON_Complex64_t *A, *B, *C, *Cinit;

    bump  = run_arg_get_double( args, "bump", bump );
    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Create the matrices */
    A = malloc( LDA*An*sizeof(CHAMELEON_Complex64_t) );
    B = malloc( LDB*An*sizeof(CHAMELEON_Complex64_t) );
    C = malloc( LDC*N *sizeof(CHAMELEON_Complex64_t) );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
    CHAMELEON_zplrnt( Am, An, B, LDB, seedB );
    CHAMELEON_zplgsy( bump, uplo, N, C, LDC, seedC );

    /* Calculate the product */
    testing_start( &test_data );
    hres = CHAMELEON_zsyr2k( uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zher2k( K, N ) );

    /* Check the solution */
    if ( check ) {
        Cinit = malloc( LDC*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplgsy( bump, uplo, N, Cinit, LDC, seedC );

        hres += check_zsyrk_std( args, ChamSymmetric, uplo, trans, N, K, alpha, A, B, LDA, beta, Cinit, C, LDC );

        free( Cinit );
    }

    free( A );
    free( B );
    free( C );

    return hres;
}

testing_t   test_zsyr2k;
const char *zsyr2k_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "n",    "k",
                                "lda",    "ldb",   "ldc",   "alpha", "beta", "seedA",
                                "seedB",  "seedC", "bump",  NULL };
const char *zsyr2k_output[] = { NULL };
const char *zsyr2k_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsyr2k_init( void ) __attribute__( ( constructor ) );
void
testing_zsyr2k_init( void )
{
    test_zsyr2k.name   = "zsyr2k";
    test_zsyr2k.helper = "Symmetrix matrix-matrix rank 2k update";
    test_zsyr2k.params = zsyr2k_params;
    test_zsyr2k.output = zsyr2k_output;
    test_zsyr2k.outchk = zsyr2k_outchk;
    test_zsyr2k.fptr_desc = testing_zsyr2k_desc;
    test_zsyr2k.fptr_std  = testing_zsyr2k_std;
    test_zsyr2k.next   = NULL;

    testing_register( &test_zsyr2k );
}
