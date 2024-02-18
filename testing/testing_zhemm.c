/**
 *
 * @file testing_zhemm.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm testing
 *
 * @version 1.3.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2023-07-05
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if defined(CHAMELEON_TESTINGS_VENDOR) || !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/cblas.h>
#endif

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zhemm_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    int                   P     = parameters_getvalue_int( "P" );
    cham_side_t           side  = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( side == ChamLeft ) ? M : N ) );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    int                   LDC   = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );
    double                bump  = testing_dalea();
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;
    void        *ws = NULL;

    bump  = run_arg_get_double( args, "bump", bump );
    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the side */
    if ( side == ChamLeft ) {
        Am = M;
    }
    else {
        Am = N;
    }

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, Am, Am, Am );
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, N, M, N );
    parameters_desc_create( "C", &descC, ChamComplexDouble, nb, nb, LDC, N, M, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( bump, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    if ( async ) {
        ws = CHAMELEON_zhemm_WS_Alloc( side, uplo, descA, descB, descC );
    }

    /* Calculates the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zhemm_Tile_Async( side, uplo, alpha, descA, descB, beta, descC,
                                           ws, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zhemm_Tile( side, uplo, alpha, descA, descB, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zhemm( side, M, N ) );

    if ( ws != NULL ) {
        CHAMELEON_zhemm_WS_Free( ws );
    }

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descCinit, seedC );

        hres +=
            check_zsymm( args, ChamHermitian, side, uplo, alpha, descA, descB, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );
    parameters_desc_destroy( &descC );

    return hres;
}
#endif

int
testing_zhemm_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int                   api   = parameters_getvalue_int( "api" );
    int                   nb    = run_arg_get_nb( args );
    cham_side_t           side  = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( side == ChamLeft ) ? M : N ) );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    int                   LDC   = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );
    double                bump  = testing_dalea();

    /* Descriptors */
    int                    An;
    CHAMELEON_Complex64_t *A, *B, *C;

    bump  = run_arg_get_double( args, "bump", bump );
    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the side */
    An = ( side == ChamLeft ) ? M : N;

    /* Create the matrices */
    A = malloc( LDA*An*sizeof(CHAMELEON_Complex64_t) );
    B = malloc( LDB*N *sizeof(CHAMELEON_Complex64_t) );
    C = malloc( LDC*N *sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe( bump, uplo, An, A, LDA, seedA );
    CHAMELEON_zplrnt( M, N, B, LDB, seedB );
    CHAMELEON_zplrnt( M, N, C, LDC, seedC );

    /* Calculates the product */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    cblas_zhemm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, M, N,
                        CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), C, LDC );
    testing_stop( &test_data, flops_zhemm( side, M, N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zhemm( side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC );
        break;
    case 2:
        CHAMELEON_cblas_zhemm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, M, N,
                               CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), C, LDC );
        break;
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zhemm( side, M, N ) );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Complex64_t *Cinit;
        Cinit = malloc( LDC*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplrnt( M, N, Cinit, LDC, seedC );

        hres += check_zsymm_std( args, ChamHermitian, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, Cinit, C, LDC );

        free( Cinit );
    }
#endif

    free( A );
    free( B );
    free( C );

    (void)check;
    return hres;
}

testing_t   test_zhemm;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zhemm_params[] = { "side", "uplo",  "m",     "n",     "lda",  "ldb",
                               "ldc",    "alpha", "beta", "seedA", "seedB", "seedC", "bump", NULL };
#else
const char *zhemm_params[] = { "mtxfmt", "nb",    "side", "uplo",  "m",     "n",     "lda",  "ldb",
                               "ldc",    "alpha", "beta", "seedA", "seedB", "seedC", "bump", NULL };
#endif
const char *zhemm_output[] = { NULL };
const char *zhemm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zhemm_init( void ) __attribute__( ( constructor ) );
void
testing_zhemm_init( void )
{
    test_zhemm.name   = "zhemm";
    test_zhemm.helper = "Hermitian matrix-matrix multiply";
    test_zhemm.params = zhemm_params;
    test_zhemm.output = zhemm_output;
    test_zhemm.outchk = zhemm_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zhemm.fptr_desc = NULL;
#else
    test_zhemm.fptr_desc = testing_zhemm_desc;
#endif
    test_zhemm.fptr_std  = testing_zhemm_std;
    test_zhemm.next   = NULL;

    testing_register( &test_zhemm );
}
