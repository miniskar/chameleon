/**
 *
 * @file testing_zgram.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram testing
 *
 * @version 1.2.0
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zgram( int _N )
{
    cham_fixdbl_t flops = 0.;
    cham_fixdbl_t N     = _N;

#if defined(PRECISION_z) || defined(PRECISION_c)
    /* 5 multiplications and 3 addition per element */
    flops = ( 5. * 6. + 6. ) * N * N;
#else
    flops = ( 5. + 3. ) * N * N;
#endif

    return flops;
}

int
testing_zgram_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         async  = parameters_getvalue_int( "async" );
    intptr_t    mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int         nb     = run_arg_get_int( args, "nb", 320 );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         P      = parameters_getvalue_int( "P" );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         LDA    = run_arg_get_int( args, "LDA", N );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA;
    void        *ws = NULL;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

    /* Fill the matrix with random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    if ( async ) {
        ws = CHAMELEON_zgram_WS_Alloc( descA );
    }

    /* Compute the gram matrix transformation */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgram_Tile_Async( uplo, descA, ws,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgram_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgram( N ) );

    if ( ws != NULL ) {
        CHAMELEON_zgemm_WS_Free( ws );
    }

    CHAMELEON_Desc_Destroy( &descA );

    (void)check;
    return hres;
}

int
testing_zgram_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         nb    = run_arg_get_int( args, "nb", 320 );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", random() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

    /* Fill the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N,A, LDA, seedA );

    /* Compute the gram matrix transformation */
    testing_start( &test_data );
    hres = CHAMELEON_zgram( uplo, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgram( N ) );

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zgram;
const char *zgram_params[] = { "mtxfmt", "nb", "uplo", "n", "n", "lda", "seedA", NULL };
const char *zgram_output[] = { NULL };
const char *zgram_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgram_init( void ) __attribute__( ( constructor ) );
void
testing_zgram_init( void )
{
    test_zgram.name   = "zgram";
    test_zgram.helper = "General Gram matrix transformation";
    test_zgram.params = zgram_params;
    test_zgram.output = zgram_output;
    test_zgram.outchk = zgram_outchk;
    test_zgram.fptr_desc = testing_zgram_desc;
    test_zgram.fptr_std  = testing_zgram_std;
    test_zgram.next   = NULL;

    testing_register( &test_zgram );
}
