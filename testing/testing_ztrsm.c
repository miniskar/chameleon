/**
 *
 * @file testing_ztrsm.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm testing
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_ztrsm_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int                   async  = parameters_getvalue_int( "async" );
    intptr_t              mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int                   nb     = run_arg_get_int( args, "nb", 320 );
    int                   P      = parameters_getvalue_int( "P" );
    cham_trans_t          trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t           side   = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t           diag   = run_arg_get_diag( args, "diag", ChamNonUnit );
    int                   N      = run_arg_get_int( args, "N", 1000 );
    int                   M      = run_arg_get_int( args, "M", N );
    int                   Ak     = ( side == ChamLeft ) ? M : N;
    int                   LDA    = run_arg_get_int( args, "LDA", Ak );
    int                   LDB    = run_arg_get_int( args, "LDB", N );
    CHAMELEON_Complex64_t alpha  = testing_zalea();
    int                   seedA  = run_arg_get_int( args, "seedA", random() );
    int                   seedB  = run_arg_get_int( args, "seedB", random() );
    int                   Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descB, *descBinit;

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create(
        &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    /* We bump a little bit the diagonal to make it stable */
    CHAMELEON_zplgsy_Tile( 2., uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Calculates the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_ztrsm_Tile_Async( side, uplo, trans, diag, alpha, descA, descB,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_ztrsm_Tile( side, uplo, trans, diag, alpha, descA, descB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztrsm( side, M, N ) );

    /* Checks the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descBinit, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descBinit, seedB );

        hres += check_ztrmm( args, CHECK_TRSM, side, uplo, trans, diag,
                             alpha, descA, descB, descBinit );

        CHAMELEON_Desc_Destroy( &descBinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );

    return hres;
}

int
testing_ztrsm_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int                   nb    = run_arg_get_int( args, "nb", 320 );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t           side  = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t           diag  = run_arg_get_diag( args, "diag", ChamNonUnit );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   An    = ( side == ChamLeft ) ? M : N;
    int                   LDA   = run_arg_get_int( args, "LDA", An );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *B, *Binit;

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( LDA*An*sizeof(CHAMELEON_Complex64_t) );
    B = malloc( LDB*N *sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    /* We bump a little bit the diagonal to make it stable */
    CHAMELEON_zplgsy( 2., uplo, N, A, LDA, seedA );
    CHAMELEON_zplrnt( M, N, B, LDB, seedB );

    /* Calculates the product */
    testing_start( &test_data );
    hres = CHAMELEON_ztrsm( side, uplo, trans, diag, M, N, alpha, A, LDA, B, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztrsm( side, M, N ) );

    /* Checks the solution */
    if ( check ) {
        Binit = malloc( LDB*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplrnt( M, N, Binit, LDB, seedB );

        hres += check_ztrmm_std( args, CHECK_TRSM, side, uplo, trans, diag, M, N, alpha, A, LDA, B, Binit, LDB );

        free( Binit );
    }

    free( A );
    free( B );

    return hres;
}

testing_t   test_ztrsm;
const char *ztrsm_params[] = { "mtxfmt", "nb",  "side", "uplo",  "trans", "diag",  "m",
                               "n",      "lda", "ldb",  "alpha", "seedA", "seedB", NULL };
const char *ztrsm_output[] = { NULL };
const char *ztrsm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrsm_init( void ) __attribute__( ( constructor ) );
void
testing_ztrsm_init( void )
{
    test_ztrsm.name   = "ztrsm";
    test_ztrsm.helper = "Triangular matrix solve";
    test_ztrsm.params = ztrsm_params;
    test_ztrsm.output = ztrsm_output;
    test_ztrsm.outchk = ztrsm_outchk;
    test_ztrsm.fptr_desc = testing_ztrsm_desc;
    test_ztrsm.fptr_std  = testing_ztrsm_std;
    test_ztrsm.next   = NULL;

    testing_register( &test_ztrsm );
}
