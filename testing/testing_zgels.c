/**
 *
 * @file testing_zgels.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgels testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include "../control/common.h"
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zgels( cham_trans_t trans, int M, int N, int NRHS )
{
    cham_fixdbl_t flops = 0.;
    (void)trans;
    (void)M;
    (void)N;
    (void)NRHS;
    return flops;
}

int
testing_zgels_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_int( args, "nb", 320 );
    int          ib     = run_arg_get_int( args, "ib", 48 );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          maxMN  = chameleon_max( M, N );
    int          NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int          LDA    = run_arg_get_int( args, "LDA", M );
    int          LDB    = run_arg_get_int( args, "LDB", maxMN );
    int          RH     = run_arg_get_int( args, "qra", 4 );
    int          seedA  = run_arg_get_int( args, "seedA", random() );
    int          seedB  = run_arg_get_int( args, "seedB", random() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descX, *descT;

    /* Make sure trans is only Notrans or ConjTrans */
    trans = ( trans == ChamNoTrans ) ? trans : ChamConjTrans;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, maxMN, NRHS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Computes the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgels_Tile_Async( trans, descA, descT, descX,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgels_Tile( trans, descA, descT, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgels( trans, M, N, NRHS ) );

    if ( check ) {
        CHAM_desc_t *descA0, *descB;
        CHAM_desc_t *subX, *subB;

        CHAMELEON_Desc_Create(
            &descA0, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
        CHAMELEON_Desc_Create(
            &descB, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, maxMN, NRHS, P, Q );

        CHAMELEON_zplrnt_Tile( descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        if ( trans == ChamNoTrans ) {
            subX = chameleon_desc_submatrix( descX, 0, 0, N, NRHS );
            subB = chameleon_desc_submatrix( descB, 0, 0, M, NRHS );
        }
        else {
            subX = chameleon_desc_submatrix( descX, 0, 0, M, NRHS );
            subB = chameleon_desc_submatrix( descB, 0, 0, N, NRHS );
        }

        /* Check the factorization and the residual */
        hres = check_zgels( args, trans, descA0, subX, subB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );

        free( subB );
        free( subX );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descX );

    return hres;
}

testing_t   test_zgels;
const char *zgels_params[] = { "mtxfmt", "nb",  "ib",  "trans", "m",     "n", "k",
                               "lda",    "ldb", "qra", "seedA", "seedB", NULL };
const char *zgels_output[] = { NULL };
const char *zgels_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgels_init( void ) __attribute__( ( constructor ) );
void
testing_zgels_init( void )
{
    test_zgels.name   = "zgels";
    test_zgels.helper = "Linear least squares with general matrix";
    test_zgels.params = zgels_params;
    test_zgels.output = zgels_output;
    test_zgels.outchk = zgels_outchk;
    test_zgels.fptr_desc = testing_zgels_desc;
    test_zgels.fptr_std  = NULL;
    test_zgels.next   = NULL;

    testing_register( &test_zgels );
}
