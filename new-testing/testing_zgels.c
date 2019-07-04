/**
 *
 * @file testing_zgels.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgels testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-08-28
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testing_zauxiliary.h"
#include "testing_zcheck.h"
#include "flops.h"
#include "../control/common.h"

static cham_fixdbl_t
flops_zgels( cham_trans_t trans, int M, int N, int NRHS )
{
    cham_fixdbl_t flops = 0.;
    return flops;
}

int
testing_zgels( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA, *descX, *descT;

    /* Reads arguments */
    int          nb    = run_arg_get_int( args, "nb", 320 );
    int          ib    = run_arg_get_int( args, "ib", 48 );
    int          P     = parameters_getvalue_int( "P" );
    cham_trans_t trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N     = run_arg_get_int( args, "N", 1000 );
    int          M     = run_arg_get_int( args, "M", N );
    int          maxMN = chameleon_max( M, N );
    int          NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int          LDA   = run_arg_get_int( args, "LDA", M );
    int          LDB   = run_arg_get_int( args, "LDB", maxMN );
    int          RH    = run_arg_get_int( args, "qra", 4 );
    int          seedA = run_arg_get_int( args, "seedA", random() );
    int          seedB = run_arg_get_int( args, "seedB", random() );
    int          Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgels( trans, M, N, NRHS );

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
        &descA, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, maxMN, NRHS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Computes the solution */
    START_TIMING( t );
    hres = CHAMELEON_zgels_Tile( trans, descA, descT, descX );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", gflops );

    if ( check ) {
        CHAM_desc_t *descA0, *descB;
        CHAM_desc_t *subX, *subB;

        CHAMELEON_Desc_Create(
            &descA0, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
        CHAMELEON_Desc_Create(
            &descB, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, maxMN, NRHS, P, Q );

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
        hres = check_zgels( trans, descA0, subX, subB );  // A(Ax-B) == 0?

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );

        free( subB );
        free( subX );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descX );

    run_id++;
    return hres;
}

testing_t   test_zgels;
const char *zgels_params[] = { "nb",  "ib",  "trans", "m",     "n",     "k",
                               "lda", "ldb", "qra",    "seedA", "seedB", NULL };
const char *zgels_output[] = { NULL };
const char *zgels_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgels_init( void ) __attribute__( ( constructor ) );
void
testing_zgels_init( void )
{
    test_zgels.name        = "zgels";
    test_zgels.helper      = "zgels";
    test_zgels.params      = zgels_params;
    test_zgels.output      = zgels_output;
    test_zgels.outchk      = zgels_outchk;
    test_zgels.params_list = "nb;ib;P;trans;m;n;k;lda;ldb;rh;seedA;seedB";
    test_zgels.fptr        = testing_zgels;
    test_zgels.next        = NULL;

    testing_register( &test_zgels );
}
