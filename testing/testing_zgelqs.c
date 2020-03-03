/**
 *
 * @file testing_zgelqs.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqs testing
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include "flops.h"
#include "control/common.h"

static cham_fixdbl_t
flops_zgelqs()
{
    cham_fixdbl_t flops = 0.;
    return flops;
}

int
testing_zgelqs( run_arg_list_t *args, int check )
{
    static int   run_id = 0;
    int          hres   = 0;
    CHAM_desc_t *descA1, *descA2, *descB1, *descB2, *descT, *descQ, *descX;

    /* Reads arguments */
    int    nb    = run_arg_get_int( args, "nb", 320 );
    int    ib    = run_arg_get_int( args, "ib", 48 );
    int    P     = parameters_getvalue_int( "P" );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    LDB   = run_arg_get_int( args, "LDB", M );
    int    RH    = run_arg_get_int( args, "qra", 0 );
    int    seedA = run_arg_get_int( args, "seedA", random() );
    int    seedB = run_arg_get_int( args, "seedB", random() );
    int    Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgelqs();

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( M >= N ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: The LQ solution is performed only when N > M\n" );
        }
        return -1;
    }

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
        &descX, NULL, ChamComplexDouble, nb, nb, nb * nb, LDB, K, 0, 0, M, K, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    hres = CHAMELEON_zgelqf_Tile( descA, descT );

    /* Checks the factorisation, orthogonality and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, NULL );
        CHAM_desc_t *subX   = chameleon_desc_submatrix( descX, 0, 0, N, NRHS );
        CHAM_desc_t *subB   = chameleon_desc_submatrix( descB, 0, 0, M, NRHS );

        CHAMELEON_zplrnt_Tile( descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres += check_zsolve( args, ChamGeneral, ChamNoTrans, ChamUpperLower, descA0, subX, subB );

        free( subB );
        free( subX );
        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );
    CHAMELEON_Desc_Destroy( &descT );

    run_id++;
    return hres;
}

testing_t   test_zgelqs;
const char *zgelqs_params[] = { "nb",  "ib", "m",     "n",     "k", "lda",
                                "ldb", "qra", "seedA", "seedB", NULL };
const char *zgelqs_output[] = { NULL };
const char *zgelqs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgelqs_init( void ) __attribute__( ( constructor ) );
void
testing_zgelqs_init( void )
{
    test_zgelqs.name        = "zgelqs";
    test_zgelqs.helper      = "General LQ solve";
    test_zgelqs.params      = zgelqs_params;
    test_zgelqs.output      = zgelqs_output;
    test_zgelqs.outchk      = zgelqs_outchk;
    test_zgelqs.params_list = "nb;ib;P;m;n;k;lda;ldb;rh;seedA;seedB";
    test_zgelqs.fptr        = testing_zgelqs;
    test_zgelqs.next        = NULL;

    testing_register( &test_zgelqs );
}
