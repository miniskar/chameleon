/**
 *
 * @file testing_zgeqrs.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrs testing
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
#include <chameleon/flops.h>
#include "control/common.h"

int
testing_zgeqrs( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA, *descX, *descT;

    /* Reads arguments */
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_int( args, "nb", 320 );
    int      ib     = run_arg_get_int( args, "ib", 48 );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      LDB    = run_arg_get_int( args, "LDB", M );
    int      RH     = run_arg_get_int( args, "qra", 0 );
    int      seedA  = run_arg_get_int( args, "seedA", random() );
    int      seedB  = run_arg_get_int( args, "seedB", random() );
    int      Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgeqrs( M, N, NRHS );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( N > M ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: The QR solution is performed only when M >= N\n" );
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
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, M, NRHS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    hres = CHAMELEON_zgeqrf_Tile( descA, descT );

    START_TIMING( t );
    hres = CHAMELEON_zgeqrs_Tile( descA, descT, descX );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

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

    return hres;
}

testing_t   test_zgeqrs;
const char *zgeqrs_params[] = { "mtxfmt", "nb", "ib", "m",     "n",     "k", "lda",
                                "ldb", "qra", "seedA", "seedB", NULL };
const char *zgeqrs_output[] = { NULL };
const char *zgeqrs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgeqrs_init( void ) __attribute__( ( constructor ) );
void
testing_zgeqrs_init( void )
{
    test_zgeqrs.name        = "zgeqrs";
    test_zgeqrs.helper      = "General QR solve";
    test_zgeqrs.params      = zgeqrs_params;
    test_zgeqrs.output      = zgeqrs_output;
    test_zgeqrs.outchk      = zgeqrs_outchk;
    test_zgeqrs.fptr        = testing_zgeqrs;
    test_zgeqrs.next        = NULL;

    testing_register( &test_zgeqrs );
}
