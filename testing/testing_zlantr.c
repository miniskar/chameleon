/**
 *
 * @file testing_zlantr.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlantr testing
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-14
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if defined(CHAMELEON_TESTINGS_VENDOR)
#include <coreblas/lapacke.h>
#include <coreblas.h>
#endif

static cham_fixdbl_t
flops_zlantr( cham_normtype_t ntype, cham_uplo_t uplo, int M, int N )
{
    cham_fixdbl_t flops   = 0.;
    cham_fixdbl_t coefabs = 1.;
    cham_fixdbl_t size;
#if defined(PRECISION_z) || defined(PRECISION_c)
    coefabs = 3.;
#endif

    switch ( uplo ) {
        case ChamUpper:
            if ( N > M ) {
                size = ( M * ( M + 1 ) / 2 ) + M * ( N - M );
            }
            else {
                size = N * ( N + 1 ) / 2;
            }
            break;
        case ChamLower:
            if ( M > N ) {
                size = ( N * ( N + 1 ) / 2 ) + N * ( M - N );
            }
            else {
                size = M * ( M + 1 ) / 2;
            }
            break;
        case ChamUpperLower:
        default:
            size = M * N;
    }
    flops = size * coefabs;

    switch ( ntype ) {
        case ChamOneNorm:
            flops += N;
            break;
        case ChamInfNorm:
            flops += M;
            break;
        case ChamMaxNorm:
            break;
        case ChamFrobeniusNorm:
            flops += size;
            break;
        default:;
    }
    (void)flops;
    return size;
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zlantr_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int             async     = parameters_getvalue_int( "async" );
    intptr_t        mtxfmt    = parameters_getvalue_int( "mtxfmt" );
    int             nb        = run_arg_get_int( args, "nb", 320 );
    int             P         = parameters_getvalue_int( "P" );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t     uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t     diag      = run_arg_get_diag( args, "diag", ChamNonUnit );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             M         = run_arg_get_int( args, "M", N );
    int             LDA       = run_arg_get_int( args, "LDA", M );
    int             seedA     = run_arg_get_int( args, "seedA", random() );
    int             Q         = parameters_compute_q( P );

    /* Descriptors */
    double       norm;
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the norm */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlantr_Tile_Async( norm_type, uplo, diag, descA, &norm,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        norm = CHAMELEON_zlantr_Tile( norm_type, uplo, diag, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlantr( norm_type, uplo, M, N ) );

    /* Checks the solution */
    if ( check ) {
        hres = check_znorm( args, ChamTriangular, norm_type, uplo, diag, norm, descA );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}
#endif

int
testing_zlantr_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int             nb        = run_arg_get_int( args, "nb", 320 );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t     uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t     diag      = run_arg_get_diag( args, "diag", ChamNonUnit );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             M         = run_arg_get_int( args, "M", N );
    int             LDA       = run_arg_get_int( args, "LDA", M );
    int             seedA     = run_arg_get_int( args, "seedA", random() );

    /* Descriptors */
    double norm;
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Calculates the norm */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    norm = LAPACKE_zlantr( LAPACK_COL_MAJOR, chameleon_lapack_const( norm_type ), uplo, diag, M, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlantr( norm_type, uplo, M, N ) );
#else
    testing_start( &test_data );
    norm = CHAMELEON_zlantr( norm_type, uplo, diag, M, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlantr( norm_type, uplo, M, N ) );

    /* Checks the solution */
    if ( check ) {
        hres = check_znorm_std( args, ChamTriangular, norm_type, uplo, diag, norm, M, N, A, LDA );
    }
#endif

    free( A );

    return hres;
}

testing_t   test_zlantr;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zlantr_params[] = { "norm", "uplo",  "diag",
                                "m",      "n",  "lda",  "seedA", NULL };
#else
const char *zlantr_params[] = { "mtxfmt", "nb", "norm", "uplo",  "diag",
                                "m",      "n",  "lda",  "seedA", NULL };
#endif
const char *zlantr_output[] = { NULL };
const char *zlantr_outchk[] = { "||A||", "||B||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlantr_init( void ) __attribute__( ( constructor ) );
void
testing_zlantr_init( void )
{
    test_zlantr.name   = "zlantr";
    test_zlantr.helper = "Triangular matrix norm";
    test_zlantr.params = zlantr_params;
    test_zlantr.output = zlantr_output;
    test_zlantr.outchk = zlantr_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zlantr.fptr_desc = NULL;
#else
    test_zlantr.fptr_desc = testing_zlantr_desc;
#endif
    test_zlantr.fptr_std  = testing_zlantr_std;
    test_zlantr.next   = NULL;

    testing_register( &test_zlantr );
}
