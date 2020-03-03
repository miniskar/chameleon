/**
 *
 * @file core_ztile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CPU kernel interface from CHAM_tile_t layout to the real one.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/coreblas_ztile.h"

void
TCORE_dzasum( cham_store_t       storev,
              cham_uplo_t        uplo,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_dzasum( storev, uplo, M, N, A->mat, A->ld, work );
}

int
TCORE_zaxpy( int                   M,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             int                   incA,
             CHAM_tile_t *         B,
             int                   incB )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zaxpy( M, alpha, A->mat, incA, B->mat, incB );
}

int
TCORE_zgeadd( cham_trans_t          trans,
              int                   M,
              int                   N,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         B )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgeadd( trans, M, N, alpha, A->mat, A->ld, beta, B->mat, B->ld );
}

int
TCORE_zgelqt( int                    M,
              int                    N,
              int                    IB,
              CHAM_tile_t *          A,
              CHAM_tile_t *          T,
              CHAMELEON_Complex64_t *TAU,
              CHAMELEON_Complex64_t *WORK )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgelqt( M, N, IB, A->mat, A->ld, T->mat, T->ld, TAU, WORK );
}

void
TCORE_zgemm( cham_trans_t          transA,
             cham_trans_t          transB,
             int                   M,
             int                   N,
             int                   K,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             const CHAM_tile_t *   B,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    CORE_zgemm(
        transA, transB, M, N, K, alpha, A->mat, A->ld, B->mat, B->ld, beta, C->mat, C->ld );
}

int
TCORE_zgeqrt( int                    M,
              int                    N,
              int                    IB,
              CHAM_tile_t *          A,
              CHAM_tile_t *          T,
              CHAMELEON_Complex64_t *TAU,
              CHAMELEON_Complex64_t *WORK )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgeqrt( M, N, IB, A->mat, A->ld, T->mat, T->ld, TAU, WORK );
}

int
TCORE_zgessm( int M, int N, int K, int IB, const int *IPIV, const CHAM_tile_t *L, CHAM_tile_t *A )
{
    assert( L->format & CHAMELEON_TILE_FULLRANK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgessm( M, N, K, IB, IPIV, L->mat, L->ld, A->mat, A->ld );
}

int
TCORE_zgessq( cham_store_t storev, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sclssq )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( sclssq->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgessq( storev, M, N, A->mat, A->ld, sclssq->mat );
}

int
TCORE_zgetrf( int M, int N, CHAM_tile_t *A, int *IPIV, int *INFO )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgetrf( M, N, A->mat, A->ld, IPIV, INFO );
}

int
TCORE_zgetrf_incpiv( int M, int N, int IB, CHAM_tile_t *A, int *IPIV, int *INFO )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgetrf_incpiv( M, N, IB, A->mat, A->ld, IPIV, INFO );
}

int
TCORE_zgetrf_nopiv( int M, int N, int IB, CHAM_tile_t *A, int *INFO )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgetrf_nopiv( M, N, IB, A->mat, A->ld, INFO );
}

void
TCORE_zhe2ge( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    CORE_zhe2ge( uplo, M, N, A->mat, A->ld, B->mat, B->ld );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_zhemm( cham_side_t           side,
             cham_uplo_t           uplo,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             const CHAM_tile_t *   B,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    CORE_zhemm( side, uplo, M, N, alpha, A->mat, A->ld, B->mat, B->ld, beta, C->mat, C->ld );
}

void
TCORE_zherk( cham_uplo_t        uplo,
             cham_trans_t       trans,
             int                N,
             int                K,
             double             alpha,
             const CHAM_tile_t *A,
             double             beta,
             CHAM_tile_t *      C )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    CORE_zherk( uplo, trans, N, K, alpha, A->mat, A->ld, beta, C->mat, C->ld );
}

void
TCORE_zher2k( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   N,
              int                   K,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              const CHAM_tile_t *   B,
              double                beta,
              CHAM_tile_t *         C )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    CORE_zher2k( uplo, trans, N, K, alpha, A->mat, A->ld, B->mat, B->ld, beta, C->mat, C->ld );
}
#endif

int
TCORE_zherfb( cham_uplo_t            uplo,
              int                    N,
              int                    K,
              int                    IB,
              int                    NB,
              const CHAM_tile_t *    A,
              const CHAM_tile_t *    T,
              CHAM_tile_t *          C,
              CHAMELEON_Complex64_t *WORK,
              int                    ldwork )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zherfb(
        uplo, N, K, IB, NB, A->mat, A->ld, T->mat, T->ld, C->mat, C->ld, WORK, ldwork );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
int
TCORE_zhessq( cham_store_t       storev,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      sclssq )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( sclssq->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zhessq( storev, uplo, N, A->mat, A->ld, sclssq->mat );
}
#endif

void
TCORE_zlacpy( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlacpy( uplo, M, N, A->mat, A->ld, B->mat, B->ld );
}

void
TCORE_zlange( cham_normtype_t    norm,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlange( norm, M, N, A->mat, A->ld, work, normA );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_zlanhe( cham_normtype_t    norm,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlanhe( norm, uplo, N, A->mat, A->ld, work, normA );
}
#endif

void
TCORE_zlansy( cham_normtype_t    norm,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlansy( norm, uplo, N, A->mat, A->ld, work, normA );
}

void
TCORE_zlantr( cham_normtype_t    norm,
              cham_uplo_t        uplo,
              cham_diag_t        diag,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlantr( norm, uplo, diag, M, N, A->mat, A->ld, work, normA );
}

int
TCORE_zlascal( cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zlascal( uplo, m, n, alpha, A->mat, A->ld );
}

void
TCORE_zlaset( cham_uplo_t           uplo,
              int                   n1,
              int                   n2,
              CHAMELEON_Complex64_t alpha,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         A )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlaset( uplo, n1, n2, alpha, beta, A->mat, A->ld );
}

void
TCORE_zlaset2( cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlaset2( uplo, n1, n2, alpha, A->mat, A->ld );
}

int
TCORE_zlatro( cham_uplo_t        uplo,
              cham_trans_t       trans,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      B )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zlatro( uplo, trans, M, N, A->mat, A->ld, B->mat, B->ld );
}

void
TCORE_zlauum( cham_uplo_t uplo, int N, CHAM_tile_t *A )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zlauum( uplo, N, A->mat, A->ld );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_zplghe( double                 bump,
              int                    m,
              int                    n,
              CHAM_tile_t *          tileA,
              int                    bigM,
              int                    m0,
              int                    n0,
              unsigned long long int seed )
{
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    CORE_zplghe( bump, m, n, tileA->mat, tileA->ld, bigM, m0, n0, seed );
}
#endif

void
TCORE_zplgsy( CHAMELEON_Complex64_t  bump,
              int                    m,
              int                    n,
              CHAM_tile_t *          tileA,
              int                    bigM,
              int                    m0,
              int                    n0,
              unsigned long long int seed )
{
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    CORE_zplgsy( bump, m, n, tileA->mat, tileA->ld, bigM, m0, n0, seed );
}

void
TCORE_zplrnt( int                    m,
              int                    n,
              CHAM_tile_t *          tileA,
              int                    bigM,
              int                    m0,
              int                    n0,
              unsigned long long int seed )
{
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    CORE_zplrnt( m, n, tileA->mat, tileA->ld, bigM, m0, n0, seed );
}

void
TCORE_zpotrf( cham_uplo_t uplo, int n, CHAM_tile_t *A, int *INFO )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_zpotrf( uplo, n, A->mat, A->ld, INFO );
}

int
TCORE_zssssm( int                M1,
              int                N1,
              int                M2,
              int                N2,
              int                K,
              int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *L1,
              const CHAM_tile_t *L2,
              const int *        IPIV )
{
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( L1->format & CHAMELEON_TILE_FULLRANK );
    assert( L2->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zssssm( M1,
                        N1,
                        M2,
                        N2,
                        K,
                        IB,
                        A1->mat,
                        A1->ld,
                        A2->mat,
                        A2->ld,
                        L1->mat,
                        L1->ld,
                        L2->mat,
                        L2->ld,
                        IPIV );
}

void
TCORE_zsymm( cham_side_t           side,
             cham_uplo_t           uplo,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             const CHAM_tile_t *   B,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    CORE_zsymm( side, uplo, M, N, alpha, A->mat, A->ld, B->mat, B->ld, beta, C->mat, C->ld );
}

void
TCORE_zsyrk( cham_uplo_t           uplo,
             cham_trans_t          trans,
             int                   N,
             int                   K,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    CORE_zsyrk( uplo, trans, N, K, alpha, A->mat, A->ld, beta, C->mat, C->ld );
}

void
TCORE_zsyr2k( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   N,
              int                   K,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              const CHAM_tile_t *   B,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         C )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    CORE_zsyr2k( uplo, trans, N, K, alpha, A->mat, A->ld, B->mat, B->ld, beta, C->mat, C->ld );
}

int
TCORE_zsyssq( cham_store_t       storev,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      sclssq )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( sclssq->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zsyssq( storev, uplo, N, A->mat, A->ld, sclssq->mat );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
int
TCORE_zsytf2_nopiv( cham_uplo_t uplo, int n, CHAM_tile_t *A )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zsytf2_nopiv( uplo, n, A->mat, A->ld );
}
#endif

int
TCORE_ztplqt( int                    M,
              int                    N,
              int                    L,
              int                    IB,
              CHAM_tile_t *          A,
              CHAM_tile_t *          B,
              CHAM_tile_t *          T,
              CHAMELEON_Complex64_t *WORK )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztplqt( M, N, L, IB, A->mat, A->ld, B->mat, B->ld, T->mat, T->ld, WORK );
}

int
TCORE_ztpmlqt( cham_side_t            side,
               cham_trans_t           trans,
               int                    M,
               int                    N,
               int                    K,
               int                    L,
               int                    IB,
               const CHAM_tile_t *    V,
               const CHAM_tile_t *    T,
               CHAM_tile_t *          A,
               CHAM_tile_t *          B,
               CHAMELEON_Complex64_t *WORK )
{
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztpmlqt( side,
                         trans,
                         M,
                         N,
                         K,
                         L,
                         IB,
                         V->mat,
                         V->ld,
                         T->mat,
                         T->ld,
                         A->mat,
                         A->ld,
                         B->mat,
                         B->ld,
                         WORK );
}

int
TCORE_ztpmqrt( cham_side_t            side,
               cham_trans_t           trans,
               int                    M,
               int                    N,
               int                    K,
               int                    L,
               int                    IB,
               const CHAM_tile_t *    V,
               const CHAM_tile_t *    T,
               CHAM_tile_t *          A,
               CHAM_tile_t *          B,
               CHAMELEON_Complex64_t *WORK )
{
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztpmqrt( side,
                         trans,
                         M,
                         N,
                         K,
                         L,
                         IB,
                         V->mat,
                         V->ld,
                         T->mat,
                         T->ld,
                         A->mat,
                         A->ld,
                         B->mat,
                         B->ld,
                         WORK );
}

int
TCORE_ztpqrt( int                    M,
              int                    N,
              int                    L,
              int                    IB,
              CHAM_tile_t *          A,
              CHAM_tile_t *          B,
              CHAM_tile_t *          T,
              CHAMELEON_Complex64_t *WORK )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztpqrt( M, N, L, IB, A->mat, A->ld, B->mat, B->ld, T->mat, T->ld, WORK );
}

int
TCORE_ztradd( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   M,
              int                   N,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         B )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztradd( uplo, trans, M, N, alpha, A->mat, A->ld, beta, B->mat, B->ld );
}

void
TCORE_ztrasm( cham_store_t       storev,
              cham_uplo_t        uplo,
              cham_diag_t        diag,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_ztrasm( storev, uplo, diag, M, N, A->mat, A->ld, work );
}

void
TCORE_ztrmm( cham_side_t           side,
             cham_uplo_t           uplo,
             cham_trans_t          transA,
             cham_diag_t           diag,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             CHAM_tile_t *         B )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    CORE_ztrmm( side, uplo, transA, diag, M, N, alpha, A->mat, A->ld, B->mat, B->ld );
}

void
TCORE_ztrsm( cham_side_t           side,
             cham_uplo_t           uplo,
             cham_trans_t          transA,
             cham_diag_t           diag,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             CHAM_tile_t *         B )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    CORE_ztrsm( side, uplo, transA, diag, M, N, alpha, A->mat, A->ld, B->mat, B->ld );
}

int
TCORE_ztrssq( cham_uplo_t        uplo,
              cham_diag_t        diag,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      sclssq )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( sclssq->format & CHAMELEON_TILE_FULLRANK );
    double *W = sclssq->mat;
    return CORE_ztrssq( uplo, diag, M, N, A->mat, A->ld, W, W + 1 );
}

void
TCORE_ztrtri( cham_uplo_t uplo, cham_diag_t diag, int N, CHAM_tile_t *A, int *info )
{
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    CORE_ztrtri( uplo, diag, N, A->mat, A->ld, info );
}

int
TCORE_ztsmlq_hetra1( cham_side_t            side,
                     cham_trans_t           trans,
                     int                    m1,
                     int                    n1,
                     int                    m2,
                     int                    n2,
                     int                    k,
                     int                    ib,
                     CHAM_tile_t *          A1,
                     CHAM_tile_t *          A2,
                     const CHAM_tile_t *    V,
                     const CHAM_tile_t *    T,
                     CHAMELEON_Complex64_t *WORK,
                     int                    ldwork )
{
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztsmlq_hetra1( side,
                               trans,
                               m1,
                               n1,
                               m2,
                               n2,
                               k,
                               ib,
                               A1->mat,
                               A1->ld,
                               A2->mat,
                               A2->ld,
                               V->mat,
                               V->ld,
                               T->mat,
                               T->ld,
                               WORK,
                               ldwork );
}

int
TCORE_ztsmqr_hetra1( cham_side_t            side,
                     cham_trans_t           trans,
                     int                    m1,
                     int                    n1,
                     int                    m2,
                     int                    n2,
                     int                    k,
                     int                    ib,
                     CHAM_tile_t *          A1,
                     CHAM_tile_t *          A2,
                     const CHAM_tile_t *    V,
                     const CHAM_tile_t *    T,
                     CHAMELEON_Complex64_t *WORK,
                     int                    ldwork )
{
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztsmqr_hetra1( side,
                               trans,
                               m1,
                               n1,
                               m2,
                               n2,
                               k,
                               ib,
                               A1->mat,
                               A1->ld,
                               A2->mat,
                               A2->ld,
                               V->mat,
                               V->ld,
                               T->mat,
                               T->ld,
                               WORK,
                               ldwork );
}

int
TCORE_ztstrf( int                    M,
              int                    N,
              int                    IB,
              int                    NB,
              CHAM_tile_t *          U,
              CHAM_tile_t *          A,
              CHAM_tile_t *          L,
              int *                  IPIV,
              CHAMELEON_Complex64_t *WORK,
              int                    LDWORK,
              int *                  INFO )
{
    assert( U->format & CHAMELEON_TILE_FULLRANK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( L->format & CHAMELEON_TILE_FULLRANK );
    return CORE_ztstrf(
        M, N, IB, NB, U->mat, U->ld, A->mat, A->ld, L->mat, L->ld, IPIV, WORK, LDWORK, INFO );
}

int
TCORE_zunmlq( cham_side_t            side,
              cham_trans_t           trans,
              int                    M,
              int                    N,
              int                    K,
              int                    IB,
              const CHAM_tile_t *    V,
              const CHAM_tile_t *    T,
              CHAM_tile_t *          C,
              CHAMELEON_Complex64_t *WORK,
              int                    LDWORK )
{
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zunmlq(
        side, trans, M, N, K, IB, V->mat, V->ld, T->mat, T->ld, C->mat, C->ld, WORK, LDWORK );
}

int
TCORE_zunmqr( cham_side_t            side,
              cham_trans_t           trans,
              int                    M,
              int                    N,
              int                    K,
              int                    IB,
              const CHAM_tile_t *    V,
              const CHAM_tile_t *    T,
              CHAM_tile_t *          C,
              CHAMELEON_Complex64_t *WORK,
              int                    LDWORK )
{
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zunmqr(
        side, trans, M, N, K, IB, V->mat, V->ld, T->mat, T->ld, C->mat, C->ld, WORK, LDWORK );
}

int
TCORE_zgram( cham_uplo_t        uplo,
             int                M,
             int                N,
             int                Mt,
             int                Nt,
             const CHAM_tile_t *Di,
             const CHAM_tile_t *Dj,
             const CHAM_tile_t *D,
             CHAM_tile_t *      A )
{
    assert( Di->format & CHAMELEON_TILE_FULLRANK );
    assert( Dj->format & CHAMELEON_TILE_FULLRANK );
    assert( D->format & CHAMELEON_TILE_FULLRANK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CORE_zgram(
        uplo, M, N, Mt, Nt, Di->mat, Di->ld, Dj->mat, Dj->ld, D->mat, A->mat, A->ld );
}
