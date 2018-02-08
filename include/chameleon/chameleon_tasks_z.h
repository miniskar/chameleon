/**
 *
 * @file chameleon_tasks_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t elementary tasks header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_tasks_z_h_
#define _chameleon_tasks_z_h_

/**
 *  Declarations of QUARK wrappers (called by CHAMELEON) - alphabetical order
 */
void InsertTask_dzasum( const RUNTIME_option_t *options,
                        MORSE_enum storev, MORSE_enum uplo, int M, int N,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn );
void InsertTask_zgeadd( const RUNTIME_option_t *options,
                        MORSE_enum trans, int m, int n, int nb,
                        MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_Complex64_t beta,  const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_zlascal( const RUNTIME_option_t *options,
                         MORSE_enum uplo,
                         int m, int n, int nb,
                         MORSE_Complex64_t alpha,
                         const MORSE_desc_t *A, int Am, int An, int lda );
void InsertTask_zbrdalg( const RUNTIME_option_t *options,
                         MORSE_enum uplo,
                         int N, int NB,
                         const MORSE_desc_t *A,
                         const MORSE_desc_t *C, int Cm, int Cn,
                         const MORSE_desc_t *S, int Sm, int Sn,
                         int i, int j, int m, int grsiz, int BAND,
                         int *PCOL, int *ACOL, int *MCOL );
void InsertTask_zgelqt( const RUNTIME_option_t *options,
                        int m, int n, int ib, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_zgemm( const RUNTIME_option_t *options,
                       MORSE_enum transA, MORSE_enum transB,
                       int m, int n, int k, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                       MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zgemm2( const RUNTIME_option_t *options,
                        MORSE_enum transA, MORSE_enum transB,
                        int m, int n, int k, int nb,
                        MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                        MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zgemm_f2( const RUNTIME_option_t *options,
                          MORSE_enum transA, MORSE_enum transB,
                          int m, int n, int k, int nb,
                          MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                          const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                          MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc,
                          const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                          const MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2 );
void InsertTask_zgemm_p2( const RUNTIME_option_t *options,
                          MORSE_enum transA, MORSE_enum transB,
                          int m, int n, int k, int nb,
                          MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                          const MORSE_Complex64_t **B, int ldb,
                          MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zgemm_p2f1( const RUNTIME_option_t *options,
                            MORSE_enum transA, MORSE_enum transB,
                            int m, int n, int k, int nb,
                            MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                            const MORSE_Complex64_t **B, int ldb,
                            MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc,
                            const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1 );
void InsertTask_zgemm_p3( const RUNTIME_option_t *options,
                          MORSE_enum transA, MORSE_enum transB,
                          int m, int n, int k, int nb,
                          MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                          const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                          MORSE_Complex64_t beta, MORSE_Complex64_t **C, int ldc );
void InsertTask_zgeqrt( const RUNTIME_option_t *options,
                        int m, int n, int ib, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_zgessm( const RUNTIME_option_t *options,
                        int m, int n, int k, int ib, int nb,
                        int *IPIV,
                        const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                        const MORSE_desc_t *D, int Dm, int Dn, int ldd,
                        const MORSE_desc_t *A, int Am, int An, int lda );
void InsertTask_zgessq( const RUNTIME_option_t *options,
                        int m, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void InsertTask_zgetrf( const RUNTIME_option_t *options,
                        int m, int n, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        int *IPIV,

                        MORSE_bool check_info, int iinfo );
void InsertTask_zgetrf_incpiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const MORSE_desc_t *A, int Am, int An, int lda,
                               const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                               int *IPIV,
                               MORSE_bool check_info, int iinfo );
void InsertTask_zgetrf_nopiv( const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const MORSE_desc_t *A, int Am, int An, int lda, int iinfo );
void InsertTask_zgetrf_reclap( const RUNTIME_option_t *options,
                               int m, int n, int nb,
                               const MORSE_desc_t *A, int Am, int An, int lda,
                               int *IPIV,

                               MORSE_bool check_info, int iinfo,
                               int nbthread );
void InsertTask_zgetrf_rectil( const RUNTIME_option_t *options,
                               const MORSE_desc_t A, const MORSE_desc_t *Amn, int Amnm, int Amnn, int size,
                               int *IPIV,

                               MORSE_bool check_info, int iinfo,
                               int nbthread );
void InsertTask_zgetrip( const RUNTIME_option_t *options,
                         int m, int n, const MORSE_desc_t *A, int Am, int An, int szeA );
void InsertTask_zgetrip_f1( const RUNTIME_option_t *options,
                            int m, int n, const MORSE_desc_t *A, int Am, int An, int szeA,
                            const MORSE_desc_t *fake, int fakem, int faken, int szeF, int paramF );
void InsertTask_zgetrip_f2( const RUNTIME_option_t *options,
                            int m, int n, const MORSE_desc_t *A, int Am, int An, int szeA,
                            const MORSE_desc_t *fake1, int fake1m, int fake1n, int szeF1, int paramF1,
                            const MORSE_desc_t *fake2, int fake2m, int fake2n, int szeF2, int paramF2 );
void InsertTask_zhe2ge( const RUNTIME_option_t *options,
                        MORSE_enum uplo,
                        int m, int n, int mb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_zhemm( const RUNTIME_option_t *options,
                       MORSE_enum side, MORSE_enum uplo,
                       int m, int n, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                       MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zhegst( const RUNTIME_option_t *options,
                        int itype, MORSE_enum uplo, int N,
                        const MORSE_desc_t *A, int Am, int An, int LDA,
                        const MORSE_desc_t *B, int Bm, int Bn, int LDB,
                        int iinfo );
void InsertTask_zherk( const RUNTIME_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       double alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       double beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zher2k( const RUNTIME_option_t *options,
                        MORSE_enum uplo, MORSE_enum trans,
                        int n, int k, int nb,
                        MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int LDB,
                        double beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zherfb( const RUNTIME_option_t *options,
                        MORSE_enum uplo,
                        int n, int k, int ib, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                        const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zlacpy( const RUNTIME_option_t *options,
                        MORSE_enum uplo, int m, int n, int mb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_zlacpyx( const RUNTIME_option_t *options,
                         MORSE_enum uplo, int m, int n, int mb,
                         int displA, const MORSE_desc_t *A, int Am, int An, int lda,
                         int displB, const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_zlange( const RUNTIME_option_t *options,
                        MORSE_enum norm, int M, int N, int NB,
                        const MORSE_desc_t *A, int Am, int An, int LDA,
                        const MORSE_desc_t *B, int Bm, int Bn );
void InsertTask_zlange_max( const RUNTIME_option_t *options,
                            const MORSE_desc_t *A, int Am, int An,
                            const MORSE_desc_t *B, int Bm, int Bn );
void InsertTask_zhessq( const RUNTIME_option_t *options,
                        MORSE_enum uplo, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void InsertTask_zlanhe( const RUNTIME_option_t *options,
                        MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                        const MORSE_desc_t *A, int Am, int An, int LDA,
                        const MORSE_desc_t *B, int Bm, int Bn );
void InsertTask_zlansy( const RUNTIME_option_t *options,
                        MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                        const MORSE_desc_t *A, int Am, int An, int LDA,
                        const MORSE_desc_t *B, int Bm, int Bn );
void InsertTask_zlantr( const RUNTIME_option_t *options,
                        MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag,
                        int M, int N, int NB,
                        const MORSE_desc_t *A, int Am, int An, int LDA,
                        const MORSE_desc_t *B, int Bm, int Bn );
void InsertTask_zlaset( const RUNTIME_option_t *options,
                        MORSE_enum uplo, int n1, int n2, MORSE_Complex64_t alpha,
                        MORSE_Complex64_t beta, const MORSE_desc_t *tileA, int tileAm, int tileAn, int ldtilea );
void InsertTask_zlaset2( const RUNTIME_option_t *options,
                         MORSE_enum uplo, int n1, int n2, MORSE_Complex64_t alpha,
                         const MORSE_desc_t *tileA, int tileAm, int tileAn, int ldtilea );
void InsertTask_zlaswp( const RUNTIME_option_t *options,
                        int n, const MORSE_desc_t *A, int Am, int An, int lda,
                        int i1,  int i2, int *ipiv, int inc );
void InsertTask_zlaswp_f2( const RUNTIME_option_t *options,
                           int n, const MORSE_desc_t *A, int Am, int An, int lda,
                           int i1,  int i2, int *ipiv, int inc,
                           const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                           const MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2 );
void InsertTask_zlaswp_ontile( const RUNTIME_option_t *options,
                               const MORSE_desc_t descA, const MORSE_desc_t *A, int Am, int An,
                               int i1,  int i2, int *ipiv, int inc, MORSE_Complex64_t *fakepanel );
void InsertTask_zlaswp_ontile_f2( const RUNTIME_option_t *options,
                                  const MORSE_desc_t descA, const MORSE_desc_t *A, int Am, int An,
                                  int i1,  int i2, int *ipiv, int inc,
                                  const MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                                  const MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2 );
void InsertTask_zlaswpc_ontile( const RUNTIME_option_t *options,
                                const MORSE_desc_t descA, const MORSE_desc_t *A, int Am, int An,
                                int i1,  int i2, int *ipiv, int inc, MORSE_Complex64_t *fakepanel );
void InsertTask_zlatro( const RUNTIME_option_t *options,
                        MORSE_enum uplo, MORSE_enum trans, int m, int n, int mb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_zlauum( const RUNTIME_option_t *options,
                        MORSE_enum uplo, int n, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda );
void InsertTask_zplghe( const RUNTIME_option_t *options,
                        double bump, int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed );
void InsertTask_zplgsy( const RUNTIME_option_t *options,
                        MORSE_Complex64_t bump, int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed );
void InsertTask_zplrnt( const RUNTIME_option_t *options,
                        int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed );
void InsertTask_zpotrf( const RUNTIME_option_t *options,
                        MORSE_enum uplo, int n, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,

                        int iinfo );
void InsertTask_zshift( const RUNTIME_option_t *options,
                        int s, int m, int n, int L,
                        MORSE_Complex64_t *A );
void InsertTask_zshiftw( const RUNTIME_option_t *options,
                         int s, int cl, int m, int n, int L,
                         const MORSE_desc_t *A, int Am, int An, MORSE_Complex64_t *W );
void InsertTask_zssssm( const RUNTIME_option_t *options,
                        int m1, int n1, int m2, int n2, int k, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *L1, int L1m, int L1n, int ldl1,
                        const MORSE_desc_t *L2, int L2m, int L2n, int ldl2,
                        const int *IPIV );
void InsertTask_zsymm( const RUNTIME_option_t *options,
                       MORSE_enum side, MORSE_enum uplo,
                       int m, int n, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                       MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zsyrk( const RUNTIME_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zsyr2k( const RUNTIME_option_t *options,
                        MORSE_enum uplo, MORSE_enum trans,
                        int n, int k, int nb,
                        MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int LDB,
                        MORSE_Complex64_t beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zsyssq( const RUNTIME_option_t *options,
                        MORSE_enum uplo, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void InsertTask_zsytrf_nopiv( const RUNTIME_option_t *options,
                              MORSE_enum uplo, int n, int nb,
                              const MORSE_desc_t *A, int Am, int An, int lda,
                              int iinfo );
void InsertTask_zswpab( const RUNTIME_option_t *options,
                        int i, int n1, int n2,
                        const MORSE_desc_t *A, int Am, int An, int szeA );
void InsertTask_zswptr_ontile( const RUNTIME_option_t *options,
                               const MORSE_desc_t descA, const MORSE_desc_t *Aij, int Aijm, int Aijn,
                               int i1,  int i2, int *ipiv, int inc,
                               const MORSE_desc_t *Akk, int Akkm, int Akkn, int ldak );
void InsertTask_ztplqt( const RUNTIME_option_t *options,
                        int m, int n, int l, int ib, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztpmlqt( const RUNTIME_option_t *options,
                         MORSE_enum side, MORSE_enum trans,
                         int M, int N, int K, int L, int ib, int nb,
                         const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                         const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                         const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_ztpmqrt( const RUNTIME_option_t *options,
                         MORSE_enum side, MORSE_enum trans,
                         int m, int n, int k, int l, int ib, int nb,
                         const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                         const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                         const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_ztpqrt( const RUNTIME_option_t *options,
                        int m, int n, int l, int ib, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztrdalg( const RUNTIME_option_t *options,
                         MORSE_enum uplo,
                         int N, int NB,
                         const MORSE_desc_t *A,
                         const MORSE_desc_t *C, int Cm, int Cn,
                         const MORSE_desc_t *S, int Sm, int Sn,
                         int i, int j, int m, int grsiz, int BAND,
                         int *PCOL, int *ACOL, int *MCOL );
void InsertTask_ztradd( const RUNTIME_option_t *options,
                        MORSE_enum uplo, MORSE_enum trans, int m, int n, int nb,
                        MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_Complex64_t beta,  const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_ztrasm( const RUNTIME_option_t *options,
                        MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *B, int Bm, int Bn );
void InsertTask_ztrmm( const RUNTIME_option_t *options,
                       MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                       int m, int n, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_ztrmm_p2( const RUNTIME_option_t *options,
                          MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                          int m, int n, int nb,
                          MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                          MORSE_Complex64_t **B, int ldb );
void InsertTask_ztrsm( const RUNTIME_option_t *options,
                       MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                       int m, int n, int nb,
                       MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb );
void InsertTask_ztrssq( const RUNTIME_option_t *options,
                        MORSE_enum uplo, MORSE_enum diag,
                        int m, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void InsertTask_ztrtri( const RUNTIME_option_t *options,
                        MORSE_enum uplo, MORSE_enum diag, int n, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,

                        int iinfo );
void InsertTask_ztslqt( const RUNTIME_option_t *options,
                        int m, int n, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztsmlq( const RUNTIME_option_t *options,
                        MORSE_enum side, MORSE_enum trans,
                        int m1, int n1, int m2, int n2, int k, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztsmlq_hetra1( const RUNTIME_option_t *options,
                               MORSE_enum side, MORSE_enum trans,
                               int m1, int n1, int m2, int n2, int k, int ib, int nb,
                               const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                               const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                               const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                               const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztsmqr( const RUNTIME_option_t *options,
                        MORSE_enum side, MORSE_enum trans,
                        int m1, int n1, int m2, int n2, int k, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztsmqr_hetra1( const RUNTIME_option_t *options,
                               MORSE_enum side, MORSE_enum trans,
                               int m1, int n1, int m2, int n2, int k, int ib, int nb,
                               const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                               const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                               const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                               const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztsqrt( const RUNTIME_option_t *options,
                        int m, int n, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_ztstrf( const RUNTIME_option_t *options,
                        int m, int n, int ib, int nb,
                        const MORSE_desc_t *U, int Um, int Un, int ldu,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                        int *IPIV,

                        MORSE_bool check_info, int iinfo );
void InsertTask_zttmqr( const RUNTIME_option_t *options,
                        MORSE_enum side, MORSE_enum trans,
                        int m1, int n1, int m2, int n2, int k, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_zttqrt( const RUNTIME_option_t *options,
                        int m, int n, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_zttmlq( const RUNTIME_option_t *options,
                        MORSE_enum side, MORSE_enum trans,
                        int m1, int n1, int m2, int n2, int k, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_zttlqt( const RUNTIME_option_t *options,
                        int m, int n, int ib, int nb,
                        const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                        const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt );
void InsertTask_zpamm( const RUNTIME_option_t *options,
                       int op, MORSE_enum side, MORSE_enum storev,
                       int m, int n, int k, int l,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       const MORSE_desc_t *W, int Wm, int Wn, int ldw );
void InsertTask_zplssq( const RUNTIME_option_t *options,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                        const MORSE_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn );
void InsertTask_zplssq2( const RUNTIME_option_t *options,
                         const MORSE_desc_t *RESULT, int RESULTm, int RESULTn );
void InsertTask_zunmlq( const RUNTIME_option_t *options,
                        MORSE_enum side, MORSE_enum trans,
                        int m, int n, int ib,  int nb, int k,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                        const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zunmqr( const RUNTIME_option_t *options,
                        MORSE_enum side, MORSE_enum trans,
                        int m, int n, int k, int ib, int nb,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                        const MORSE_desc_t *C, int Cm, int Cn, int ldc );
void InsertTask_zbuild( const RUNTIME_option_t *options,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        void *user_data, void* user_build_callback );

#endif
