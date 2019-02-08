/**
 *
 * @file chameleon_tasks_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2018-11-07
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_tasks_z_h_
#define _chameleon_tasks_z_h_

/**
 *  Declarations of QUARK wrappers (called by CHAMELEON) - alphabetical order
 */
void INSERT_TASK_dzasum( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int M, int N,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zaxpy( const RUNTIME_option_t *options,
                        int M, CHAMELEON_Complex64_t alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB );
void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_zlascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *A, int Am, int An, int lda );
void INSERT_TASK_zbrdalg( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int N, int NB,
                          const CHAM_desc_t *A,
                          const CHAM_desc_t *C, int Cm, int Cn,
                          const CHAM_desc_t *S, int Sm, int Sn,
                          int i, int j, int m, int grsiz, int BAND,
                          int *PCOL, int *ACOL, int *MCOL );
void INSERT_TASK_zgelqt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt );
void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zgemm2( const RUNTIME_option_t *options,
                         cham_trans_t transA, cham_trans_t transB,
                         int m, int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                         CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zgemm_f2( const RUNTIME_option_t *options,
                           cham_trans_t transA, cham_trans_t transB,
                           int m, int n, int k, int nb,
                           CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                           const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                           CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc,
                           const CHAM_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                           const CHAM_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2 );
void INSERT_TASK_zgemm_p2( const RUNTIME_option_t *options,
                           cham_trans_t transA, cham_trans_t transB,
                           int m, int n, int k, int nb,
                           CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                           const CHAMELEON_Complex64_t **B, int ldb,
                           CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zgemm_p2f1( const RUNTIME_option_t *options,
                             cham_trans_t transA, cham_trans_t transB,
                             int m, int n, int k, int nb,
                             CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                             const CHAMELEON_Complex64_t **B, int ldb,
                             CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc,
                             const CHAM_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1 );
void INSERT_TASK_zgemm_p3( const RUNTIME_option_t *options,
                           cham_trans_t transA, cham_trans_t transB,
                           int m, int n, int k, int nb,
                           CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                           const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                           CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t **C, int ldc );
void INSERT_TASK_zgeqrt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt );
void INSERT_TASK_zgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln, int ldl,
                         const CHAM_desc_t *D, int Dm, int Dn, int ldd,
                         const CHAM_desc_t *A, int Am, int An, int lda );
void INSERT_TASK_zgessq( const RUNTIME_option_t *options,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo );
void INSERT_TASK_zgetrf_incpiv( const RUNTIME_option_t *options,
                                int m, int n, int ib, int nb,
                                const CHAM_desc_t *A, int Am, int An, int lda,
                                const CHAM_desc_t *L, int Lm, int Ln, int ldl,
                                int *IPIV,
                                cham_bool_t check_info, int iinfo );
void INSERT_TASK_zgetrf_nopiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int lda, int iinfo );
void INSERT_TASK_zgetrf_reclap( const RUNTIME_option_t *options,
                                int m, int n, int nb,
                                const CHAM_desc_t *A, int Am, int An, int lda,
                                int *IPIV,

                                cham_bool_t check_info, int iinfo,
                                int nbthread );
void INSERT_TASK_zgetrf_rectil( const RUNTIME_option_t *options,
                                const CHAM_desc_t A, const CHAM_desc_t *Amn, int Amnm, int Amnn, int size,
                                int *IPIV,

                                cham_bool_t check_info, int iinfo,
                                int nbthread );
void INSERT_TASK_zgetrip( const RUNTIME_option_t *options,
                          int m, int n, const CHAM_desc_t *A, int Am, int An, int szeA );
void INSERT_TASK_zgetrip_f1( const RUNTIME_option_t *options,
                             int m, int n, const CHAM_desc_t *A, int Am, int An, int szeA,
                             const CHAM_desc_t *fake, int fakem, int faken, int szeF, int paramF );
void INSERT_TASK_zgetrip_f2( const RUNTIME_option_t *options,
                             int m, int n, const CHAM_desc_t *A, int Am, int An, int szeA,
                             const CHAM_desc_t *fake1, int fake1m, int fake1n, int szeF1, int paramF1,
                             const CHAM_desc_t *fake2, int fake2m, int fake2n, int szeF2, int paramF2 );
void INSERT_TASK_zhe2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_zhemm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zhegst( const RUNTIME_option_t *options,
                         int itype, cham_uplo_t uplo, int N,
                         const CHAM_desc_t *A, int Am, int An, int LDA,
                         const CHAM_desc_t *B, int Bm, int Bn, int LDB,
                         int iinfo );
void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                        double beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zher2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int LDB,
                         double beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zherfb( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                         const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int mb,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_zlange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An, int LDA,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlange_max( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zhessq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zlanhe( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An, int LDA,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlansy( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An, int LDA,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An, int LDA,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha,
                         CHAMELEON_Complex64_t beta, const CHAM_desc_t *tileA, int tileAm, int tileAn, int ldtilea );
void INSERT_TASK_zlaset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *tileA, int tileAm, int tileAn, int ldtilea );
void INSERT_TASK_zlaswp( const RUNTIME_option_t *options,
                         int n, const CHAM_desc_t *A, int Am, int An, int lda,
                         int i1,  int i2, int *ipiv, int inc );
void INSERT_TASK_zlaswp_f2( const RUNTIME_option_t *options,
                            int n, const CHAM_desc_t *A, int Am, int An, int lda,
                            int i1,  int i2, int *ipiv, int inc,
                            const CHAM_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                            const CHAM_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2 );
void INSERT_TASK_zlaswp_ontile( const RUNTIME_option_t *options,
                                const CHAM_desc_t descA, const CHAM_desc_t *A, int Am, int An,
                                int i1,  int i2, int *ipiv, int inc, CHAMELEON_Complex64_t *fakepanel );
void INSERT_TASK_zlaswp_ontile_f2( const RUNTIME_option_t *options,
                                   const CHAM_desc_t descA, const CHAM_desc_t *A, int Am, int An,
                                   int i1,  int i2, int *ipiv, int inc,
                                   const CHAM_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                                   const CHAM_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2 );
void INSERT_TASK_zlaswpc_ontile( const RUNTIME_option_t *options,
                                 const CHAM_desc_t descA, const CHAM_desc_t *A, int Am, int An,
                                 int i1,  int i2, int *ipiv, int inc, CHAMELEON_Complex64_t *fakepanel );
void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_zlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda );
void INSERT_TASK_zplghe( const RUNTIME_option_t *options,
                         double bump, int m, int n, const CHAM_desc_t *A, int Am, int An, int lda,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An, int lda,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplrnt( const RUNTIME_option_t *options,
                         int m, int n, const CHAM_desc_t *A, int Am, int An, int lda,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,

                         int iinfo );
void INSERT_TASK_zshift( const RUNTIME_option_t *options,
                         int s, int m, int n, int L,
                         CHAMELEON_Complex64_t *A );
void INSERT_TASK_zshiftw( const RUNTIME_option_t *options,
                          int s, int cl, int m, int n, int L,
                          const CHAM_desc_t *A, int Am, int An, CHAMELEON_Complex64_t *W );
void INSERT_TASK_zssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                         const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                         const CHAM_desc_t *L1, int L1m, int L1n, int ldl1,
                         const CHAM_desc_t *L2, int L2m, int L2n, int ldl2,
                         const int *IPIV );
void INSERT_TASK_zsymm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zsyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zsyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int LDB,
                         CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zsytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An, int lda,
                               int iinfo );
void INSERT_TASK_zswpab( const RUNTIME_option_t *options,
                         int i, int n1, int n2,
                         const CHAM_desc_t *A, int Am, int An, int szeA );
void INSERT_TASK_zswptr_ontile( const RUNTIME_option_t *options,
                                const CHAM_desc_t descA, const CHAM_desc_t *Aij, int Aijm, int Aijn,
                                int i1,  int i2, int *ipiv, int inc,
                                const CHAM_desc_t *Akk, int Akkm, int Akkn, int ldak );
void INSERT_TASK_ztplqt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt );
void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                          const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                          const CHAM_desc_t *A, int Am, int An, int lda,
                          const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int m, int n, int k, int l, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                          const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                          const CHAM_desc_t *A, int Am, int An, int lda,
                          const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt );
void INSERT_TASK_ztrdalg( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int N, int NB,
                          const CHAM_desc_t *A,
                          const CHAM_desc_t *C, int Cm, int Cn,
                          const CHAM_desc_t *S, int Sm, int Sn,
                          int i, int j, int m, int grsiz, int BAND,
                          int *PCOL, int *ACOL, int *MCOL );
void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_ztrasm( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_ztrmm_p2( const RUNTIME_option_t *options,
                           cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                           int m, int n, int nb,
                           CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                           CHAMELEON_Complex64_t **B, int ldb );
void INSERT_TASK_ztrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,

                         int iinfo );
void INSERT_TASK_ztsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                                const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                                const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                                const CHAM_desc_t *T, int Tm, int Tn, int ldt );
void INSERT_TASK_ztsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                                const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                                const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                                const CHAM_desc_t *T, int Tm, int Tn, int ldt );
void INSERT_TASK_ztstrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un, int ldu,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *L, int Lm, int Ln, int ldl,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo );
void INSERT_TASK_zpamm( const RUNTIME_option_t *options,
                        int op, cham_side_t side, cham_store_t storev,
                        int m, int n, int k, int l,
                        const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                        const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                        const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                        const CHAM_desc_t *W, int Wm, int Wn, int ldw );
void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn );
void INSERT_TASK_zplssq2( const RUNTIME_option_t *options,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn );
void INSERT_TASK_zunmlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int ib,  int nb, int k,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                         const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zunmqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                         const CHAM_desc_t *C, int Cm, int Cn, int ldc );
void INSERT_TASK_zbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An, int lda,
                         void *user_data, void* user_build_callback );


/**
 * Keep these insert_task for retro-compatibility
 */
static inline void
INSERT_TASK_ztslqt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    return INSERT_TASK_ztplqt( options, m, n, 0, ib, nb,
                               A1, A1m, A1n, lda1,
                               A2, A2m, A2n, lda2,
                               T,  Tm,  Tn,  ldt );
}

static inline void
INSERT_TASK_ztsqrt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    return INSERT_TASK_ztpqrt( options, m, n, 0, ib, nb,
                               A1, A1m, A1n, lda1,
                               A2, A2m, A2n, lda2,
                               T,  Tm,  Tn,  ldt );
}

static inline void
INSERT_TASK_zttlqt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    return INSERT_TASK_ztplqt( options, m, n, n, ib, nb,
                               A1, A1m, A1n, lda1,
                               A2, A2m, A2n, lda2,
                               T,  Tm,  Tn,  ldt );
}

static inline void
INSERT_TASK_zttqrt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    return INSERT_TASK_ztpqrt( options, m, n, m, ib, nb,
                               A1, A1m, A1n, lda1,
                               A2, A2m, A2n, lda2,
                               T,  Tm,  Tn,  ldt );
}

static inline void
INSERT_TASK_ztsmlq( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    (void)m1;
    (void)n1;
    return INSERT_TASK_ztpmlqt( options, side, trans, m2, n2, k, 0, ib, nb,
                                V, Vm, Vn, ldv, T, Tm, Tn, ldt,
                                A1, A1m, A1n, lda1, A2, A2m, A2n, lda2 );
}

static inline void
INSERT_TASK_ztsmqr( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    (void)m1;
    (void)n1;
    return INSERT_TASK_ztpmqrt( options, side, trans, m2, n2, k, 0, ib, nb,
                                V, Vm, Vn, ldv, T, Tm, Tn, ldt,
                                A1, A1m, A1n, lda1, A2, A2m, A2n, lda2 );
}

static inline void
INSERT_TASK_zttmlq( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    (void)m1;
    (void)n1;
    return INSERT_TASK_ztpmlqt( options, side, trans, m2, n2, k, n2, ib, nb,
                                V, Vm, Vn, ldv, T, Tm, Tn, ldt,
                                A1, A1m, A1n, lda1, A2, A2m, A2n, lda2 );
}

static inline void
INSERT_TASK_zttmqr( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                    const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                    const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt )
{
    (void)m1;
    (void)n1;
    return INSERT_TASK_ztpmqrt( options, side, trans, m2, n2, k, m2, ib, nb,
                                V, Vm, Vn, ldv, T, Tm, Tn, ldt,
                                A1, A1m, A1n, lda1, A2, A2m, A2n, lda2 );
}

#endif /* _chameleon_tasks_z_h_ */
