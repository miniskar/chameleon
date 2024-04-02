/**
 *
 * @file compute_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational functions header
 *
 * @version 1.3.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Matthieu Kuhn
 * @author Lionel Eyraud-Dubois
 * @date 2023-09-08
 * @precisions normal z -> c d s
 *
 */
#ifndef _compute_z_h_
#define _compute_z_h_

/**
 * @brief Data structure to handle the GEMM workspaces
 */
struct chameleon_pzgemm_s {
    cham_gemm_t alg;
    CHAM_desc_t WA;
    CHAM_desc_t WB;
};

/**
 * @brief Data structure to handle the GETRF workspaces with partial pivoting
 */
struct chameleon_pzgetrf_s {
    cham_getrf_t alg;
    int          ib;         /**< Internal blocking parameter */
    int          batch_size; /**< Batch size for the panel    */
    CHAM_desc_t  U;
    CHAM_desc_t  Up;
};

/**
 * @brief Data structure to handle the Centering-Scaled workspaces
 */
struct chameleon_pzcesca_s {
    CHAM_desc_t Wgcol;
    CHAM_desc_t Wgrow;
    CHAM_desc_t Wgelt;
    CHAM_desc_t Wdcol;
    CHAM_desc_t Wdrow;
};

/**
 * @brief Data structure to handle the GRAM workspaces
 */
struct chameleon_pzgram_s {
    CHAM_desc_t Wcol;
    CHAM_desc_t Welt;
};

/**
 *  Declarations of internal sequential functions
 */
int chameleon_zshift(CHAM_context_t *chamctxt, int m, int n, CHAMELEON_Complex64_t *A,
                     int nprob, int me, int ne, int L,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

/**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 */
#if defined(PRECISION_z) || defined(PRECISION_d)
void chameleon_pzgered( cham_uplo_t uplo, double prec, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgerst( cham_uplo_t uplo, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
#endif
int chameleon_pzgebrd( int genD, cham_job_t jobu, cham_job_t jobvt,
                       CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                       CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *VT, int LDVT,
                       double *E, double *S, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgelqf( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgelqfrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgenm2( double tol, const CHAM_desc_t *A, double *result, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgemm( struct chameleon_pzgemm_s *options, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgepdf_qdwh( cham_mtxtype_t trans, CHAM_desc_t *descU, CHAM_desc_t *descH, gepdf_info_t *info, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgepdf_qr( int genD, int doqr, int optid, const libhqr_tree_t *qrtreeT, const libhqr_tree_t *qrtreeB, CHAM_desc_t *A1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *D1, CHAM_desc_t *Q1, CHAM_desc_t *A2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *D2, CHAM_desc_t *Q2, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgeqrf( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgeqrfrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf( struct chameleon_pzgetrf_s *ws, CHAM_desc_t *A, CHAM_ipiv_t *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzgetrf_incpiv(CHAM_desc_t *A, CHAM_desc_t *L, CHAM_desc_t *D, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf_nopiv(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf_reclap(CHAM_desc_t *A, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgetrf_rectil(CHAM_desc_t *A, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzhegst(int itype, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzhemm( struct chameleon_pzgemm_s *ws,cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzherk(cham_uplo_t uplo, cham_trans_t trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzher2k(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzhetrd_he2hb(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *E, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlacpy(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlag2c(CHAM_desc_t *A, CHAM_desc_t *SB, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlange_generic( cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                                double *result, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlansy_generic( cham_normtype_t norm, cham_uplo_t uplo, cham_trans_t trans,
                                CHAM_desc_t *A, double *result,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlascal(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlaset( cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlaset2(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha,                          CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlaswp(CHAM_desc_t *B, int *IPIV, int inc, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlaswpc(CHAM_desc_t *B, int *IPIV, int inc, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzlatms( cham_dist_t idist, unsigned long long int seed, cham_sym_t sym, double *D, int mode, double cond, double dmax, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzlauum(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzplghe(double bump, cham_uplo_t uplo, CHAM_desc_t *A, int bigM, int m0, int n0, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzplgsy(CHAMELEON_Complex64_t bump, cham_uplo_t uplo, CHAM_desc_t *A, int bigM, int m0, int n0, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzplrnt(CHAM_desc_t *A, int bigM, int m0, int n0, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzplrnk(int K, CHAM_desc_t *C, unsigned long long int seedA, unsigned long long int seedB, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzpotrf(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzpotrimm(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzshift(int, int, int, CHAMELEON_Complex64_t *, int *, int, int, int, RUNTIME_sequence_t*, RUNTIME_request_t*);
void chameleon_pzsymm( struct chameleon_pzgemm_s *ws,cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pzsyrk(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta,  CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzsyr2k(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzsytrf(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztile2band(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *descAB, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztpgqrt( int kt, int L, CHAM_desc_t *V2, CHAM_desc_t *T2, CHAM_desc_t *Q1, CHAM_desc_t *Q2, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pztpqrt( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pztradd(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrmm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrsmpl(CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *L, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrsmrv(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *W, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztrtri(cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungbr(cham_side_t side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungbrrh(cham_side_t side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungqr( int genD, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungqrrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunglq( int genD, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunglqrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungtr(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmqr( int genD, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmqrrh( int genD, int BS, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmlq( int genD, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmlqrh( int genD, int BS, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzbuild( cham_uplo_t uplo, CHAM_desc_t *A, void *user_data, void* user_build_callback, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

int chameleon_pzgelqf_step( int genD, int k, int ib,
                            CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                            RUNTIME_option_t *options, RUNTIME_sequence_t *sequence );
int chameleon_pzgeqrf_step( int genD, int k, int ib,
                            CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                            RUNTIME_option_t *options, RUNTIME_sequence_t *sequence );
int  chameleon_pzgelqf_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                   const libhqr_tree_t *qrtree, int *tiles,
                                   CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                   RUNTIME_option_t *options, RUNTIME_sequence_t *sequence );
int  chameleon_pzgeqrf_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                   const libhqr_tree_t *qrtree, int *tiles,
                                   CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                   RUNTIME_option_t *options, RUNTIME_sequence_t *sequence );
void chameleon_pzungqr_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                   const libhqr_tree_t *qrtree, int nbtiles, int *tiles,
                                   CHAM_desc_t *A, CHAM_desc_t *Q,
                                   CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                   RUNTIME_option_t *options, RUNTIME_sequence_t *sequence );
void chameleon_pzgelqf_param( int genD, int K, const libhqr_tree_t *qrtree,
                              CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzgeqrf_param( int genD, int K, const libhqr_tree_t *qrtree,
                              CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmlq_param( int genD, const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunmqr_param( int genD, const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzunglq_param( int genD, const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *Q,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pzungqr_param( int genD, int K, const libhqr_tree_t *qrtree,
                              CHAM_desc_t *A, CHAM_desc_t *Q,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void chameleon_pztpgqrt_param( int genD, cham_uplo_t uplo, int kt, const libhqr_tree_t *qrtree,
                               CHAM_desc_t *V2, CHAM_desc_t *Q1, CHAM_desc_t *Q2,
                               CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *DD,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void chameleon_pztpqrt_param( int genD, cham_uplo_t uplo, int K,
                              const libhqr_tree_t *qrtree,
                              CHAM_desc_t *ATop, CHAM_desc_t *A,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

/**
 * Centered-Scaled function prototypes
 */
void chameleon_pzcesca( struct chameleon_pzcesca_s *ws, int center, int scale, cham_store_t axis, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
/**
 * Gram function prototypes
 */
void chameleon_pzgram( struct chameleon_pzgram_s *ws, cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

/**
 *  LAPACK/Tile Descriptor accesses
 */
#define ChamDescInput  1
#define ChamDescOutput 2
#define ChamDescInout  (ChamDescInput | ChamDescOutput)

/**
 *  Macro for matrix conversion / Lapack interface
 */
static inline int
chameleon_zdesc_alloc_diag( CHAM_desc_t *descA, int nb, int m, int n, int p, int q ) {
    int diag_m = chameleon_min( m, n );
    return chameleon_desc_init( descA, CHAMELEON_MAT_ALLOC_TILE,
                                ChamComplexDouble, nb, nb, nb*nb,
                                diag_m, nb, 0, 0, diag_m, nb, p, q,
                                chameleon_getaddr_diag,
                                chameleon_getblkldd_ccrb,
                                chameleon_getrankof_2d_diag, NULL );
}

#define chameleon_zdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free) \
    {                                                                   \
        int rc;                                                         \
        rc = chameleon_desc_init( &(descA), CHAMELEON_MAT_ALLOC_GLOBAL, \
                                  ChamComplexDouble, (mb), (nb), ((mb)*(nb)), \
                                  (m), (n), (i), (j), (m), (n), 1, 1,   \
                                  NULL, NULL, NULL, NULL );             \
        if ( rc != CHAMELEON_SUCCESS ) {                                \
            {free;}                                                     \
            return rc;                                                  \
        }                                                               \
    }

/**
 * @brief Create a copy of a descriptor restricted to a smaller size.
 * @param[in]  descIn  The input descriptor from which the structure should be copied.
 * @param[out] descOut The output descriptor that is a copy of the input one with allocation on the fly.
 * @param[in]  m       The number of rows of the output descriptor.
 * @param[in]  n       The number of columns of the output descriptor.
 * @return CHAMELEON_SUCCESS on success, the associated error on failure.
 */
static inline int
chameleon_zdesc_copy_and_restrict( const CHAM_desc_t *descIn,
                                   CHAM_desc_t *descOut,
                                   int m, int n )
{
    int rc;
    rc = chameleon_desc_init( descOut, CHAMELEON_MAT_ALLOC_TILE,
                              ChamComplexDouble, descIn->mb, descIn->nb, descIn->mb * descIn->nb,
                              m, n, 0, 0, m, n, descIn->p, descIn->q,
                              descIn->get_blkaddr,
                              descIn->get_blkldd,
                              descIn->get_rankof_init, descIn->get_rankof_init_arg );
    return rc;
}

/**
 * @brief Internal function to convert the lapack format to tile format in
 * LAPACK interface calls
 */
static inline int
chameleon_zlap2tile( CHAM_context_t *chamctxt,
                     CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                     int mode, cham_uplo_t uplo,
                     CHAMELEON_Complex64_t *A, int mb, int nb, int lm, int ln, int m, int n,
                     RUNTIME_sequence_t *seq, RUNTIME_request_t *req )
{
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        /* Initialize the Lapack descriptor */
        chameleon_desc_init( descAl, A, ChamComplexDouble, mb, nb, (mb)*(nb),
                            lm, ln, 0, 0, m, n, 1, 1,
                             chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );
        descAl->styp = ChamCM;

        /* Initialize the tile descriptor */
        chameleon_desc_init( descAt, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, mb, nb, (mb)*(nb),
                             lm, ln, 0, 0, m, n, 1, 1,
                             chameleon_getaddr_ccrb, chameleon_getblkldd_ccrb, NULL, NULL );

        if ( mode & ChamDescInput ) {
            chameleon_pzlacpy( uplo, descAl, descAt, seq, req );
        }
    }
    else {
        /* Initialize the tile descriptor */
        chameleon_desc_init( descAt, A, ChamComplexDouble, mb, nb, (mb)*(nb),
                             lm, ln, 0, 0, m, n, 1, 1,
                             chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );
    }
    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to convert back the tile format to the lapack format
 * in LAPACK interface calls
 */
static inline int
chameleon_ztile2lap( CHAM_context_t *chamctxt, CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                     int mode, cham_uplo_t uplo, RUNTIME_sequence_t *seq, RUNTIME_request_t *req )
{
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        if ( mode & ChamDescOutput ) {
            chameleon_pzlacpy( uplo, descAt, descAl, seq, req );
        }
        RUNTIME_desc_flush( descAl, seq );
    }
    RUNTIME_desc_flush( descAt, seq );

    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to cleanup the temporary data from the layout
 * conversions in LAPACK interface calls
 */
static inline void
chameleon_ztile2lap_cleanup( CHAM_context_t *chamctxt, CHAM_desc_t *descAl, CHAM_desc_t *descAt )
{
    (void)chamctxt;
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        chameleon_desc_destroy( descAl );
    }
    chameleon_desc_destroy( descAt );
}

#endif /* _compute_z_h_ */
