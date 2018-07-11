/**
 *
 * @file compute_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational functions header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
/**
 *  LAPACK/Tile Descriptor accesses
 */
#define ChamDescInput  1
#define ChamDescOutput 2
#define ChamDescInout  (ChamDescInput | ChamDescOutput)

/**
 *  Macro for matrix conversion / Lapack interface
 */
#define morse_zdesc_alloc_diag( descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = morse_desc_init_diag(                                       \
        ChamComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    morse_desc_mat_alloc( &(descA) );                                   \
    RUNTIME_desc_create( &(descA) );

#define morse_zdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = morse_desc_init(                                            \
        ChamComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                           \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return CHAMELEON_ERR_OUT_OF_RESOURCES;                              \
    }                                                                   \
    RUNTIME_desc_create( &(descA) );

/**
 *  Declarations of internal sequential functions
 */
int morse_zshift(CHAM_context_t *morse, int m, int n, CHAMELEON_Complex64_t *A,
                  int nprob, int me, int ne, int L,
                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

/**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 */
void morse_pzbarrier_pnl2tl(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzbarrier_row2tl(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzbarrier_tl2pnl(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzbarrier_tl2row(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgebrd_gb2bd(cham_uplo_t uplo, CHAM_desc_t *A, double *D, double *E, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgebrd_ge2gb(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgelqf(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgelqfrh(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgemm(cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgeqrf(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgeqrfrh(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgetrf_incpiv(CHAM_desc_t *A, CHAM_desc_t *L, CHAM_desc_t *D, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgetrf_nopiv(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgetrf_reclap(CHAM_desc_t *A, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgetrf_rectil(CHAM_desc_t *A, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzhegst(int itype, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzhemm(cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzherk(cham_uplo_t uplo, cham_trans_t trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzher2k(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzhetrd_he2hb(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *E, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlacpy(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlag2c(CHAM_desc_t *A, CHAM_desc_t *SB, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlange(cham_normtype_t norm, CHAM_desc_t *A, double *result, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlanhe(cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A, double *result, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlansy(cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A, double *result, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlantr(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, double *result, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlascal(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlaset( cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlaset2(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha,                          CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlaswp(CHAM_desc_t *B, int *IPIV, int inc, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlaswpc(CHAM_desc_t *B, int *IPIV, int inc, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzlauum(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzplghe(double bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void morse_pzplgsy(CHAMELEON_Complex64_t bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void morse_pzplrnt(CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void morse_pzpotrf(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzpotrimm(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzshift(int, int, int, CHAMELEON_Complex64_t *, int *, int, int, int, RUNTIME_sequence_t*, RUNTIME_request_t*);
void morse_pzsymm(cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzsyrk(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta,  CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzsyr2k(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzsytrf(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pztile2band(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *descAB, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pztpgqrt( int L, CHAM_desc_t *V1, CHAM_desc_t *T1, CHAM_desc_t *V2, CHAM_desc_t *T2, CHAM_desc_t *Q1, CHAM_desc_t *Q2, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void morse_pztpqrt( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
void morse_pztradd(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pztrmm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pztrsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pztrsmpl(CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *L, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pztrsmrv(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *W, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pztrtri(cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzungbr(cham_side_t side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzungbrrh(cham_side_t side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzungqr(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzungqrrh(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D,int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunglq(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunglqrh(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzungtr(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunmqr(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunmqrrh(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunmlq(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunmlqrh(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D, int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzbuild( cham_uplo_t uplo, CHAM_desc_t *A, void *user_data, void* user_build_callback, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

void morse_pzgelqf_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzgeqrf_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunmlq_param(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans,
                         CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunmqr_param(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans,
                         CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzunglq_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *Q,
                         CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
void morse_pzungqr_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *Q,
                         CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);


/**
 * @brief Internal function to convert the lapack format to tile format in
 * LAPACK interface calls
 */
static inline int
morse_zlap2tile( CHAM_context_t *morse,
                 CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                 int mode, cham_uplo_t uplo,
                 CHAMELEON_Complex64_t *A, int mb, int nb, int lm, int ln, int m, int n,
                 RUNTIME_sequence_t *seq, RUNTIME_request_t *req )
{
    /* Initialize the Lapack descriptor */
    *descAl = morse_desc_init_user( ChamComplexDouble, mb, nb, (mb)*(nb),
                                    lm, ln, 0, 0, m, n, 1, 1,
                                    morse_getaddr_cm, morse_getblkldd_cm, NULL  );
    descAl->mat = A;
    descAl->styp = ChamCM;

    /* Initialize the tile descriptor */
    *descAt = morse_desc_init( ChamComplexDouble, mb, nb, (mb)*(nb),
                               lm, ln, 0, 0, m, n, 1, 1 );

    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        if ( morse_desc_mat_alloc( descAt ) ) {
            morse_error( "morse_zlap2tile", "morse_desc_mat_alloc() failed");
            return CHAMELEON_ERR_OUT_OF_RESOURCES;
        }

        RUNTIME_desc_create( descAl );
        RUNTIME_desc_create( descAt );

        if ( mode & ChamDescInput ) {
            morse_pzlacpy( uplo, descAl, descAt, seq, req );
        }
    }
    else {
        morse_fatal_error( "morse_zlap2tile", "INPLACE translation not supported yet");
        descAt->mat = A;

        RUNTIME_desc_create( descAl );
        RUNTIME_desc_create( descAt );

        if ( mode & ChamDescInput ) {
            /* CHAMELEON_zgecfi_Async( lm, ln, A, ChamCM, mb, nb, */
            /*                     ChamCCRB, mb, nb, seq, req ); */
        }
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }

    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to convert back the tile format to the lapack format
 * in LAPACK interface calls
 */
static inline int
morse_ztile2lap( CHAM_context_t *morse, CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                 int mode, cham_uplo_t uplo, RUNTIME_sequence_t *seq, RUNTIME_request_t *req )
{
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        if ( mode & ChamDescOutput ) {
            morse_pzlacpy( uplo, descAt, descAl, seq, req );
        }
    }
    else {
        morse_fatal_error( "morse_ztile2lap", "INPLACE translation not supported yet");
        if ( mode & ChamDescOutput ) {
            /* CHAMELEON_zgecfi_Async( descAl->lm, descAl->ln, descAl->mat, */
            /*                     ChamCCRB, descAl->mb, descAl->nb,   */
            /*                     ChamCM, descAl->mb, descAl->nb, seq, req ); */
        }
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    RUNTIME_desc_flush( descAl, seq );
    RUNTIME_desc_flush( descAt, seq );

    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to cleanup the temporary data from the layout
 * conversions in LAPACK interface calls
 */
static inline void
morse_ztile2lap_cleanup( CHAM_context_t *morse, CHAM_desc_t *descAl, CHAM_desc_t *descAt )
{
    if ( CHAMELEON_TRANSLATION == ChamOutOfPlace ) {
        morse_desc_mat_free( descAt );
    }
    RUNTIME_desc_destroy( descAl );
    RUNTIME_desc_destroy( descAt );
}
