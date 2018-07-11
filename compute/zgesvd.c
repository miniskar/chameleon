/**
 *
 * @file zgesvd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesvd wrappers
 *
 * @version 1.0.0
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/lapacke.h>
#endif

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgesvd - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors. The SVD is written
 *
 *       A = U * SIGMA * transpose(V)
 *
 *  where SIGMA is an M-by-N matrix which is zero except for its
 *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 *  are the singular values of A; they are real and non-negative, and
 *  are returned in descending order.  The first min(m,n) columns of
 *  U and V are the left and right singular vectors of A.
 *
 *  Note that the routine returns V**T, not V.
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = ChamVec   = 'A'(lapack):  all M columns of U are returned
 *                        in array U;
 *          = ChamNoVec = 'N':  no columns of U (no left singular vectors)
 *                        are computed.
 *          = ChamSVec  = 'S': the first min(m,n) columns of U (the left
 *                        singular vectors) are returned in the array U;
 *                        NOT SUPPORTTED YET
 *          = ChamOVec  = 'O': the first min(m,n) columns of U (the left
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**H.
 *          Intended usage:
 *          = ChamVec   = 'A'(lapack): all N rows of V**H are returned
 *                        in the array VT;
 *          = ChamNoVec = 'N': no rows of V**H (no right singular vectors)
 *                        are computed.
 *          = ChamSVec  = 'S': the first min(m,n) rows of V**H (the right
 *                        singular vectors) are returned in the array VT;
 *                        NOT SUPPORTTED YET
 *          = ChamOVec  = 'O': the first min(m,n) rows of V**H (the right
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 *          Note: jobu and jobvt cannot both be ChamOVec.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**H (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] S
 *          The double precision singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by CHAMELEON_Alloc_Workspace_zgesvd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1; if
 *          JOBU = 'S' or 'A', LDU >= M.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**H;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**H (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[in] LDVT
 *         The leading dimension of the array VT.  LDVT >= 1; if
 *         JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgesvd_Tile
 * @sa CHAMELEON_zgesvd_Tile_Async
 * @sa CHAMELEON_cgesvd
 * @sa CHAMELEON_dgesvd
 * @sa CHAMELEON_sgesvd
 *
 */
int CHAMELEON_zgesvd( cham_job_t jobu, cham_job_t jobvt,
                  int M, int N,
                  CHAMELEON_Complex64_t *A, int LDA,
                  double *S,
                  CHAM_desc_t *descT,
                  CHAMELEON_Complex64_t *U, int LDU,
                  CHAMELEON_Complex64_t *VT, int LDVT )
{
    int NB;
    int status;
    CHAM_context_t  *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zgesvd", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobu != ChamNoVec  && jobu != ChamVec) {
        morse_error("CHAMELEON_zgesvd", "illegal value of jobu");
        return -1;
    }
    if (jobvt != ChamNoVec && jobvt != ChamVec) {
        morse_error("CHAMELEON_zgesvd", "illegal value of jobvt");
        return -2;
    }
    if (M < 0) {
        morse_error("CHAMELEON_zgesvd", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        morse_error("CHAMELEON_zgesvd", "illegal value of N");
        return -4;
    }
    if (LDA < chameleon_max(1, M)) {
        morse_error("CHAMELEON_zgesvd", "illegal value of LDA");
        return -6;
    }
    if (LDU < 1) {
        morse_error("CHAMELEON_zgesvd", "illegal value of LDU");
        return -9;
    }
    if (LDVT < 1) {
        morse_error("CHAMELEON_zgesvd", "illegal value of LDVT");
        return -11;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0) {
        return CHAMELEON_SUCCESS;
    }

    /* Tune NB & IB depending on M & N; Set NBNB */
    status = morse_tune(CHAMELEON_FUNC_ZGESVD, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zgesvd", "morse_tune() failed");
        return status;
    }

    /* Set MT, NT */
    NB = CHAMELEON_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, ChamDescInout, ChamUpperLower,
                     A, NB, NB,  LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgesvd_Tile_Async( jobu, jobvt, &descAt, S, descT, U, LDU, VT, LDVT, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgesvd_Tile - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Tile equivalent of CHAMELEON_zgesvd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = ChamVec   = 'A'(lapack):  all M columns of U are returned
 *                        in array U;
 *          = ChamNoVec = 'N':  no columns of U (no left singular vectors)
 *                        are computed.
 *          = ChamSVec  = 'S': the first min(m,n) columns of U (the left
 *                        singular vectors) are returned in the array U;
 *                        NOT SUPPORTTED YET
 *          = ChamOVec  = 'O': the first min(m,n) columns of U (the left
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**H.
 *          Intended usage:
 *          = ChamVec   = 'A'(lapack): all N rows of V**H are returned
 *                        in the array VT;
 *          = ChamNoVec = 'N': no rows of V**H (no right singular vectors)
 *                        are computed.
 *          = ChamSVec  = 'S': the first min(m,n) rows of V**H (the right
 *                        singular vectors) are returned in the array VT;
 *                        NOT SUPPORTTED YET
 *          = ChamOVec  = 'O': the first min(m,n) rows of V**H (the right
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 *          Note: jobu and jobvt cannot both be ChamOVec.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**H (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[out] S
 *          The singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[in, out] T
 *          On entry, descriptor as return by CHAMELEON_Alloc_Workspace_zgesvd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1; if
 *          JOBU = 'S' or 'A', LDU >= M.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**H;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**H (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[in] LDVT
 *         The leading dimension of the array VT.  LDVT >= 1; if
 *         JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
 *
 *******************************************************************************
 *
 * @return
 *          \return CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgesvd
 * @sa CHAMELEON_zgesvd_Tile_Async
 * @sa CHAMELEON_cgesvd_Tile
 * @sa CHAMELEON_dgesvd_Tile
 * @sa CHAMELEON_sgesvd_Tile
 *
 */
int CHAMELEON_zgesvd_Tile( cham_job_t jobu, cham_job_t jobvt,
                       CHAM_desc_t *A,
                       double *S,
                       CHAM_desc_t *T,
                       CHAMELEON_Complex64_t *U, int LDU,
                       CHAMELEON_Complex64_t *VT, int LDVT )
{
    CHAM_context_t *morse;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zgesvd_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    CHAMELEON_zgesvd_Tile_Async( jobu, jobvt, A, S, T, U, LDU, VT, LDVT, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( T, sequence );

    morse_sequence_wait( morse, sequence );
    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zgesvd_Tile_Async - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Non-blocking equivalent of CHAMELEON_zgesvd_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgesvd
 * @sa CHAMELEON_zgesvd_Tile
 * @sa CHAMELEON_cgesvd_Tile_Async
 * @sa CHAMELEON_dgesvd_Tile_Async
 * @sa CHAMELEON_sgesvd_Tile_Async
 *
 */
int CHAMELEON_zgesvd_Tile_Async( cham_job_t jobu, cham_job_t jobvt,
                             CHAM_desc_t *A,
                             double *S,
                             CHAM_desc_t *T,
                             CHAMELEON_Complex64_t *U, int LDU,
                             CHAMELEON_Complex64_t *VT, int LDVT,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_desc_t descA;
    CHAM_desc_t descT;
    CHAM_desc_t descUl, descUt;
    CHAM_desc_t descVTl, descVTt;
    CHAM_desc_t descAB;
    CHAM_desc_t D, *Dptr = NULL;
    CHAM_desc_t *subA, *subT, *subUVT;
    double *E;
    int M, N, MINMN, NB, LDAB;
    int KL, KU, uplo, nru, ncvt;
    int info = 0;
    char gbbrd_vect;

    CHAM_context_t *morse;
    morse = morse_context_self();

    if (morse == NULL) {
        morse_fatal_error("CHAMELEON_zgesvd_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("CHAMELEON_zgesvd_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("CHAMELEON_zgesvd_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return morse_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zgesvd_Tile_Async", "invalid first descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (morse_desc_check(T) != CHAMELEON_SUCCESS) {
        morse_error("CHAMELEON_zgesvd_Tile_Async", "invalid fourth descriptor");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    /* Check input arguments */
    if (jobu != ChamNoVec  && jobu != ChamVec) {
        morse_error("CHAMELEON_zgesvd_Tile_Async", "illegal value of jobu");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    if (jobvt != ChamNoVec && jobvt != ChamVec) {
        morse_error("CHAMELEON_zgesvd_Tile_Async", "illegal value of jobvt");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    if (descA.nb != descA.mb) {
        morse_error("CHAMELEON_zgesvd_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    M     = descA.m;
    N     = descA.n;
    MINMN = chameleon_min(M, N);
    NB    = descA.mb;
    LDAB  = NB + 1;
    uplo  = M >= N ? ChamUpper : ChamLower;

#if defined(CHAMELEON_COPY_DIAG)
    {
        morse_zdesc_alloc(D, A->mb, A->nb, A->m, A->n, 0, 0, A->m, A->n, );
        Dptr = &D;
    }
#endif
    /* Reduction to band */
    morse_pzgebrd_ge2gb( &descA, &descT, Dptr,
                         sequence, request );

    /* Allocate band structure */
    morse_zdesc_alloc_diag( descAB,
                            LDAB, NB,
                            LDAB, MINMN,
                            0, 0,
                            LDAB, MINMN,
                            1, 1 );

    /* Convert matrix to band form */
    morse_pztile2band( uplo,
                       &descA, &descAB,
                       sequence, request );

    E = malloc( MINMN * sizeof(double) );
    if (E == NULL) {
        morse_error("CHAMELEON_zheevd_Tile_Async", "malloc(E) failed");
        free(E);
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }
    memset(E, 0, MINMN * sizeof(double) );

#if !defined(CHAMELEON_SIMULATION)
    /* NCC = 0, C = NULL, we do not update any matrix with new singular vectors */
    /* On exit, AB = U (S +~ E) VT */
    if (uplo == ChamUpper){
        KL = 0;
        KU = NB;
    }
    else{
        KL = NB;
        KU = 0;
    }

    /* Manage the case where only singular values are required */
    if (jobu == ChamNoVec) {
        nru = 0;
        if (jobvt == ChamNoVec) {
            gbbrd_vect = 'N';
            ncvt = 0;
        }
        else {
            gbbrd_vect = 'P';
            ncvt = N;
        }
    }
    else {
        nru = M;
        if (jobvt == ChamNoVec) {
            gbbrd_vect = 'Q';
            ncvt = 0;
        }
        else {
            gbbrd_vect = 'B';
            ncvt = N;
        }
    }

    morse_sequence_wait( morse, sequence );

    info = LAPACKE_zgbbrd( LAPACK_COL_MAJOR,
                           gbbrd_vect,
                           M, N,
                           0, KL, KU,
                           (CHAMELEON_Complex64_t *) descAB.mat, LDAB,
                           S, E,
                           U, LDU,
                           VT, LDVT,
                           NULL, 1 );
    if (info != 0) {
        fprintf(stderr, "CHAMELEON_zgesvd_Tile_Async: LAPACKE_zgbbrd = %d\n", info );
    }
#else
    morse_sequence_wait( morse, sequence );
#endif /* !defined(CHAMELEON_SIMULATION) */

    morse_desc_mat_free( &descAB );

    subA = NULL;
    subT = NULL;
    subUVT = NULL;

    if ( jobu != ChamNoVec ) {
        morse_zlap2tile( morse, &descUl, &descUt, ChamDescInout, ChamUpperLower,
                         U, NB, NB, LDU, M, M, M, sequence, request );

        if ( M < N ){
            subA   = morse_desc_submatrix(&descA,  descA.mb,  0, descA.m -descA.mb,  descA.n-descA.nb);
            subUVT = morse_desc_submatrix(&descUt, descUt.mb, 0, descUt.m-descUt.mb, descUt.n);
            subT   = morse_desc_submatrix(&descT,  descT.mb,  0, descT.m -descT.mb,  descT.n-descT.nb);

            morse_pzunmqr( ChamLeft, ChamNoTrans,
                           subA, subUVT, subT, Dptr,
                           sequence, request );

            free(subA); free(subUVT); free(subT);
        }
        else {
            morse_pzunmqr( ChamLeft, ChamNoTrans,
                           &descA, &descUt, &descT, Dptr,
                           sequence, request );
        }

        morse_ztile2lap( morse, &descUl, &descUt,
                         ChamDescInout, ChamUpperLower, sequence, request );
    }

    if ( jobvt != ChamNoVec ) {
        morse_zlap2tile( morse, &descVTl, &descVTt, ChamDescInout, ChamUpperLower,
                         VT, NB, NB, LDVT, N, N, N, sequence, request );

        if ( M < N ){
            morse_pzunmlq( ChamRight, ChamNoTrans,
                           &descA, &descVTt, &descT, Dptr,
                           sequence, request );
        }
        else {
            subA   = morse_desc_submatrix(&descA,   0, descA.nb,   descA.m-descA.mb, descA.n  -descA.nb  );
            subUVT = morse_desc_submatrix(&descVTt, 0, descVTt.nb, descVTt.m,        descVTt.n-descVTt.nb);
            subT   = morse_desc_submatrix(&descT,   0, descT.nb,   descT.m-descT.mb, descT.n  -descT.nb  );

            morse_pzunmlq( ChamRight, ChamNoTrans,
                           subA, subUVT, subT, Dptr,
                           sequence, request );

            free(subA); free(subUVT); free(subT);
        }

        morse_ztile2lap( morse, &descVTl, &descVTt,
                         ChamDescInout, ChamUpperLower, sequence, request );
    }
    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    if ( jobu != ChamNoVec ) {
        morse_ztile2lap_cleanup( morse, &descUl,  &descUt  );
    }
    if ( jobvt != ChamNoVec ) {
        morse_ztile2lap_cleanup( morse, &descVTl, &descVTt );
    }

    /* Solve the bidiagonal SVD problem */
    /* On exit, U and VT are updated with bidiagonal matrix singular vectors */
#if !defined(CHAMELEON_SIMULATION)
    info = LAPACKE_zbdsqr( LAPACK_COL_MAJOR, 'U',
                           MINMN, ncvt, nru, 0,
                           S, E,
                           VT, LDVT, U, LDU, NULL, 1 );
    if (info != 0) {
        fprintf(stderr, "CHAMELEON_zgesvd_Tile_Async: LAPACKE_zbdsqr = %d\n", info );
    }
#endif /* !defined(CHAMELEON_SIMULATION) */

    free(E);
    if ( Dptr ) {
        morse_desc_mat_free( Dptr );
    }
    (void)D;
    return CHAMELEON_SUCCESS;
}
