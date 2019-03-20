/**
 *
 * @file cuda_ztpmqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztpmqrt GPU kernel
 *
 * @version 0.9.2
 * @author Florent Pruvost
 * @date 2016-12-22
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * @brief Applies a complex orthogonal matrix Q.
 *
 * The matrix Q is obtained from a "triangular-pentagonal" complex block
 * reflector H to a general complex matrix C, which consists of two blocks A and
 * B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg ChamLeft  : apply Q or Q**H from the Left;
 *         @arg ChamRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg ChamNoTrans   :  No transpose, apply Q;
 *         @arg ChamConjTrans :  ConjTranspose, apply Q**H.
 *
 * @param[in] M
 *         The number of rows of the tile B. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile B. N >= 0.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] L
 *          The number of rows of the upper trapezoidal part of V.
 *          K >= L >= 0.  See Further Details.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTPQRT in the first k rows of its array argument V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] A
 *         A is COMPLEX*16 array, dimension (LDA,N) if side = ChamLeft
 *         or (LDA,K) if SIDE = ChamRight
 *         On entry, the K-by-N or M-by-K matrix A.
 *         On exit, A is overwritten by the corresponding block of
 *         Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
 *
 * @param[in] LDA
 *         The leading dimension of the array A. LDA >= max(1,M).
 *         If side = ChamLeft,  LDA >= max(1,K);
 *         If side = Chamright, LDA >= max(1,M).
 *
 * @param[in,out] B
 *         On entry, the M-by-N tile B.
 *         On exit, B is overwritten by the corresponding block of
 *         Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
 *
 * @param[in] LDB
 *         The leading dimension of the tile B. LDB >= max(1,M).
 *
 * @param[out] WORK
 *         Workspace array of size LDWORK-by-NB.
 *         LDWORK = N if side = ChamLeft, or  M if side = ChamRight.
 *
 *******************************************************************************
 *
 * @par Further Details:
 * =====================
 *
 *  The columns of the pentagonal matrix V contain the elementary reflectors
 *  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
 *  trapezoidal block V2:
 *
 *        V = [V1]
 *            [V2].
 *
 *  The size of the trapezoidal block V2 is determined by the parameter L,
 *  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L
 *  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular;
 *  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
 *
 *  If side = ChamLeft:  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K.
 *                            [B]
 *
 *  If side = ChamRight: C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K.
 *
 *  The complex orthogonal matrix Q is formed from V and T.
 *
 *  If trans='N' and side='L', C is on exit replaced with Q * C.
 *
 *  If trans='C' and side='L', C is on exit replaced with Q**H * C.
 *
 *  If trans='N' and side='R', C is on exit replaced with C * Q.
 *
 *  If trans='C' and side='R', C is on exit replaced with C * Q**H.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CUDA_ztpmqrt( cham_side_t side, cham_trans_t trans,
              int M, int N, int K, int L, int IB,
              const cuDoubleComplex *V, int LDV,
              const cuDoubleComplex *T, int LDT,
                    cuDoubleComplex *A, int LDA,
                    cuDoubleComplex *B, int LDB,
                    cuDoubleComplex *WORK, int lwork,
              CUBLAS_STREAM_PARAM )
{
    int m1, n1;

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        cudablas_error(1, "Illegal value of side");
        return -1;
    }

    if ( side == ChamLeft ) {
        m1 = K;
        n1 = N;
    }
    else {
        m1 = M;
        n1 = K;
    }

    /* TS case */
    if (L == 0) {
        CUDA_ztsmqr( side, trans, m1, n1, M, N, K, IB,
                     A, LDA, B, LDB, V, LDV, T, LDT,
                     WORK, lwork,
                     CUBLAS_STREAM_VALUE );
    }
    /* TT case */
    else  if( L == M ) {
        CUDA_zttmqr( side, trans, m1, n1, M, N, K, IB,
                     A, LDA, B, LDB, V, LDV, T, LDT,
                     WORK, lwork,
                     CUBLAS_STREAM_VALUE );
    }
    else {
        cudablas_error(-6, "TPMQRT not available on GPU for general cases yet\n" );
        return -6;
    }

    return CHAMELEON_SUCCESS;
}
