/**
 *
 * @file testing_zgesvd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesvd testing
 *
 * @version 1.0.0
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <chameleon.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <coreblas.h>
#include "testing_zauxiliary.h"

static int check_orthogonality( cham_side_t, int, int, CHAMELEON_Complex64_t*, int, double);
static int check_reduction(int, int, double*, CHAMELEON_Complex64_t*, int, CHAMELEON_Complex64_t*, int, CHAMELEON_Complex64_t*, int, double);
static int check_solution(int, double*, double*, double);

int testing_zgesvd(int argc, char **argv)
{
    int tree = 0;

    if ( argc < 1 ){
        goto usage;
    } else {
        tree = atoi(argv[0]);
    }

    /* Check for number of arguments*/
    if ( ((tree == 0) && (argc != 4)) ||
         ((tree != 0) && (argc != 5)) ){
      usage:
        USAGE("GESVD", "MODE M N LDA [RH]",
              "   - MODE : 0: flat, 1: tree (RH needed)\n"
              "   - M    : number of rows of the matrix A\n"
              "   - N    : number of columns of the matrix A\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - RH   : Size of each subdomains\n");
        return -1;
    }

    int M   = atoi(argv[1]);
    int N   = atoi(argv[2]);
    int LDA = atoi(argv[3]);
    int rh;
    if ( tree ) {
        rh = atoi(argv[4]);
        CHAMELEON_Set(CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder);
        CHAMELEON_Set(CHAMELEON_HOUSEHOLDER_SIZE, rh);
    }

    if (LDA < M){
        printf("LDA should be >= M !\n");
        return -1;
    }

    double eps  = LAPACKE_dlamch_work('e');
    double dmax = 1.0;
    cham_job_t jobu  = ChamVec;
    cham_job_t jobvt = ChamVec;
    int info_orthou    = 0;
    int info_orthovt   = 0;
    int info_solution  = 0;
    int info_reduction = 0;
    int minMN = min(M, N);
    int mode  = 4;
    double rcond = (double) minMN;
    int INFO=-1;

    CHAMELEON_Complex64_t *A1   = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
    double *S1              = (double *)           malloc(minMN*sizeof(double));
    double *S2              = (double *)           malloc(minMN*sizeof(double));
    CHAMELEON_Complex64_t *work = (CHAMELEON_Complex64_t *)malloc(3*max(M, N)* sizeof(CHAMELEON_Complex64_t));
    CHAMELEON_Complex64_t *A2 = NULL;
    CHAMELEON_Complex64_t *U  = NULL;
    CHAMELEON_Complex64_t *VT = NULL;
    CHAM_desc_t *T;

    /* Check if unable to allocate memory */
    if ( (!A1) || (!S1) || (!S2) || (!work) ) {
        free(A1); free(work);
        free(S1); free(S2);
        printf("Out of Memory \n ");
        return -2;
    }

    /* TODO: check problem with workspace!!! */
    CHAMELEON_Alloc_Workspace_zgesvd(M, N, &T, 1, 1);

    /*----------------------------------------------------------
    *  TESTING ZGESVD
    */
    /* Initialize A1 */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, M, N,
                         morse_lapack_const(ChamDistUniform), ISEED,
                         morse_lapack_const(ChamNonsymPosv), S1, mode, rcond,
                         dmax, M, N,
                         morse_lapack_const(ChamNoPacking), A1, LDA, work );
    free(work);

    /* Copy A1 for check */
    if ( (jobu == ChamVec) && (jobvt == ChamVec) ) {
        A2 = (CHAMELEON_Complex64_t *)malloc(LDA*N*sizeof(CHAMELEON_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N, A1, LDA, A2, LDA);
    }
    if ( jobu == ChamVec ) {
        U = (CHAMELEON_Complex64_t *)malloc(M*M*sizeof(CHAMELEON_Complex64_t));
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', M, M, 0., 1., U, M);
    }
    if ( jobvt == ChamVec ) {
        VT = (CHAMELEON_Complex64_t *)malloc(N*N*sizeof(CHAMELEON_Complex64_t));
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', N, N, 0., 1., VT, N);
    }

    /* CHAMELEON ZGESVD */
    INFO = CHAMELEON_zgesvd(jobu, jobvt, M, N, A1, LDA, S2, T, U, M, VT, N);
    if( INFO != 0 ){
        printf(" CHAMELEON_zgesvd returned with error code %d\n",INFO);
        goto fin;
    }

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZGESVD ROUTINE -------  \n");
    printf("        Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");


    /* Check the orthogonality, reduction and the singular values */
    if ( jobu == ChamVec )
        info_orthou = check_orthogonality(ChamLeft, M, M, U, M, eps);

    if ( jobvt == ChamVec )
        info_orthovt = check_orthogonality(ChamRight, N, N, VT, N, eps);

    if ( (jobu == ChamVec) && (jobvt == ChamVec) )
        info_reduction = check_reduction(M, N, S2, A2, LDA, U, M, VT, N, eps);

    info_solution = check_solution(minMN, S1, S2, eps);

    if ( (info_solution == 0) & (info_orthou == 0) &
         (info_orthovt == 0) & (info_reduction == 0) ) {
        if (M >= N) {
           printf("***************************************************\n");
           printf(" ---- TESTING ZGESVD .. M >= N ........... PASSED !\n");
           printf("***************************************************\n");
        }
        else {
           printf("***************************************************\n");
           printf(" ---- TESTING ZGESVD .. M < N ............ PASSED !\n");
           printf("***************************************************\n");
        }
    }
    else {
        if (M >= N) {
           printf("************************************************\n");
           printf(" - TESTING ZGESVD .. M >= N .. FAILED !\n");
           printf("************************************************\n");
        }
        else {
           printf("************************************************\n");
           printf(" - TESTING ZGESVD .. M < N .. FAILED !\n");
           printf("************************************************\n");
        }
    }

fin:
    if ( A2 != NULL ) free(A2);
    if ( U  != NULL ) free(U);
    if ( VT != NULL ) free(VT);
    free(A1); free(S1); free(S2);
    CHAMELEON_Dealloc_Workspace(&T);

    return 0;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of U VT
 */
static int check_orthogonality(cham_side_t side, int M, int N, CHAMELEON_Complex64_t *Q, int LDQ, double eps)
{
    double  done =  1.0;
    double  mdone  = -1.0;
    double  normQ, result;
    int     info_ortho;
    int     minMN = min(M, N);
    double *work = (double *)malloc(minMN*sizeof(double));

    /* Build the idendity matrix */
    CHAMELEON_Complex64_t *Id = (CHAMELEON_Complex64_t *) malloc(minMN*minMN*sizeof(CHAMELEON_Complex64_t));
    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN);

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, done, Q, LDQ, mdone, Id, N);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,   M, N, done, Q, LDQ, mdone, Id, M);

    normQ = LAPACKE_zlansy_work(LAPACK_COL_MAJOR, morse_lapack_const(ChamInfNorm), 'U', minMN, Id, minMN, work);

    if (getenv("CHAMELEON_TESTING_VERBOSE"))
        printf( "||Q||_oo=%e\n", normQ );

    result = normQ / (M * eps);
    if (side == ChamLeft)
    {
        printf(" ======================================================\n");
        printf(" ||Id-U'*U||_oo / (M*eps)            : %e \n",  result );
        printf(" ======================================================\n");
    }
    else
    {
        printf(" ======================================================\n");
        printf(" ||Id-VT'*VT||_oo / (M*eps)          : %e \n",  result );
        printf(" ======================================================\n");
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Orthogonality is suspicious ! \n");
        info_ortho = 1;
    }
    else {
        printf("-- Orthogonality is CORRECT ! \n");
        info_ortho = 0;
    }

    free(work); free(Id);

    return info_ortho;
}

/*------------------------------------------------------------
 *  Check the bidiagonal reduction
 */
static int check_reduction(int M, int N, double *S, CHAMELEON_Complex64_t *A, int LDA,
                           CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *VT, int LDVT, double eps )
{
    CHAMELEON_Complex64_t zone  =  1.0;
    CHAMELEON_Complex64_t mzone = -1.0;
    double Anorm, Rnorm, result;
    int info_reduction;
    int i;
    int maxMN = max(M, N);
    int minMN = min(M, N);

    CHAMELEON_Complex64_t *TEMP     = (CHAMELEON_Complex64_t *)malloc(minMN*minMN*sizeof(CHAMELEON_Complex64_t));
    CHAMELEON_Complex64_t *Residual = (CHAMELEON_Complex64_t *)malloc(M    *N    *sizeof(CHAMELEON_Complex64_t));
    double *work = (double *)malloc(maxMN*sizeof(double));

    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, morse_lapack_const(ChamUpperLower), M, N, A, LDA, Residual, M);

    if ( M >= N ) {
        /* Compute TEMP =  SIGMA * Vt */
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', N, N, VT, LDVT, TEMP, N);
        for (i = 0; i < minMN; i++){
            cblas_zdscal(N, S[i], TEMP + i, N);
        }

        /* Compute Residual = A - U * (SIGMA * VT) */
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N,
                    CBLAS_SADDR(mzone), U,        LDU,
                                        TEMP,     minMN,
                    CBLAS_SADDR(zone),  Residual, M);
    }
    else {
        /* Compute TEMP =  U * SIGMA */
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, M, U, LDU, TEMP, M);
        for (i = 0; i < minMN; i++){
            cblas_zdscal(M, S[i], TEMP + i*M, 1);
        }

        /* Compute Residual = A - (U * SIGMA) * VT */
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, M,
                    CBLAS_SADDR(mzone), TEMP,     M,
                                        VT,       LDVT,
                    CBLAS_SADDR(zone),  Residual, M);
    }

    /* Compute the norms */
    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(ChamOneNorm), M, N, Residual, M,   work);
    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(ChamOneNorm), M, N, A,        LDA, work);

    if (getenv("CHAMELEON_TESTING_VERBOSE"))
        printf( "||A||_oo=%e\n||A - U*SIGMA*VT||_oo=%e\n", Anorm, Rnorm );

    result = Rnorm / ( Anorm * maxMN * eps);
    printf(" ======================================================\n");
    printf(" ||A-U*SIGMA*V'||_oo/(||A||_oo.N.eps) : %e \n",  result );
    printf(" ======================================================\n");

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Reduction is suspicious ! \n");
        info_reduction = 1;
    }
    else {
        printf("-- Reduction is CORRECT ! \n");
        info_reduction = 0;
    }

    free(TEMP);
    free(Residual);
    free(work);

    return info_reduction;
}

/*------------------------------------------------------------
 *  Check the eigenvalues
 */
static int check_solution(int N, double *E1, double *E2, double eps)
{
    int info_solution, i;
    double resid;
    double maxtmp;
    double maxel = fabs( fabs(E1[0]) - fabs(E2[0]) );
    double maxeig = max( fabs(E1[0]), fabs(E2[0]) );
    for (i = 1; i < N; i++){
        resid   = fabs(fabs(E1[i])-fabs(E2[i]));
        maxtmp  = max(fabs(E1[i]), fabs(E2[i]));

        /* Update */
        maxeig = max(maxtmp, maxeig);
        maxel  = max(resid,  maxel );
    }

    maxel = maxel / (maxeig * N * eps);
    printf(" ======================================================\n");
    printf(" | S - singularcomputed | / (|S| * N * eps) : %e \n",  maxel );
    printf(" ======================================================\n");

    if ( isnan(maxel) || isinf(maxel) || (maxel > 100) ) {
        printf("-- The singular values are suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The singular values are CORRECT ! \n");
        info_solution = 0;
    }
    return info_solution;
}
