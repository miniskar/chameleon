/**
 *
 * @file starpu/codelet_zgram.c
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram StarPU codelet
 *
 * @version 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2019-04-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgram_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int m, n, mt, nt;
    double *Di;
    int lddi;
    double *Dj;
    int lddj;
    double *D;
    double *A;
    int lda;

    Di = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    Dj = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    D  = (double *)STARPU_MATRIX_GET_PTR(descr[2]);
    A  = (double *)STARPU_MATRIX_GET_PTR(descr[3]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &m, &n, &mt, &nt, &lddi, &lddj, &lda);
    CORE_zgram( uplo,
                m, n, mt, nt,
                Di, lddi,
                Dj, lddj,
                D,
                A, lda);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgram, 4, cl_zgram_cpu_func)

void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din, int lddi,
                        const CHAM_desc_t *Dj, int Djm, int Djn, int lddj,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An, int lda)
{
  struct starpu_codelet *codelet = &cl_zgram;
  void (*callback)(void*) = options->profiling ? cl_zgram_callback : NULL;

  CHAMELEON_BEGIN_ACCESS_DECLARATION;
  CHAMELEON_ACCESS_R(Di, Dim, Din);
  CHAMELEON_ACCESS_R(Dj, Djm, Djn);
  CHAMELEON_ACCESS_R(D, Dm, Dn);
  CHAMELEON_ACCESS_RW(A, Am, An);
  CHAMELEON_END_ACCESS_DECLARATION;

  starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,                      sizeof(int),
        STARPU_VALUE,    &m,                         sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_VALUE,    &mt,                        sizeof(int),
        STARPU_VALUE,    &nt,                        sizeof(int),
        STARPU_R,        RTBLKADDR(Di, double, Dim, Din),
        STARPU_VALUE,    &lddi,                      sizeof(int),
        STARPU_R,        RTBLKADDR(Dj, double, Djm, Djn),
        STARPU_VALUE,    &lddj,                      sizeof(int),
        STARPU_R,        RTBLKADDR(D, double, Dm, Dn),
        STARPU_RW,       RTBLKADDR(A, double, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgram",
#endif
        0);
}
