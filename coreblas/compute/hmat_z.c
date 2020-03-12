/**
 *
 * @file core_zhmat.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2019-2019 Universidad Jaume I. All rights reserved.
 *
 * @brief Chameleon CPU kernel interface from CHAM_tile_t layout to the real one.
 *
 * @version 1.0.0
 * @author Rocio Carratala-Saez
 * @author Mathieu Faverge
 * @date 2019-12-02
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/hmat.h"

/**
 * @brief The hmat interface to the C++ functions
 */
hmat_interface_t hmat_zinterface;

void __hmat_zinit() __attribute__(( constructor ));
void __hmat_zinit() {
    hmat_init_default_interface( &hmat_zinterface, HMAT_DOUBLE_COMPLEX );
}

void __hmat_zfini() __attribute__(( destructor ));
void __hmat_zfini() {
    hmat_zinterface.finalize();
}

int hmat_zpotrf( hmat_matrix_t *A ) {
    hmat_factorization_context_t ctx_facto;
    hmat_factorization_context_init( &ctx_facto );
    ctx_facto.factorization = hmat_factorization_llt;
    hmat_zinterface.factorize_generic( A, &ctx_facto );
    return 0;
}

int hmat_zgetrf( hmat_matrix_t *A ) {
    hmat_factorization_context_t ctx_facto;
    hmat_factorization_context_init( &ctx_facto );
    ctx_facto.factorization = hmat_factorization_lu;
    hmat_zinterface.factorize_generic( A, &ctx_facto );
    return 0;
}

int hmat_zgemm( char transA, char transB, void* alpha, hmat_matrix_t* A,
		hmat_matrix_t* B, void* beta, hmat_matrix_t* C ) {
    return hmat_zinterface.gemm( transA, transB, alpha, A, B, beta, C );
}

int hmat_zgemv( char transA, void* alpha, hmat_matrix_t* A,
		void* B, void* beta, void* C, int nrhs ) {
    return hmat_zinterface.gemm_scalar( transA, alpha, A, B, beta, C, nrhs );
}

int hmat_ztrsm( char side, char uplo, char trans, char diag, int m, int n,
		void* alpha, hmat_matrix_t* A, int is_b_hmat, void* B ) {
    return hmat_zinterface.trsm( side, uplo, trans, diag, m, n, alpha, A, is_b_hmat, B );
}

hmat_matrix_t *hmat_zread( void *buffer ) {
    hmat_buffer_comm_t buffer_struct = { 0, (char*)buffer };
    hmat_matrix_t *hmat = hmat_zinterface.read_struct( (hmat_iostream)(&buffer_comm_read),
                                                       (void *)(&buffer_struct) );
    hmat_zinterface.read_data( hmat, (hmat_iostream)(&buffer_comm_read), (void *)(&buffer_struct) );
    return hmat;
}

size_t hmat_zsize( hmat_matrix_t *hmat ) {
    size_t size = 0;
    hmat_zinterface.write_struct( hmat, (hmat_iostream)(&buffer_comm_size), (void *)(&size) );
    hmat_zinterface.write_data(   hmat, (hmat_iostream)(&buffer_comm_size), (void *)(&size) );
    return size;
}

void hmat_zwrite( hmat_matrix_t *hmat, char *ptr ) {
    hmat_buffer_comm_t buffer_struct = { 0, ptr };
    hmat_zinterface.write_struct( hmat, (hmat_iostream)(&buffer_comm_write), (void *)(&buffer_struct) );
    hmat_zinterface.write_data(   hmat, (hmat_iostream)(&buffer_comm_write), (void *)(&buffer_struct) );
}

void hmat_zdestroy( hmat_matrix_t *hmat ) {
    hmat_zinterface.destroy( hmat );
}
