/**
 *
 * @file vendor_ztesting.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief MKL CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.1.0
 * @author Alycia Lisito
 * @date 2022-02-25
 * @precisions normal z -> c d s
 *
 */
#include "testings.h"

int main (int argc, char **argv) {

    int ncores, human, i, niter;
    int nowarmup;
    int rc, info, check = 0;
    int run_id = 0;
    char *func_name;
    char *input_file;
    run_list_t *runlist;
    testing_t * test;
    run_list_elt_t *run, *next;

    /* Reads the arguments from command line */
    parameters_parser( argc, argv );
    input_file = parameters_getvalue_str( "file" );
    if ( input_file != NULL ) {
        parameters_read_file( input_file );
        free(input_file);
    }
    ncores    = parameters_getvalue_int( "threads"  );
    human     = parameters_getvalue_int( "human"    );
    func_name = parameters_getvalue_str( "op"       );
    niter     = parameters_getvalue_int( "niter"    );
    nowarmup  = parameters_getvalue_int( "nowarmup" );

    rc = CHAMELEON_Init( ncores, 0 );
    if ( rc != CHAMELEON_SUCCESS ) {
        fprintf( stderr, "CHAMELEON_Init failed and returned %d.\n", rc );
        info = 1;
        goto end;
    }

    /* Set ncores to the right value */
    if ( ncores == -1 ) {
        parameter_t *param;
        param = parameters_get( 't' );
        param->value.ival = CHAMELEON_GetThreadNbr();
    }

    /* Binds the right function to be called and builds the parameters combinations */
    test = testing_gettest( argv[0], func_name );
    free(func_name);
    test_fct_t fptr = test->fptr_std;
    if ( fptr == NULL ) {
        fprintf( stderr, "The vendor API is not available for function %s\n", func_name );
        info = 1;
        goto end;
    }

    /* Generate the cartesian product of the parameters */
    runlist = run_list_generate( test->params );

    /* Executes the tests */
    run_print_header( test, check, human );
    run = runlist->head;

    /* Warmup */
    if ( !nowarmup ) {
        run_arg_list_t copy = run_arg_list_copy( &(run->args) );
        fptr( &copy, check );
        run_arg_list_destroy( &copy );
    }

    /* Perform all runs */
    while ( run != NULL ) {
        for(i=0; i<niter; i++) {
            run_arg_list_t copy = run_arg_list_copy( &(run->args) );
            rc = fptr( &copy, check );

            /* If rc < 0, we skipped the test */
            if ( rc >= 0 ) {
                run_arg_add_int( &copy, "RETURN", rc );
                run_print_line( test, &copy, check, human, run_id );
                run_id++;
                info += rc;
            }
            run_arg_list_destroy( &copy );
        }

        /* Move to next run */
        next = run->next;
        run_list_destroy( run );
        run = next;
    }

    free( runlist );

  end:
    ;/* OpenMP end */
    CHAMELEON_Finalize();
    parameters_destroy();

    return info;
}
