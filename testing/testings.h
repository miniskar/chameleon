/**
 *
 * @file testings.h
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 ***
 *
 * @brief Chameleon auxiliary routines for testing structures
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-22
 *
 */
#ifndef _testings_h_
#define _testings_h_

#include <chameleon.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#if defined(CHAMELEON_HAVE_GETOPT_LONG)
#include <getopt.h>
#else
struct option;
#endif

#define STR_MAX_LENGTH 256

typedef enum valtype_ {
    TestValInt,
    TestValFloat,
    TestValDouble,
    TestValComplex32,
    TestValComplex64,
    TestTrans,
    TestUplo,
    TestDiag,
    TestSide,
    TestJob,
    TestNormtype,
    TestString,
} valtype_e;

/* Define to avoid conversion */
#define TestValFixdbl TestValDouble

union val_u {
    int                   ival;
    cham_trans_t          trans;
    cham_uplo_t           uplo;
    cham_diag_t           diag;
    cham_side_t           side;
    cham_job_t            job;
    cham_normtype_t       ntype;
    CHAMELEON_Complex64_t zval;
    CHAMELEON_Complex32_t cval;
    double                dval;
    float                 sval;
    char                 *str;
};
typedef union val_u val_t;

/**
 * @brief Defines a values that can be taken by any parameter
 */
struct vallist_s;
typedef struct vallist_s vallist_t;

struct vallist_s {
    val_t      value; /**< Value of the parameter                       */
    vallist_t *next;  /**< Pointer to the next parameter                */
};

#define PARAM_OPTION (1 << 0)
#define PARAM_INPUT  (1 << 1)
#define PARAM_OUTPUT (1 << 2)

typedef val_t (*read_fct_t)( const char * );
typedef char *(*sprint_fct_t)( val_t, int, int, char * );

/**
 * @brief Defines the possible testing parameters used to define the tests
 */
typedef struct parameter_s {
    const char   *name;    /**< Name of the option                   */
    const char   *helper;  /**< Helper string for the usage function */
    int           shname;  /**< The associated short option          */
    int           flags;   /**< Flags to define if it is: an option, an input, an output */
    int           has_arg; /**< Defines the number of arguments: 0- no arguments, 1- A single argument for all test, 2- multiple arguments possible */
    int           psize;   /**< Number of characters to printed in the human readable format */
    valtype_e     valtype; /**< Type of the argument if has_arg > 1, unread otherwise   */
    val_t         value;   /**< Default value if has_arg < 2                            */
    vallist_t    *vallist; /**< List of values for input parameters (has_arg >= 2)      */
    read_fct_t    read;    /**< Parser for one argument                                 */
    sprint_fct_t  sprint;  /**< Sprint function for one argument                        */
} parameter_t;

/**
 * @brief Defines a single parameter for one test case. The parameters are
 * chained together through a NULL terminated list.
 */
struct run_arg_s;
typedef struct run_arg_s run_arg_t;

struct run_arg_s {
    parameter_t *param; /**< Pointer to the testing parameter information */
    val_t        value; /**< Value of the parameter                       */
    run_arg_t   *next;  /**< Pointer to the next parameter                */
};

typedef struct run_arg_list_s {
    run_arg_t *head;
    run_arg_t *tail;
} run_arg_list_t;

/**
 * @brief Defines a single run to perform. The tests are chained together by a
 * null terminated list.
 */
struct run_list_elt_s;
typedef struct run_list_elt_s run_list_elt_t;

struct run_list_elt_s {
    run_arg_list_t  args; /**< List of parameters defining the run */
    run_list_elt_t *next; /**< Pointer to the next run description */
};

typedef struct run_list_s {
    run_list_elt_t *head; /**< List of parameters defining the run */
    run_list_elt_t *tail; /**< Pointer to the next run description */
} run_list_t;

/**
 * @brief Defines the possible testing available for the given precision
 *
 * This is a NULL terminated list.
 */
struct testing_;
typedef struct testing_ testing_t;
typedef int (*test_fct_t)( run_arg_list_t *, int );

typedef struct testing_ {
    const char  *name;        /**< Name of the operation tested                               */
    const char  *helper;      /**< Helper of the function tested for the usage                */
    const char **params;      /**< Parameters used by the operation                           */
    const char **output;      /**< Additional output parameters generated by the function     */
    const char **outchk;      /**< Check specific output parameters generated by the function */
    test_fct_t   fptr_std;    /**< Function performing the standard test                      */
    test_fct_t   fptr_desc;   /**< Function performing the descriptor test                    */
    testing_t   *next;        /**< Pointer to following test                                  */
} testing_t;

val_t pread_int      ( const char *str );
val_t pread_float    ( const char *str );
val_t pread_double   ( const char *str );
val_t pread_complex32( const char *str );
val_t pread_complex64( const char *str );
val_t pread_trans    ( const char *str );
val_t pread_uplo     ( const char *str );
val_t pread_diag     ( const char *str );
val_t pread_side     ( const char *str );
val_t pread_job      ( const char *str );
val_t pread_norm     ( const char *str );
val_t pread_string   ( const char *str );

#define pread_fixdbl pread_double

void print_usage( const char* prog_name );

char *sprint_int      ( val_t val, int human, int nbchar, char *str_in );
char *sprint_float    ( val_t val, int human, int nbchar, char *str_in );
char *sprint_double   ( val_t val, int human, int nbchar, char *str_in );
char *sprint_complex32( val_t val, int human, int nbchar, char *str_in );
char *sprint_complex64( val_t val, int human, int nbchar, char *str_in );
char *sprint_trans    ( val_t val, int human, int nbchar, char *str_in );
char *sprint_uplo     ( val_t val, int human, int nbchar, char *str_in );
char *sprint_diag     ( val_t val, int human, int nbchar, char *str_in );
char *sprint_side     ( val_t val, int human, int nbchar, char *str_in );
char *sprint_job      ( val_t val, int human, int nbchar, char *str_in );
char *sprint_norm     ( val_t val, int human, int nbchar, char *str_in );
char *sprint_string   ( val_t val, int human, int nbchar, char *str_in );
char *sprint_check    ( val_t val, int human, int nbchar, char *str_in );

#define sprint_fixdbl sprint_double

float                 testing_salea();
double                testing_dalea();
CHAMELEON_Complex32_t testing_calea();
CHAMELEON_Complex64_t testing_zalea();

const run_arg_t *run_arg_get_byname( const run_arg_list_t *arglist, const char *name );

int                   run_arg_get_int      ( run_arg_list_t *arglist, const char *name, int                   defval );
float                 run_arg_get_float    ( run_arg_list_t *arglist, const char *name, float                 defval );
double                run_arg_get_double   ( run_arg_list_t *arglist, const char *name, double                defval );
CHAMELEON_Complex32_t run_arg_get_complex32( run_arg_list_t *arglist, const char *name, CHAMELEON_Complex32_t defval );
CHAMELEON_Complex64_t run_arg_get_complex64( run_arg_list_t *arglist, const char *name, CHAMELEON_Complex64_t defval );
cham_trans_t          run_arg_get_trans    ( run_arg_list_t *arglist, const char *name, cham_trans_t          defval );
cham_uplo_t           run_arg_get_uplo     ( run_arg_list_t *arglist, const char *name, cham_uplo_t           defval );
cham_diag_t           run_arg_get_diag     ( run_arg_list_t *arglist, const char *name, cham_diag_t           defval );
cham_side_t           run_arg_get_side     ( run_arg_list_t *arglist, const char *name, cham_side_t           defval );
cham_job_t            run_arg_get_job      ( run_arg_list_t *arglist, const char *name, cham_job_t            defval );
cham_normtype_t       run_arg_get_ntype    ( run_arg_list_t *arglist, const char *name, cham_normtype_t       defval );

int run_arg_add_int   ( run_arg_list_t *arglist, const char *name, int    defval );
int run_arg_add_double( run_arg_list_t *arglist, const char *name, double defval );
#define run_arg_add_fixdbl run_arg_add_double

void run_print_header( const testing_t *test, int check, int human );
void run_print_line( const testing_t *test, const run_arg_list_t *arglist,
                     int check, int human, int id );

run_arg_list_t run_arg_list_copy( const run_arg_list_t *arglist );
void run_arg_list_destroy( run_arg_list_t *arglist );

void         parameters_read( parameter_t *param, const char  *values );
void         parameters_read_file( const char  *filename );
parameter_t *parameters_getbyname( const char *name );
void         parameters_addvalues( parameter_t *param, const char  *values );
int          parameters_getvalue_int( const char *name );
int          parameters_compute_q( int p );
parameter_t *parameters_get( int shname );
int          parameters_compute_q( int p );
void         parameters_getopt_init( char *optstring, struct option **longopts );
parameter_t *parameters_get( int shname );
int          parameters_getvalue_int( const char *name );
char *       parameters_getvalue_str( const char *name );
parameter_t *parameters_getbyname( const char *name );
void         parameters_parser( int argc, char **argv );
void         parameters_destroy( );

run_list_t *run_list_generate( const char **params );
void        run_list_destroy( run_list_elt_t *run );

/**
 * @brief Define the data associated to a single run of a testing
 */
struct testing_;
typedef struct testing_ testing_t;
typedef int (*test_fct_t)( run_arg_list_t *, int );

typedef struct testdata_ {
    run_arg_list_t     *args;     /**< The parameters of the test           */
    int                 hres;     /**< The returned value of the test       */
    cham_fixdbl_t       texec;    /**< The execution time of test           */
    cham_fixdbl_t       tsub;     /**< The task submission tome of the test */
    RUNTIME_sequence_t *sequence; /**< The sequence to run the test if splitsub */
    RUNTIME_request_t   request;  /**< The request to run the test if splitsub  */
} testdata_t;

/**
 * @brief Structure to store the read parameters for a quicker access
 */
typedef struct testing_options_ {
    /* Static parameters */
    int api;
    int async;
    int check;
    int forcegpu;
    int generic;
    int gpus;
    int human;
    int mtxfmt;
    int niter;
    int nowarmup;
    int P;
    int profile;
    int splitsub;
    int threads;
    int trace;
    char *file;
    char *op;
    /* Additionnal information to exchange between the main and the testings */
    int run_id;
} testing_options_t;

void       testing_register( testing_t *test );
testing_t *testing_gettest( const char *prog_name, const char *func_name );
void       testing_start( testdata_t *tdata );
void       testing_stop( testdata_t *tdata, cham_fixdbl_t flops );
void       testing_options_init( testing_options_t *options );

#endif /* _testings_h_ */
