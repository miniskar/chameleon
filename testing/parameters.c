/**
 *
 * @file parameters.c
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
 * @date 2022-02-22
 *
 */
#include "testings.h"

/**
 * @brief Defines all the parameters of the testings
 */
static parameter_t parameters[] = {
    /* Name, helper, shname, flags, has_arg, psize, valtype, value, vallist, read, sprint */
    { "id", "Id of the run", 0, PARAM_OUTPUT, 0, 3, TestValInt, {0}, NULL, NULL, sprint_int },

    { NULL, "Options", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "help",     "Show this help",                           'h', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "human",    "Enable human readable mode",               'H', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "niter",    "Perform multiple iteration per test",      'l', PARAM_OPTION, 1, 0, TestValInt, {1}, NULL, pread_int, sprint_int },
    { "nowarmup", "Disable the warmup run to load libraries", -31, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { "check",    "Enable checking of the result",            'c', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "trace",    "Enable the trace generation",              -30, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "mtxfmt",   "Change the way the matrix is stored (0: global, 1: tiles, 2: OOC)", -32, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 1, 6, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "profile",  "Display the kernel profiling",             -33, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "forcegpu", "Force kernels on GPU",                     -34, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "async",    "Switch to the Async interface",                        's', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "splitsub", "Split the task submission and execution stages",       'S', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "generic",  "Switch to the non optimized generic algorithms",       -35, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "api",      "Select the API to test (0: Descriptors, 1: Standard)", -36, PARAM_OPTION, 1, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
#endif

    { NULL, "Machine parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "threads", "Number of CPU workers per node",      't', PARAM_OPTION | PARAM_OUTPUT, 1, 7, TestValInt, {-1}, NULL, pread_int, sprint_int },
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { "gpus",    "Number of GPU workers per node",      'g', PARAM_OPTION | PARAM_OUTPUT, 1, 4, TestValInt, { 0}, NULL, pread_int, sprint_int },
    { "P",       "Rows (P) in the PxQ process grid",    'P', PARAM_OPTION | PARAM_OUTPUT, 1, 2, TestValInt, { 1}, NULL, pread_int, sprint_int },
    { "Q",       "Columns (Q) in the PxQ process grid", 'Q', PARAM_OUTPUT,                1, 2, TestValInt, { 1}, NULL, pread_int, sprint_int },
#endif

    { NULL, "Main input parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "op",   "Operation to test/time ('-o help' to get the list of functions)", 'o', PARAM_OPTION | PARAM_OUTPUT, 1, 1, TestString, {0}, NULL, pread_string, sprint_string },
    { "file", "Input file",                                                      'f', PARAM_OPTION,                1, 1, TestString, {0}, NULL, pread_string, sprint_string },

    { NULL, "Matrix definition parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "m",    "Dimension M of the operation",    'm', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "n",    "Dimension N of the operation",    'n', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "k",    "Dimension K of the operation",    'k', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "nrhs", "Dimension NRHS of the operation", 'r', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "nb", "Tile size nb",       'b', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { "ib", "Inner tile size ib", 'i', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 2, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "l1", "Size of the first level of recursion",  '1', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "l2", "Size of the second level of recursion", '2', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "l3", "Size of the third level of recursion",  '3', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
#endif

    { "lda", "Leading dimension of the matrix A", 'A', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "ldb", "Leading dimension of the matrix B", 'B', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "ldc", "Leading dimension of the matrix C", 'C', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "seedA", "Seed for the matrix A random generation", 'X', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "seedB", "Seed for the matrix B random generation", 'Y', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "seedC", "Seed for the matrix C random generation", 'Z', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },

    { NULL, "Matrix generation numerical parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "bump",  "Bump value to make a matrix diagonal dominant",           'z', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },
    { "mode",  "Mode that specifies the eigen/singular values in xlatms", -30, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2,  4, TestValInt,       {0}, NULL, pread_int,       sprint_int    },
    { "cond",  "Conditional number of the matrix used by xlatms",         -31, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValDouble,    {0}, NULL, pread_double,    sprint_double    },

    { NULL, "Operation specific parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "trans",  "Value of the trans parameter",  -11, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transA", "Value of the transA parameter", -12, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transB", "Value of the transB parameter", -13, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "uplo",   "Value of the uplo parameter",   -14, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestUplo,     {0}, NULL, pread_uplo,  sprint_uplo  },
    { "diag",   "Value of the diag parameter",   -15, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestDiag,     {0}, NULL, pread_diag,  sprint_diag  },
    { "side",   "Value of the side parameter",   -16, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestSide,     {0}, NULL, pread_side,  sprint_side  },
    { "norm",   "Value of the norm parameter",   -17, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestNormtype, {0}, NULL, pread_norm,  sprint_norm  },

    { NULL, "Operation specific scalar", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "alpha", "Value of the scalar alpha",                       'x', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },
    { "beta",  "Value of the scalar beta",                        'y', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },

#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { NULL, "QR/LQ parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "qra",    "Size of TS domain (=RH for householder trees)",           -20, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "qrp",    "Size of high level tree for distributed",                 -21, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "llvl",   "Tree used for low level reduction insides nodes",         -22, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "hlvl",   "Tree used for high level reduction between nodes",        -23, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "domino", "Enable/Disable the domino between upper and lower trees", -24, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 6, TestValInt, {0}, NULL, pread_int, sprint_int },
#endif

    { "tsub",          "Graph submission time in s",             999, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "time",          "Time in s",                             1000, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "gflops",        "GFlop/s",                               1001, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "RETURN",        "Result of the testing: SUCCESS/FAILED", 1002, PARAM_OUTPUT, 2,  7, TestValInt,    {0}, NULL, pread_int,    sprint_check  },
    { "||Ax-b||",      "Norm of the residual",                  1003, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||A-fact(A)||", "Norm of the residual",                  1004, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||A||",         "Norm of the matrix A",                  1005, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||B||",         "Norm of the matrix B",                  1006, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||C||",         "Norm of the matrix C",                  1007, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||R||",         "Residual norm",                         1008, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||b||",         "Norm of the vector b",                  1009, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||x||",         "Norm of the vector x",                  1010, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||Ax-b||/N/eps/(||A||||x||+||b||", "",                   1011, PARAM_OUTPUT, 2, 22, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||I-QQ'||",     "Orthonormality of Q",                   1012, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
};

void print_usage( const char* prog_name )
{
    int rank = 0;
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    if ( CHAMELEON_Initialized() ) {
        rank = CHAMELEON_Comm_rank();
    }
#endif

    if (rank == 0) {
        parameter_t *param = parameters;
        int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );
        printf( "Usage:\n"
                "  %s -o|--op operation_name [options]\n"
                "  %s -f|--file input_file [options]\n",
                prog_name, prog_name );

        for (i=0; i<nbparams; i++, param++) {
            char str[STR_MAX_LENGTH];

            /* This is not an option, we skip it */
            if ( !(param->flags & PARAM_OPTION) ) {
                continue;
            }

            /* This is an option header */
            if ( param->name == NULL ) {
                printf( "\n  %s:\n", param->helper );
                continue;
            }

            if ( param->shname > 0 ) {
                snprintf( str, STR_MAX_LENGTH, "-%c, --%s",
                          param->shname, param->name );
            }
            else {
                snprintf( str, STR_MAX_LENGTH, "    --%s",
                          param->name );
            }

            /* If an argument is needed, add " x" */
            if ( param->has_arg > 0 ) {
                int len = strlen(str);
                assert( len < (STR_MAX_LENGTH-3) );

                str[ len   ] = ' ';
                str[ len+1 ] = 'x';
                str[ len+2 ] = '\0';
            }
            printf( "    %-23s %s\n",
                    str, param->helper );
        }

        printf( "\n"
                "For example: %s -H -o gemm -t 2 -m 2000 -n 2000 -k 2000 -b 200\n"
                "  will run one gemm with three matrices of size 2000x2000 each and a tile size of 200.\n"
                "  The output will be in the human readable format\n"
                "\n", prog_name );
#if !defined(CHAMELEON_TESTINGS_VENDOR)
        printf( "Remarks about timing:\n"
                "  Timings are reported respectively as 'tsub' for the graph submission time, and 'time'\n"
                "  for the execution time.\n"
                "  By default the synchronous tile interface is used to perform the timings. 'tsub' is null.\n"
                "  If the --async option is enabled, then the asynchronous interface is called. 'tsub' reports\n"
                "  the task submission time, and 'time' the execution time that includes 'tsub'.\n"
                "  If the --splitsub option is enabled, then the asynchronous interface is called and task\n"
                "  submission is fully performed before starting the computation. 'tsub' reports the\n"
                "  task submission time, and 'time' the execution time excluding 'tsub'.\n"
                "  Note that the 'gflops' field is always computed with 'time'\n" );
#endif
    }
}

/**
 ********************************************************************************
 *
 * @brief Get the list of values associated to a given parameter
 *
 *******************************************************************************
 *
 * @param[in] name
 *          The name of the parameter we are interested in.
 *
 * @return NULL if no parameter exists with this name, otherwise the pointer to
 * the list of values associated to this parameter.
 *
 *******************************************************************************
 */
vallist_t *
parameters_getlist( const char *name )
{
    parameter_t *param = parameters_getbyname( name );
    if ( param == NULL ) {
        return NULL;
    }
    else {
        return param->vallist;
    }
}

/**
 ********************************************************************************
 *
 * @brief Parses a list in form A1, A2, ..., An and insert the values in an
 * argument list.
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter associated to the list.
 *          On exit, the list of values are added to the parameter list of
 *          possible values.
 *
 * @param[in] liststr
 *          The string that holds the list
 *
 *******************************************************************************
 */
void
parameters_read_list( parameter_t *param,
                      const char  *liststr )
{
    const char *delim = ", \n";
    char *str = strdup( liststr );
    char *token, *saveptr;
    vallist_t *previous, *current;

    /* Initialize the list items */
    previous = NULL;
    current  = param->vallist;

    /* Move to the end of the list if some parameters have already been registered */
    while( current != NULL ) {
        previous = current;
        current  = current->next;
    }

    token = strtok_r( str, delim, &saveptr );
    while ( token != NULL ) {
        assert( current == NULL );
        current = calloc( 1, sizeof(vallist_t) );

        /* Read the value */
        current->value = param->read( token );

        /* Insert at the end of the list */
        if ( previous != NULL ) {
            previous->next = current;
        }
        else {
            /* Nothing was in the list */
            param->vallist = current;
        }

        previous = current;
        current  = NULL;

        /* Move to the next token */
        token = strtok_r( NULL, delim, &saveptr );
    }

    free( str );
}

/**
 ********************************************************************************
 *
 * @brief Parses a list in form start:end[:step] and inserts the values in an
 * argument list.
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter associated to the list.
 *          On exit, the range of values are added to the parameter list of
 *          possible values.
 *
 * @param[in] rangestr
 *          The string that holds the range
 *
 * @param[in] min
 *          The minimum value available
 *
 * @param[in] max
 *          The maximum value available
 *
 *******************************************************************************
 */
void
parameters_read_intrange( parameter_t *param,
                          const char  *rangestr,
                          int min, int max )
{
    int start, end, step, count;
    vallist_t *previous, *current;

    max = (max == -1) ? INT32_MAX : max;

    count = sscanf( rangestr, "%d:%d:%d", &start, &end, &step );
    if ( count < 2 ) {
        fprintf(stderr, "Incorrect range syntax (%s): data skipped\n", rangestr );
        return;
    }
    else if (count == 2) {
        step = 1;
    }

    /* Check the range */
    if ( (start < min) || (start > max) ||
         (end   < min) || (end   > max) )
    {
        /* Try to shift to 0 to see if now we fit */
        start += min;
        end   += min;
        if ( (start < min) || (start > max) ||
             (end   < min) || (end   > max) )
        {
            fprintf( stderr, "Incorrect range values outside the possible ranges [%d:%d]",
                     min, max );
            if ( min > 0 ) {
                fprintf( stderr, " or [%d:%d]\n", 0, max-min );
            }
            else {
                fprintf( stderr, "\n" );
            }
        }
    }

    /* Initialize the list items */
    previous = NULL;
    current  = param->vallist;

    /* Move to the end of the list if some parameters have already been registered */
    while( current != NULL ) {
        previous = current;
        current  = current->next;
    }

    while ( start <= end ) {
        assert( current == NULL );
        current = calloc( 1, sizeof(vallist_t) );

        /* Read the value */
        current->value.ival = start;

        /* Insert at the end of the list */
        if ( previous != NULL ) {
            previous->next = current;
        }
        else {
            /* Nothing was in the list */
            param->vallist = current;
        }

        previous = current;
        current  = NULL;

        start += step;
    }
}

/**
 ********************************************************************************
 *
 * @brief Wrapper to parse a list or range of values associated to a parameter.
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter associated to the list.
 *          On exit, the range of values are added to the parameter list of
 *          possible values.
 *
 * @param[in] values
 *          The string that holds the range of list of values
 *
 *******************************************************************************
 */
void
parameters_read( parameter_t *param,
                 const char  *values )
{
    int range = ( strchr( values, ':' ) != NULL );

    /* If we have a ranged of integer values */
    if ( range )
    {
        switch ( param->valtype ) {
        case TestValInt:
            parameters_read_intrange( param, values, 0, -1 );
            break;
        case TestTrans:
            parameters_read_intrange( param, values, ChamNoTrans, ChamConjTrans );
            break;
        case TestUplo:
            parameters_read_intrange( param, values, ChamUpper, ChamUpperLower );
            break;
        case TestDiag:
            parameters_read_intrange( param, values, ChamNonUnit, ChamUnit );
            break;
        case TestSide:
            parameters_read_intrange( param, values, ChamLeft, ChamRight );
            break;
        case TestNormtype:
            parameters_read_intrange( param, values, ChamOneNorm, ChamMaxNorm );
            break;
        default:
            fprintf( stderr, "parameters_read: range is not available for this datatype (%d)\n",
                     param->valtype );
        }
        return;
    }

    parameters_read_list( param, values );
}

/**
 ********************************************************************************
 *
 * @brief Generic function to add value(s) to a given parameter
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter that will receive the value
 *          On exit, the value(s) (switch, list, range, ...) is/are added to the
 *          parameter list of possible values
 *
 * @param[in] values
 *          The string that holds the values (list, range, NULL if switch)
 *
 *******************************************************************************
 */
void
parameters_addvalues( parameter_t *param,
                      const char  *values )
{
    if ( param->has_arg == 0 ) {
        param->value.ival = 1;
    }
    else if ( param->has_arg == 1 ) {
        param->value = param->read( values );
    }
    else {
        parameters_read( param, values );
    }
}

/**
 ********************************************************************************
 *
 * @brief Parses an input test file.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The name of the input file.
 *
 *******************************************************************************
 */
void
parameters_read_file( const char  *filename )
{
    FILE        *fp;
    const char  *delim = " =";
    char        *saveptr;
    char        *line_read, *line;
    char        *name, *values;
    size_t       len = 256;
    parameter_t *param;

    fp = fopen( filename, "r" );
    if ( fp == NULL ) {
        fprintf( stderr, "Error reading input file %s\n", filename );
        perror("fopen");
        exit(1);
    }

    len = 256;
    line_read = malloc( len * sizeof( char ) );

    while ( getline( &line_read, &len, fp ) != -1 )
    {
        line = line_read;

        /* Ignores comments and empty lines */
        if ( (line[0] == '#' ) ||
             (line[0] == '\n') )
        {
            continue;
        }

        /* Removes possible extra spaces */
        while ( line[0] == ' ' ) {
            line++;
        }

        /* Reads the parameter name and values */
        name   = strtok_r( line, delim, &saveptr );
        values = strtok_r( NULL, "",    &saveptr );

        /* Goes for the listed values */
        while ( (values[0] == ' ') ||
                (values[0] == '=') )
        {
            values++;
        }

        param = parameters_getbyname( name );
        if ( param == NULL ) {
            fprintf( stderr, "Parameter %s is not know. We skip it\n", name );
            continue;
        }
        parameters_addvalues( param, values );
    }

    free(line_read);
    fclose(fp);
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
parameters_compute_q( int p )
{
    parameter_t *param;
    int np = CHAMELEON_Comm_size();

    if ( (np % p) != 0 ) {
        fprintf( stderr, "ERROR: The number of processes (%d) must be a multiple of P (%d)\n", np, p );
        exit(1);
    }

    param = parameters_get( 'Q' );
    param->value.ival = np / p;
    return param->value.ival;
}
#endif

void
parameters_getopt_init( char           *optstring,
                        struct option **longopts )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );
    int nboptions = 0;
    int strpos = 0;

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( !(param->flags & PARAM_OPTION) ||
             (param->name == NULL) )
        {
            continue;
        }

        nboptions++;

        if ( param->shname < 0 ) {
            continue;
        }

        optstring[strpos] = param->shname;
        strpos++;
        assert( strpos < STR_MAX_LENGTH );

        if ( param->has_arg > 0 ) {
            optstring[strpos] = ':';
            strpos++;
            assert( strpos < STR_MAX_LENGTH );
        }
    }
    optstring[strpos] = '\0';

    /* Now, let's generate the long opt if needed */
#if defined(CHAMELEON_HAVE_GETOPT_LONG)
    if ( longopts != NULL ) {
        struct option *opt;
        *longopts = calloc( nboptions+1, sizeof( struct option ) );

        opt = *longopts;
        param = parameters;

        for ( i=0; i<nboptions; i++, opt++, param++ ) {

            /* Look for a valid option */
            while ( !(param->flags & PARAM_OPTION) ||
                    (param->name == NULL) )
            {
                param++;
            }

            opt->name    = param->name;
            opt->has_arg = ( param->has_arg > 0 ) ? 1 : 0;
            opt->flag    = NULL;
            opt->val     = param->shname;
        }
    }
#endif
}

parameter_t *
parameters_get( int shname )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( shname == param->shname ) {
            return param;
        }
    }

    fprintf( stderr, "parameters_get could not find parameter %d(%c)\n", shname, shname );
    return NULL;
}

int
parameters_getvalue_int( const char *name )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            continue;
        }

        if ( param->has_arg > 1 ) {
            fprintf( stderr, "parameters_getvalue_int should not be called with parameter %s\n", name );
            return -1;
        }

        if ( param->valtype != TestValInt ) {
            fprintf( stderr, "parameters_getvalue_int has been called with a non integer parameter (%s)\n", name );
            return -1;
        }

        return param->value.ival;
    }

    fprintf( stderr, "parameters_getvalue_int could not find parameter %s\n", name );
    return -1;
}

char *
parameters_getvalue_str( const char *name )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            continue;
        }

        if ( param->has_arg > 1 ) {
            fprintf( stderr, "parameters_getvalue_str should not be called with parameter %s\n", name );
            return NULL;
        }

        if ( param->valtype != TestString ) {
            fprintf( stderr, "parameters_getvalue_str has been called with a non string parameter (%s)\n", name );
            return NULL;
        }

        return param->value.str;
    }

    fprintf( stderr, "parameters_getvalue_str could not find parameter %s\n", name );
    return NULL;
}

parameter_t *
parameters_getbyname( const char *name )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            continue;
        }

        /* if ( param->has_arg < 2 ) { */
        /*     fprintf( stderr, "parameters_getbyname should not be called with parameter %s\n", name ); */
        /*     return NULL; */
        /* } */

        return param;
    }

    fprintf( stderr, "parameters_getbyname could not find parameter %s\n", name );
    return NULL;
}

void parameters_parser( int argc, char **argv )
{
    int opt;
    char optstring[STR_MAX_LENGTH];
    struct option *longopts = NULL;
    parameter_t *param;

    parameters_getopt_init( optstring, &longopts );

#if defined(CHAMELEON_HAVE_GETOPT_LONG)
    while ((opt = getopt_long(argc, argv, optstring, longopts, NULL)) != -1)
#else
    while ((opt = getopt(argc, argv, optstring)) != -1)
#endif
    {
        switch(opt) {
        case 'h':
            print_usage(argv[0]);
            exit(0);

        case '?': /* error from getopt[_long] */
            exit(1);
            break;

        default:
            param = parameters_get( opt );
            if ( param == NULL ) {
                print_usage(argv[0]);
                exit(1);
            }
            parameters_addvalues( param, optarg );
        }
    }

    if ( longopts != NULL ) {
        free( longopts );
    }

#if !defined(CHAMELEON_TESTINGS_VENDOR)
    /* Force Async if splitsub is enabled */
    {
        int splitsub = parameters_getvalue_int( "splitsub" );

        if ( splitsub ) {
            param = parameters_get( 's' );
            if ( param == NULL ) {
                print_usage(argv[0]);
                exit(1);
            }
            parameters_addvalues( param, NULL );

#if defined(CHAMELEON_RUNTIME_SYNC)
            fprintf( stderr, "Spliting the submission and the execution stages is not possible when the option CHAMELEON_RUNTIME_SYNC is enabled\n" );
            exit(0);
#endif
        }
    }
#endif
}

void
parameters_destroy()
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );
    vallist_t *current, *next;

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->has_arg < 2 ) {
            continue;
        }

        current = param->vallist;
        while ( current != NULL )
        {
            next = current->next;
            free( current );
            current = next;
        }
    }
    return;
}