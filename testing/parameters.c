/**
 *
 * @file parameters.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 ***
 *
 * @brief Chameleon auxiliary routines for testing structures
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 *
 */
#include "testings.h"

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
    ssize_t      nbread;
    parameter_t *param;

    fp = fopen( filename, "r" );
    if ( fp == NULL ) {
        fprintf( stderr, "Error reading input file %s\n", filename );
        perror("fopen");
        exit(1);
    }

    len = 256;
    line_read = malloc( len * sizeof( char ) );

    while ( (nbread = getline( &line_read, &len, fp )) != -1 )
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

