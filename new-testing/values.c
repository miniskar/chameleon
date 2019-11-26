/**
 *
 * @file values.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 ***
 *
 * @brief Chameleon testing values toutine to read/print the parameters
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2019-07-18
 *
 */
#include "testings.h"
#include <coreblas.h>

/**
 * @brief Convert the input string to an integer
 * @param[in] str
 *    The input string
 * @return The integer read.
 */
val_t pread_int( const char *str )
{
    val_t val;
    val.ival = atoi( str );
    return val;
}

/**
 * @brief Convert the input string to a float
 * @param[in] str
 *    The input string
 * @return The float read.
 */
val_t pread_float( const char *str )
{
    val_t val;
    val.sval = strtof( str, NULL );
    return val;
}

/**
 * @brief Convert the input string to a double
 * @param[in] str
 *    The input string
 * @return The double read.
 */
val_t pread_double( const char *str )
{
    val_t val;
    val.dval = strtod( str, NULL );
    return val;
}

/**
 * @brief Convert the input string to a complex single precision
 * @param[in] str
 *    The input string
 * @return The complex single precision read.
 */
val_t pread_complex32( const char *str )
{
    float re, im;
    val_t val;
    int rc;

    rc = sscanf( str, "%e,%e", &re, &im );

    if ( rc == 2 ) {
        val.cval = re + I * im;
    }
    else if (rc == 1){
        val.cval = re;
    }
    else {
        val.cval = nan("NaN");
    }

    return val;
}

/**
 * @brief Convert the input string to a complex double precision
 * @param[in] str
 *    The input string
 * @return The complex double precision read.
 */
val_t pread_complex64( const char *str )
{
    double re, im;
    val_t val;
    int rc;

    rc = sscanf( str, "%le,%le", &re, &im );

    if ( rc == 2 ) {
        val.zval = re + I * im;
    }
    else if (rc == 1){
        val.zval = re;
    }
    else {
        val.zval = nan("NaN");
    }

    return val;
}

/**
 * @brief Convert the input string to a cham_trans_t
 * @param[in] str
 *    The input string
 * @return The cham_trans_t read.
 */
val_t pread_trans( const char *str )
{
    val_t val;
    val.trans = ChamNoTrans;

    if ( ( strcasecmp( "ChamConjTrans", str ) == 0 ) ||
         ( strcasecmp( "ConjTrans", str ) == 0 ) )
    {
        val.trans = ChamConjTrans;
    }
    else if ( ( strcasecmp( "ChamTrans", str ) == 0 ) ||
              ( strcasecmp( "Trans", str ) == 0 ) )
    {
        val.trans = ChamTrans;
    }
    else if ( ( strcasecmp( "ChamNoTrans", str ) == 0 ) ||
              ( strcasecmp( "NoTrans", str ) == 0 ) )
    {
        val.trans = ChamNoTrans;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamConjTrans) || (v == (ChamConjTrans-ChamNoTrans)) ) {
            val.trans = ChamConjTrans;
        }
        else if ( (v == ChamTrans) || (v == (ChamTrans-ChamNoTrans)) ) {
            val.trans = ChamTrans;
        }
        else {
            val.trans = ChamNoTrans;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_uplo_t
 * @param[in] str
 *    The input string
 * @return The cham_uplo_t read.
 */
val_t pread_uplo( const char *str )
{
    val_t val;
    val.uplo = ChamUpperLower;

    if ( ( strcasecmp( "ChamUpper", str ) == 0 ) ||
         ( strcasecmp( "Upper",     str ) == 0 ) )
    {
        val.uplo = ChamUpper;
    }
    else if ( ( strcasecmp( "ChamLower", str ) == 0 ) ||
              ( strcasecmp( "Lower",     str ) == 0 ) )
    {
        val.uplo = ChamLower;
    }
    else if ( ( strcasecmp( "ChamUpperLower", str ) == 0 ) ||
              ( strcasecmp( "UpperLower",     str ) == 0 ) ||
              ( strcasecmp( "General",        str ) == 0 ) )
    {
        val.uplo = ChamUpperLower;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamUpper) || (v == 0) ) {
            val.uplo = ChamUpper;
        }
        else if ( (v == ChamLower) || (v == (ChamLower-ChamUpper)) ) {
            val.uplo = ChamLower;
        }
        else {
            val.uplo = ChamUpperLower;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_diag_t
 * @param[in] str
 *    The input string
 * @return The cham_diag_t read.
 */
val_t pread_diag( const char *str )
{
    val_t val;
    val.diag = ChamNonUnit;

    if ( ( strcasecmp( "ChamNonUnit", str ) == 0 ) ||
         ( strcasecmp( "NonUnit",     str ) == 0 ) )
    {
        val.diag = ChamNonUnit;
    }
    else if ( ( strcasecmp( "ChamUnit", str ) == 0 ) ||
              ( strcasecmp( "Unit",     str ) == 0 ) )
    {
        val.diag = ChamUnit;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamUnit) || (v == (ChamUnit-ChamNonUnit)) ) {
            val.diag = ChamUnit;
        }
        else {
            val.diag = ChamNonUnit;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_side_t
 * @param[in] str
 *    The input string
 * @return The cham_side_t read.
 */
val_t pread_side( const char *str )
{
    val_t val;
    val.side = ChamLeft;

    if ( ( strcasecmp( "ChamLeft", str ) == 0 ) ||
         ( strcasecmp( "Left",     str ) == 0 ) )
    {
        val.side = ChamLeft;
    }
    else if ( ( strcasecmp( "ChamRight", str ) == 0 ) ||
              ( strcasecmp( "Right",     str ) == 0 ) )
    {
        val.side = ChamRight;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamRight) || (v == (ChamRight-ChamLeft)) ) {
            val.side = ChamRight;
        }
        else {
            val.side = ChamLeft;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_normtype_t
 * @param[in] str
 *    The input string
 * @return The cham_normtype_t read.
 */
val_t pread_norm( const char *str )
{
    val_t val;
    val.ntype = ChamOneNorm;

    if ( ( strcasecmp( "ChamOneNorm", str ) == 0 ) ||
         ( strcasecmp( "OneNorm",     str ) == 0 ) )
    {
        val.ntype = ChamOneNorm;
    }
    else if ( ( strcasecmp( "ChamFrobeniusNorm", str ) == 0 ) ||
              ( strcasecmp( "FrobeniusNorm",     str ) == 0 ) )
    {
        val.ntype = ChamFrobeniusNorm;
    }
    else if ( ( strcasecmp( "ChamInfNorm", str ) == 0 ) ||
              ( strcasecmp( "InfNorm",     str ) == 0 ) )
    {
        val.ntype = ChamInfNorm;
    }
    else if ( ( strcasecmp( "ChamMaxNorm", str ) == 0 ) ||
              ( strcasecmp( "MaxNorm",     str ) == 0 ) )
    {
        val.ntype = ChamMaxNorm;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamMaxNorm) || (v == (ChamMaxNorm-ChamOneNorm)) ) {
            val.ntype = ChamMaxNorm;
        }
        else if ( (v == ChamInfNorm) || (v == (ChamInfNorm-ChamOneNorm)) ) {
            val.ntype = ChamInfNorm;
        }
        else if ( (v == ChamFrobeniusNorm) || (v == (ChamFrobeniusNorm-ChamOneNorm)) ) {
            val.ntype = ChamFrobeniusNorm;
        }
        else {
            val.ntype = ChamOneNorm;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a string
 * @param[in] str
 *    The input string
 * @return The string read.
 */
val_t pread_string( const char *str )
{
    val_t val;
    int i, len = strlen( str );
    val.str = strdup( str );
    for( i=0; i<len; i++ ) {
        if ( (val.str[i] == '\n') ||
             (val.str[i] == ' ') )
        {
            val.str[i] = '\0';
            break;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to an integer
 * @param[in] str
 *    The input string
 * @return The integer read.
 */
char *sprint_int( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*d", nbchar, val.ival );
    }
    else {
        rc = sprintf( str_in, ";%d", val.ival );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a float
 * @param[in] str
 *    The input string
 * @return The float read.
 */
char *sprint_float( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*e", chameleon_max( 13, nbchar ), val.sval );
    }
    else {
        rc = sprintf( str_in, ";%e", val.sval );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a double
 * @param[in] str
 *    The input string
 * @return The double read.
 */
char *sprint_double( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*e", chameleon_max( 13, nbchar ), val.dval );
    }
    else {
        rc = sprintf( str_in, ";%e", val.dval );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a complex single precision
 * @param[in] str
 *    The input string
 * @return The complex single precision read.
 */
char *sprint_complex32( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " (%*e,%*e)", chameleon_max( 13, nbchar ), crealf(val.cval), chameleon_max( 13, nbchar ), cimagf(val.cval) );
    }
    else {
        rc = sprintf( str_in, ";%e,%e", crealf(val.cval), cimagf(val.cval) );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a complex double precision
 * @param[in] str
 *    The input string
 * @return The complex double precision read.
 */
char *sprint_complex64( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " (%*e,%*e)", chameleon_max( 13, nbchar ), creal(val.zval), chameleon_max( 13, nbchar ), cimag(val.zval) );
    }
    else {
        rc = sprintf( str_in, ";%e,%e", creal(val.zval), cimag(val.zval) );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_trans_t
 * @param[in] str
 *    The input string
 * @return The cham_trans_t read.
 */
char *sprint_trans( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.trans == ChamConjTrans) ? "ConjTrans" :
                      ((val.trans == ChamTrans) ? "Trans" : "NoTrans") );
    }
    else {
        rc = sprintf( str_in, ";%d", val.trans );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_uplo_t
 * @param[in] str
 *    The input string
 * @return The cham_uplo_t read.
 */
char *sprint_uplo( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.uplo == ChamUpper) ? "Upper" :
                      ((val.uplo == ChamLower) ? "Lower" : "General") );
    }
    else {
        rc = sprintf( str_in, ";%d", val.uplo );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_diag_t
 * @param[in] str
 *    The input string
 * @return The cham_diag_t read.
 */
char *sprint_diag( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.diag == ChamUnit) ? "Unit" : "NonUnit" );
    }
    else {
        rc = sprintf( str_in, ";%d", val.diag );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_side_t
 * @param[in] str
 *    The input string
 * @return The cham_side_t read.
 */
char *sprint_side( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.side == ChamLeft) ? "Left" : "Right" );
    }
    else {
        rc = sprintf( str_in, ";%d", val.side );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_normtype_t
 * @param[in] str
 *    The input string
 * @return The cham_normtype_t read.
 */
char *sprint_norm( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        char *name;
        switch( val.ntype ) {
        case ChamMaxNorm:
            name = "Max";
            break;
        case ChamOneNorm:
            name = "One";
            break;
        case ChamInfNorm:
            name = "Inf";
            break;
        case ChamFrobeniusNorm:
            name = "Frb";
            break;
        default:
            name = "ERR";
        }
        rc = sprintf( str_in, " %-*s", nbchar, name );
    }
    else {
        rc = sprintf( str_in, ";%d", val.ntype );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a string
 * @param[in] str
 *    The input string
 * @return The string read.
 */
char *sprint_check( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*s", nbchar, ( val.ival == 0 ) ? "SUCCESS" : "FAILED" );
    }
    else {
        rc = sprintf( str_in, ";%s", ( val.ival == 0 ) ? "SUCCESS" : "FAILED" );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a string
 * @param[in] str
 *    The input string
 * @return The string read.
 */
char *sprint_string( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar, val.str );
    }
    else {
        rc = sprintf( str_in, ";%s", val.str );
    }
    return str_in+rc;
}

/**
 * @brief Generate a random float
 */
float
testing_salea()
{
    float val;
    CORE_splrnt( 1, 1, &val, 1, 1, 1, 0, random() );
    return val;
}

/**
 * @brief Generate a random double
 */
double
testing_dalea()
{
    double val;
    CORE_dplrnt( 1, 1, &val, 1, 1, 1, 0, random() );
    return val;
}

/**
 * @brief Generate a random complex single precision
 */
CHAMELEON_Complex32_t
testing_calea()
{
    CHAMELEON_Complex32_t val;
    CORE_cplrnt( 1, 1, &val, 1, 1, 1, 0, random() );
    return val;
}

/**
 * @brief Generate a random complex double precision
 */
CHAMELEON_Complex64_t
testing_zalea()
{
    CHAMELEON_Complex64_t val;
    CORE_zplrnt( 1, 1, &val, 1, 1, 1, 0, random() );
    return val;
}
