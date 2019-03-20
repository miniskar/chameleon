/**
 *
 * @file global.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global coreblas variables and functions
 *
 * @version 0.9.2
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @date 2014-11-16
 *
 */
static int coreblas_gemm3m_enabled = 0;

void
set_coreblas_gemm3m_enabled( int v ) {
    coreblas_gemm3m_enabled = v;
}

int
get_coreblas_gemm3m_enabled(void) {
    return coreblas_gemm3m_enabled;
}

/**
 *  LAPACK Constants
 */
char *chameleon_lapack_constants[] =
{
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 100
    "Row",                  // 101: ChamRowMajor
    "Column",               // 102: ChamColMajor
    "", "", "", "", "", "", "", "",
    "No transpose",         // 111: ChamNoTrans
    "Transpose",            // 112: ChamTrans
    "Conjugate transpose",  // 113: ChamConjTrans
    "", "", "", "", "", "", "",
    "Upper",                // 121: ChamUpper
    "Lower",                // 122: ChamLower
    "All",                  // 123: ChamUpperLower
    "", "", "", "", "", "", "",
    "Non-unit",             // 131: ChamNonUnit
    "Unit",                 // 132: ChamUnit
    "", "", "", "", "", "", "", "",
    "Left",                 // 141: ChamLeft
    "Right",                // 142: ChamRight
    "", "", "", "", "", "", "", "",
    "",                     // 151:
    "",                     // 152:
    "",                     // 153:
    "",                     // 154:
    "",                     // 155:
    "",                     // 156:
    "Epsilon",              // 157: ChamEps
    "",                     // 158:
    "",                     // 159:
    "",                     // 160:
    "", "", "", "", "", "", "", "", "", "",
    "One norm",             // 171: ChamOneNorm
    "",                     // 172: ChamRealOneNorm
    "",                     // 173: ChamTwoNorm
    "Frobenius norm",       // 174: ChamFrobeniusNorm
    "Infinity norm",        // 175: ChamInfNorm
    "",                     // 176: ChamRealInfNorm
    "Maximum norm",         // 177: ChamMaxNorm
    "",                     // 178: ChamRealMaxNorm
    "",                     // 179
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 200
    "Uniform",              // 201: ChamDistUniform
    "Symmetric",            // 202: ChamDistSymmetric
    "Normal",               // 203: ChamDistNormal
    "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 240
    "Hermitian",            // 241 ChamHermGeev
    "Positive ev Hermitian",// 242 ChamHermPoev
    "NonSymmetric pos sv",  // 243 ChamNonsymPosv
    "Symmetric pos sv",     // 244 ChamSymPosv
    "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 290
    "No Packing",           // 291 ChamNoPacking
    "U zero out subdiag",   // 292 ChamPackSubdiag
    "L zero out superdiag", // 293 ChamPackSupdiag
    "C",                    // 294 ChamPackColumn
    "R",                    // 295 ChamPackRow
    "B",                    // 296 ChamPackLowerBand
    "Q",                    // 297 ChamPackUpeprBand
    "Z",                    // 298 ChamPackAll
    "",                     // 299

    "",                     // 300
    "No vectors",           // 301 ChamNoVec
    "Vectors needed",       // 302 ChamVec
    "I",                    // 303 ChamIvec
    "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 390
    "Forward",              // 391
    "Backward",             // 392
    "", "", "", "", "", "", "", "",
    "Columnwise",           // 401
    "Rowwise",              // 402
    "", "", "", "", "", "", "", ""  // Remember to add a coma!
};
