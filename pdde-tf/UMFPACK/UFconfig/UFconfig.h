/* ========================================================================== */
/* === UFconfig.h =========================================================== */
/* ========================================================================== */

/* Configuration file for most UFsparse packages (AMD, COLAMD, CCOLAMD,
 * CAMD, CHOLMOD, UMFPACK, and CXSparse).  KLU and BTF do not include a long
 * version.  LDL and CSparse do not, and never will (use CXSparse instead).
 *
 * UFconfig.h provides the definition of the long integer.  On most systems,
 * a C program can be compiled in LP64 mode, in which long's and pointers are
 * both 64-bits, and int's are 32-bits.  Windows 64, however, uses the LLP64
 * model, in which int's and long's are 32-bits, and long long's and pointers
 * are 64-bits.
 *
 * The UFsparse matrix packages that include long integer versions are
 * intended for the LP64 mode.  However, as a workaround for Windows 64
 * (and perhaps other systems), the long integer can be redefined.
 *
 * If _WIN64 is defined, then the __int64 type is used instead of long.
 * The long integer can also be defined at compile time.  For example, this
 * could be added to UFconfig.mk:
 *
 * CFLAGS = -O -D'UF_long=long long' -D'UF_long_max=9223372036854775801' \
 *   -D'UF_long_id="%lld"'
 *
 * This file defines UF_long as either long (on all but _WIN64) or
 * __int64 on Windows 64.  The intent is that a UF_long is always a 64-bit
 * integer in a 64-bit code.  ptrdiff_t might be a better choice than long;
 * it is always the same size as a pointer.
 *
 * This file also defines the UFSPARSE_VERSION and related definitions.
 *
 * Copyright (c) 2006, University of Florida.  No licensing restrictions
 * apply to this file or to the UFconfig directory.  Author: Timothy A. Davis.
 */

#ifndef _UFCONFIG_H
#define _UFCONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <limits.h>

/* ========================================================================== */
/* === UF_long ============================================================== */
/* ========================================================================== */

#ifndef UF_long

#ifdef _WIN64

#define UF_long __int64
#define UF_long_max _I64_MAX
#define UF_long_id "Id"

#else

#define UF_long long
#define UF_long_max LONG_MAX
#define UF_long_id "%ld"

#endif
#endif

/* ========================================================================== */
/* === UFsparse version ===================================================== */
/* ========================================================================== */

/* UFsparse is not a package itself, but a collection of packages, some of
 * which must be used together (UMFPACK requires AMD, CHOLMOD requires AMD,
 * COLAMD, CAMD, and CCOLAMD, etc).  A version number is provided here
 * for the collection itself.  The versions of packages within each version
 * of UFsparse are meant to work together.  Combining one packge from one
 * version of UFsparse, with another package from another version of UFsparse,
 * may or may not work.
 *
 * UFsparse Version 2.1.0 contains the following packages:
 *
 *  AMD		version 2.0.1
 *  CAMD	version 2.1.1
 *  COLAMD	version 2.5.1
 *  CCOLAMD	version 2.5.1
 *  CHOLMOD	version 1.2
 *  CSparse	version 2.0.2
 *  CXSparse	version 2.0.2
 *  KLU		version 0.9
 *  BTF		version 0.9
 *  LDL		version 1.3
 *  UFconfig	version number is the same as UFsparse
 *  UMFPACK	version 5.0.1
 *
 * Other package dependencies:
 *  BLAS	required by CHOLMOD and UMFPACK
 *  LAPACK	required by CHOLMOD
 *  METIS 4.0.1	required by CHOLMOD (optional)
 */

#define UFSPARSE_DATE "Aug 31, 2006"
#define UFSPARSE_VER_CODE(main,sub) ((main) * 1000 + (sub))
#define UFSPARSE_MAIN_VERSION 2
#define UFSPARSE_SUB_VERSION 1
#define UFSPARSE_SUBSUB_VERSION 0
#define UFSPARSE_VERSION \
    UFSPARSE_VER_CODE(UFSPARSE_MAIN_VERSION,UFSPARSE_SUB_VERSION)

#ifdef __cplusplus
}
#endif


#endif
