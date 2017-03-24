/* ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- */

#ifndef CSPBLAS_H
#define CSPBLAS_H

#include <cstdlib>
#include <cstddef>

enum cspblas_Trans{ NoTrans = 0, Trans = 1 };

/* SPARSE MATRIX ROUTINES*/

void cspblas_mmx(char format, enum cspblas_Trans, const ptrdiff_t n, const ptrdiff_t m, const int* Ap, const int* Ai, const double* Ax,
                 double* out, const double* in, double alpha);

void cspblas_mmxpy(char format, enum cspblas_Trans, const ptrdiff_t n, const ptrdiff_t m, const int* Ap, const int* Ai, const double* Ax,
                   double* out, const double* in, double alpha, const double* C, double beta);

void cspblas_mmxm(char format, enum cspblas_Trans, const ptrdiff_t n, const ptrdiff_t m, const int* Ap, const int* Ai, const double* Ax,
                  double* out, ptrdiff_t ldout, const double* in, ptrdiff_t ldin, double alpha, ptrdiff_t nrhs);

void cspblas_mmxmpym(char format, enum cspblas_Trans, const ptrdiff_t n, const ptrdiff_t m, const int* Ap, const int* Ai, const double* Ax,
                     double* out, ptrdiff_t ldout,
                     const double* in, ptrdiff_t ldin, double alpha,
                     const double* Y, ptrdiff_t ldY, double beta, ptrdiff_t nrhs);

#endif
