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

#include <stdlib.h>

enum cspblas_Trans{ NoTrans = 0, Trans = 1 };

/* SPARSE MATRIX ROUTINES*/

void cspblas_mmx(char format, enum cspblas_Trans, const size_t n, const size_t m, const int* Ap, const int* Ai, const double* Ax,
                 double* out, const double* in, double alpha);

void cspblas_mmxpy(char format, enum cspblas_Trans, const size_t n, const size_t m, const int* Ap, const int* Ai, const double* Ax,
                   double* out, const double* in, double alpha, const double* C, double beta);

void cspblas_mmxm(char format, enum cspblas_Trans, const size_t n, const size_t m, const int* Ap, const int* Ai, const double* Ax,
                  double* out, size_t ldout, const double* in, size_t ldin, double alpha, size_t nrhs);

void cspblas_mmxmpym(char format, enum cspblas_Trans, const size_t n, const size_t m, const int* Ap, const int* Ai, const double* Ax,
                     double* out, size_t ldout,
                     const double* in, size_t ldin, double alpha,
                     const double* Y, size_t ldY, double beta, size_t nrhs);

#endif
