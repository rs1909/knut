// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages packages root directory
//
// ------------------------------------------------------------------------- //

#ifndef CSPBLAS_H
#define CSPBLAS_H

enum cspblas_Trans{ NoTrans = 0, Trans = 1 }; 

// SPARSE MATRIX ROUTINES

void cspblas_mmx( char format, enum cspblas_Trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax, 
                  double* out, const double* in, double alpha );

void cspblas_mmxpy( char format, enum cspblas_Trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                    double* out, const double* in, double alpha, const double* C, double beta );

void cspblas_mmxm( char format, enum cspblas_Trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                   double* out, int ldout, const double* in, int ldin, double alpha, int nrhs );

void cspblas_mmxmpym( char format, enum cspblas_Trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                      double* out, int ldout,
                      const double* in, int ldin, double alpha,
                      const double* Y, int ldY, double beta, int nrhs );

// DENSE MATRIX ROUTINES

void cblas_mmx( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                double* out, const double* in, double alpha );

void cblas_mmxpy( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                  double* out, const double* in, double alpha, const double* C, double beta );

void cblas_mmxm( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                 double* out, int ldout, const double* in, int ldin, double alpha, int nrhs );

void cblas_mmxmpym( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                    double* out, int ldout,
                    const double* in, int ldin, double alpha,
                    const double* Y, int ldY, double beta, int nrhs );

#endif
