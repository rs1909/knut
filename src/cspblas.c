// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages packages root directory
//
// ------------------------------------------------------------------------- //

#include "cspblas.h"

static
inline double DnDot( const int n, double* in1, const double* in2, const double alpha )
{
	int j;
	register double out, r0, r1, r2, r3, r4, r5, r6, r7, t0, t1;
	
	out = 0.0;
	for( j = 0; j+7 < n; j+=8 )
	{
		r0 = in1[j] * in2[j];
		r1 = in1[j+1] * in2[j+1];
		r2 = in1[j+2] * in2[j+2];
		r3 = in1[j+3] * in2[j+3];
		r4 = in1[j+4] * in2[j+4];
		r5 = in1[j+5] * in2[j+5];
		r6 = in1[j+6] * in2[j+6];
		r7 = in1[j+7] * in2[j+7];
		t0 = r0 + r1 + r2 + r3;
		t1 = r4 + r5 + r6 + r7;
		out += alpha * ( t0 + t1 );
	}
	for( ; j < n; j+=1 )
	{
		out += alpha * in1[j] * in2[j];
	}
	return out;
}

static
inline void DnVCopy( const int n, double* out, const double* in, const double alpha )
{
	int j;
	for( j = 0; j+7 < n; j+=8 )
	{
		out[j]   = alpha * in[j];
		out[j+1] = alpha * in[j+1];
		out[j+2] = alpha * in[j+2];
		out[j+3] = alpha * in[j+3];
		out[j+4] = alpha * in[j+4];
		out[j+5] = alpha * in[j+5];
		out[j+6] = alpha * in[j+6];
		out[j+7] = alpha * in[j+7];
	}
	for( ; j < n; j+=1 )
	{
		out[j] = alpha * in[j];
	}
}

static
inline void SpMulR( const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                    double* out, const double* in, const double alpha, const int incIN )
{
	register int j, k;
	register double r0, r1, r2, r3, r4, r5, r6, r7, t0, t1;

	for( j = 0; j < n; j+=1 )
	{
		out[j] = 0.0;
		for( k = Ap[j]; k+7 < Ap[j+1]; k+=8 )
		{
			r0 = Ax[k]*in[incIN*Ai[k]];
			r1 = Ax[k+1]*in[incIN*Ai[k+1]];
			r2 = Ax[k+2]*in[incIN*Ai[k+2]];
			r3 = Ax[k+3]*in[incIN*Ai[k+3]];
			t0 = r0 + r1 + r2 + r3;
			r4 = Ax[k+4]*in[incIN*Ai[k+4]];
			r5 = Ax[k+5]*in[incIN*Ai[k+5]];
			r6 = Ax[k+6]*in[incIN*Ai[k+6]];
			r7 = Ax[k+7]*in[incIN*Ai[k+7]];
			t1 = r4 + r5 + r6 + r7;
			out[j] += t0 + t1;
		}
		for( ; k < Ap[j+1]; k++ )
		{
			out[j] += Ax[k]*in[incIN*Ai[k]];
		}
		out[j] *= alpha;
	}
}

static
inline void SpMulRpY( const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                      double* out, const double* in, const double alpha,
                      const double* Y, double beta, const int incIN, const int incY )
{
	register int j, k;
	register double r0, r1, r2, r3, r4, r5, r6, r7, t0, t1;

	for( j = 0; j < n; j+=1 )
	{
		out[j] = 0.0;
		for( k = Ap[j]; k+7 < Ap[j+1]; k+=8 )
		{
			r0 = Ax[k]*in[incIN*Ai[k]];
			r1 = Ax[k+1]*in[incIN*Ai[k+1]];
			r2 = Ax[k+2]*in[incIN*Ai[k+2]];
			r3 = Ax[k+3]*in[incIN*Ai[k+3]];
			t0 = r0 + r1 + r2 + r3;
			r4 = Ax[k+4]*in[incIN*Ai[k+4]];
			r5 = Ax[k+5]*in[incIN*Ai[k+5]];
			r6 = Ax[k+6]*in[incIN*Ai[k+6]];
			r7 = Ax[k+7]*in[incIN*Ai[k+7]];
			t1 = r4 + r5 + r6 + r7;
			out[j] += t0 + t1;
		}
		for( ; k < Ap[j+1]; k++ )
		{
			out[j] += Ax[k]*in[incIN*Ai[k]];
		}
		out[j] = out[j] * alpha + Y[incY*j] * beta;
	}
}

static
inline void SpMulC( const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                    double* out, const double* in, const double alpha, const int incIN )
{
	register int i, j, k;
	register double t0;
	
	for( i = 0; i+7 < m; i+=8  )
	{
		out[i] = 0.0;
		out[i+1] = 0.0;
		out[i+2] = 0.0;
		out[i+3] = 0.0;
		out[i+4] = 0.0;
		out[i+5] = 0.0;
		out[i+6] = 0.0;
		out[i+7] = 0.0;
	}
	for( ; i < m; i+=1  ) { out[i] = 0.0; }
	for( j = 0; j < n; j+=1 )
	{
		t0 = alpha*in[incIN*j];
		k = Ap[j];
		for( k = Ap[j]; k+7 < Ap[j+1]; k+=8 )
		{
			out[Ai[k]] += Ax[k] * t0;
			out[Ai[k+1]] += Ax[k+1] * t0;
			out[Ai[k+2]] += Ax[k+2] * t0;
			out[Ai[k+3]] += Ax[k+3] * t0;
			out[Ai[k+4]] += Ax[k+4] * t0;
			out[Ai[k+5]] += Ax[k+5] * t0;
			out[Ai[k+6]] += Ax[k+6] * t0;
			out[Ai[k+7]] += Ax[k+7] * t0;
		}
		for( ; k < Ap[j+1]; k+=1 )
		{
			out[Ai[k]] += Ax[k] * t0;
		}
	}
}

static
inline void SpMulCpY( const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                      double* out, const double* in, const double alpha,
                      const double* Y, const double beta, const int incIN, const int incY )
{
	register int i, j, k;
	register double t0;
	
	for( i = 0; i+7 < m; i+=8  )
	{
		out[i] = Y[incY*i] * beta;
		out[i+1] = Y[incY*(i+1)] * beta;
		out[i+2] = Y[incY*(i+2)] * beta;
		out[i+3] = Y[incY*(i+3)] * beta;
		out[i+4] = Y[incY*(i+4)] * beta;
		out[i+5] = Y[incY*(i+5)] * beta;
		out[i+6] = Y[incY*(i+6)] * beta;
		out[i+7] = Y[incY*(i+7)] * beta;
	}
	for( i = 0; i < m; i+=1  )
	{
		out[i] = Y[incY*i] * beta;
	}
	for( j = 0; j < n; j+=1 )
	{
		t0 = alpha*in[incIN*j];
		for( k = Ap[j]; k+7 < Ap[j+1]; k+=8 )
		{
			out[Ai[k]] += Ax[k] * t0;
			out[Ai[k+1]] += Ax[k+1] * t0;
			out[Ai[k+2]] += Ax[k+2] * t0;
			out[Ai[k+3]] += Ax[k+3] * t0;
			out[Ai[k+4]] += Ax[k+4] * t0;
			out[Ai[k+5]] += Ax[k+5] * t0;
			out[Ai[k+6]] += Ax[k+6] * t0;
			out[Ai[k+7]] += Ax[k+7] * t0;
		}
		for( ; k < Ap[j+1]; k+=1 )
		{
			out[Ai[k]] += Ax[k] * t0;
		}
	}
}

double cblas_dndot( const int n, double* in1, const double* in2, const double alpha )
{
	if( alpha == 1.0 ) return DnDot( n, in1, in2, 1.0 );
	else if( alpha == -1.0 ) return DnDot( n, in1, in2, -1.0 );
	else return DnDot( n, in1, in2, alpha );
}

void cblas_dnvcopy( const int n, double* out, const double* in, const double alpha )
{
	if( alpha == 1.0 ) DnVCopy( n, out, in, 1.0 );
	else if( alpha == -1.0 ) DnVCopy( n, out, in, -1.0 );
	else DnVCopy( n, out, in, alpha );
}

void cspblas_mmx( char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax, 
                  double* out, const double* in, double alpha )
{
	if( alpha == 1.0 )
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			SpMulC( n, m, Ap, Ai, Ax, out, in, 1.0, 1 );
		}else
		{
			SpMulR( n, m, Ap, Ai, Ax, out, in, 1.0, 1 );
		}
	}else
	if( alpha == -1.0 )
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			SpMulC( n, m, Ap, Ai, Ax, out, in, -1.0, 1 );
		}else
		{
			SpMulR( n, m, Ap, Ai, Ax, out, in, -1.0, 1 );
		}
	}else
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			SpMulC( n, m, Ap, Ai, Ax, out, in, alpha, 1 );
		}else
		{
			SpMulR( n, m, Ap, Ai, Ax, out, in, alpha, 1 );
		}
	}
}

void cspblas_mmxpy( char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                    double* out, const double* in, double alpha, const double* C, double beta )
{
	if( alpha == 1.0 )
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			SpMulCpY( n, m, Ap, Ai, Ax, out, in, 1.0, C, beta, 1, 1 );
		}else
		{
			SpMulRpY( n, m, Ap, Ai, Ax, out, in, 1.0, C, beta, 1, 1 );
		}
	}else
	if( alpha == -1.0 )
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			SpMulCpY( n, m, Ap, Ai, Ax, out, in, -1.0, C, beta, 1, 1 );
		}else
		{
			SpMulRpY( n, m, Ap, Ai, Ax, out, in, -1.0, C, beta, 1, 1 );
		}
	}else
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			SpMulCpY( n, m, Ap, Ai, Ax, out, in, alpha, C, beta, 1, 1 );
		}else
		{
			SpMulRpY( n, m, Ap, Ai, Ax, out, in, alpha, C, beta, 1, 1 );
		}
	}
}

void cspblas_mmxm( char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                   double* out, int ldout, const double* in, int ldin, double alpha, int nrhs )
{
	int c;
	if( alpha == 1.0 )
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulC( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulR( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, 1 );
			}
		}
	}else
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulC( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulR( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, 1 );
			}
		}
	}
}

void cspblas_mmxmpym( char format, enum cspblas_Trans trans, const int n, const int m, const int* Ap, const int* Ai, const double* Ax,
                      double* out, int ldout,
                      const double* in, int ldin, double alpha,
                      const double* Y, int ldY, double beta, int nrhs )
{
	int c;
	if( alpha == 1.0 )
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulCpY( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, &Y[c*ldY], beta, 1, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulRpY( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], 1.0, &Y[c*ldY], beta, 1, 1 );
			}
		}
	}else
	{
		if( ((format == 'C')&&(trans == NoTrans))||((format == 'R')&&(trans == Trans)) )
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulCpY( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, &Y[c*ldY], beta, 1, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				SpMulRpY( n, m, Ap, Ai, Ax, &out[c*ldout], &in[c*ldin], alpha, &Y[c*ldY], beta, 1, 1 );
			}
		}
	}
}

// -----------------------------------------
//
// DENSE MATRICES
//
// -----------------------------------------

static
inline void DnMulN( const int n, const int m, const double* A, const int lda,
                    double* out, const double* in, const double alpha, const int incIN )
{
	register int i, j;
	register double r0, r1, r2, r3, r4, r5, r6, r7;

	for( i = 0; i+3 < n; i+=4 )
	{
		out[i  ] = 0.0;
		out[i+1] = 0.0;
		out[i+2] = 0.0;
		out[i+3] = 0.0;
		for( j = 0; j+3 < m; j+=4 )
		{
			r0 = A[i   + j*lda] * in[j];
			r1 = A[i+1 + j*lda] * in[j];
			r2 = A[i+2 + j*lda] * in[j];
			r3 = A[i+3 + j*lda] * in[j];
			
			r4 = A[i   + (j+1)*lda] * in[j+1];
			r5 = A[i+1 + (j+1)*lda] * in[j+1];
			r6 = A[i+2 + (j+1)*lda] * in[j+1];
			r7 = A[i+3 + (j+1)*lda] * in[j+1];
			
			r0 += A[i   + (j+2)*lda] * in[j+2];
			r1 += A[i+1 + (j+2)*lda] * in[j+2];
			r2 += A[i+2 + (j+2)*lda] * in[j+2];
			r3 += A[i+3 + (j+2)*lda] * in[j+2];
			
			r4 += A[i   + (j+3)*lda] * in[j+3];
			r5 += A[i+1 + (j+3)*lda] * in[j+3];
			r6 += A[i+2 + (j+3)*lda] * in[j+3];
			r7 += A[i+3 + (j+3)*lda] * in[j+3];
			out[i  ] += alpha * (r0 + r4);
			out[i+1] += alpha * (r1 + r5);
			out[i+2] += alpha * (r2 + r6);
			out[i+3] += alpha * (r3 + r7);
		}
		for( ; j < m; j+=1 )
		{
			out[i  ] += alpha * A[i   + j*lda] * in[j];
			out[i+1] += alpha * A[i+1 + j*lda] * in[j];
			out[i+2] += alpha * A[i+2 + j*lda] * in[j];
			out[i+3] += alpha * A[i+3 + j*lda] * in[j];
		}
	}
	for( ; i < n; i+=1 )
	{
		out[i] = 0.0;
		for( j = 0; j+3 < m; j+=4 )
		{
			r0 = A[i + (j  )*lda] * in[j];
			r1 = A[i + (j+1)*lda] * in[j+1];
			r2 = A[i + (j+2)*lda] * in[j+2];
			r3 = A[i + (j+3)*lda] * in[j+3];
			out[i] += alpha * (r0 + r1 + r2 + r3);
		}
		for( ; j < m; j+=1 )
		{
			out[i] += alpha * A[i   + j*lda] * in[j];
		}
	}
}

static
inline void DnMulNpY( const int n, const int m, const double* A, const int lda,
                      double* out, const double* in, const double alpha,
                      const double* Y, double beta, const int incIN, const int incY )
{
	register int i, j;
	register double r0, r1, r2, r3, r4, r5, r6, r7;

	for( i = 0; i+3 < n; i+=4 )
	{
		out[i  ] = beta * Y[i];
		out[i+1] = beta * Y[i+1];
		out[i+2] = beta * Y[i+2];
		out[i+3] = beta * Y[i+3];
		for( j = 0; j+3 < m; j+=4 )
		{
			r0 = A[i   + j*lda] * in[j];
			r1 = A[i+1 + j*lda] * in[j];
			r2 = A[i+2 + j*lda] * in[j];
			r3 = A[i+3 + j*lda] * in[j];
			
			r4 = A[i   + (j+1)*lda] * in[j+1];
			r5 = A[i+1 + (j+1)*lda] * in[j+1];
			r6 = A[i+2 + (j+1)*lda] * in[j+1];
			r7 = A[i+3 + (j+1)*lda] * in[j+1];
			
			r0 += A[i   + (j+2)*lda] * in[j+2];
			r1 += A[i+1 + (j+2)*lda] * in[j+2];
			r2 += A[i+2 + (j+2)*lda] * in[j+2];
			r3 += A[i+3 + (j+2)*lda] * in[j+2];
			
			r4 += A[i   + (j+3)*lda] * in[j+3];
			r5 += A[i+1 + (j+3)*lda] * in[j+3];
			r6 += A[i+2 + (j+3)*lda] * in[j+3];
			r7 += A[i+3 + (j+3)*lda] * in[j+3];
			out[i  ] += alpha * (r0 + r4);
			out[i+1] += alpha * (r1 + r5);
			out[i+2] += alpha * (r2 + r6);
			out[i+3] += alpha * (r3 + r7);
		}
		for( ; j < m; j+=1 )
		{
			out[i  ] += alpha * A[i   + j*lda] * in[j];
			out[i+1] += alpha * A[i+1 + j*lda] * in[j];
			out[i+2] += alpha * A[i+2 + j*lda] * in[j];
			out[i+3] += alpha * A[i+3 + j*lda] * in[j];
		}
	}
	for( ; i < n; i+=1 )
	{
		out[i] = beta * Y[i];
		for( j = 0; j+3 < m; j+=4 )
		{
			r0 = A[i + (j  )*lda] * in[j];
			r1 = A[i + (j+1)*lda] * in[j+1];
			r2 = A[i + (j+2)*lda] * in[j+2];
			r3 = A[i + (j+3)*lda] * in[j+3];
			out[i] += alpha * (r0 + r1 + r2 + r3);
		}
		for( ; j < m; j+=1 )
		{
			out[i] += alpha * A[i   + j*lda] * in[j];
		}
	}
}

static
inline void DnMulT( const int n, const int m, const double* A, const int lda,
                    double* out, const double* in, const double alpha, const int incIN )
{
	register int i, j;
	register double r0, r1, r2, r3, r4, r5, r6, r7, c0, c1, c2, c3, c4;

	for( i = 0; i+3 < m; i+=4 )
	{
		out[i  ] = 0.0;
		out[i+1] = 0.0;
		out[i+2] = 0.0;
		out[i+3] = 0.0;
		for( j = 0; j+3 < n; j+=4 )
		{
			r0 = A[j + (i  )*lda] * in[j];
			r1 = A[j + (i+1)*lda] * in[j];
			r2 = A[j + (i+2)*lda] * in[j];
			r3 = A[j + (i+3)*lda] * in[j];
			
			r4 = A[j+1 + (i )*lda] * in[j+1];
			r5 = A[j+1 + (i+1)*lda] * in[j+1];
			r6 = A[j+1 + (i+2)*lda] * in[j+1];
			r7 = A[j+1 + (i+3)*lda] * in[j+1];
			
			r0 += A[j+2 + (i  )*lda] * in[j+2];
			r1 += A[j+2 + (i+1)*lda] * in[j+2];
			r2 += A[j+2 + (i+2)*lda] * in[j+2];
			r3 += A[j+2 + (i+3)*lda] * in[j+2];
			
			r4 += A[j+3 + (i  )*lda] * in[j+3];
			r5 += A[j+3 + (i+1)*lda] * in[j+3];
			r6 += A[j+3 + (i+2)*lda] * in[j+3];
			r7 += A[j+3 + (i+3)*lda] * in[j+3];
			out[i  ] += alpha * (r0 + r4);
			out[i+1] += alpha * (r1 + r5);
			out[i+2] += alpha * (r2 + r6);
			out[i+3] += alpha * (r3 + r7);
		}
		for( ; j < n; j+=1 )
		{
			out[i  ] += alpha * A[j + (i  )*lda] * in[j];
			out[i+1] += alpha * A[j + (i+1)*lda] * in[j];
			out[i+2] += alpha * A[j + (i+2)*lda] * in[j];
			out[i+3] += alpha * A[j + (i+3)*lda] * in[j];
		}
	}
	for( ; i < m; i+=1 )
	{
		out[i] = 0.0;
		for( j = 0; j+3 < n; j+=4 )
		{
			r0 = A[j   + i*lda] * in[j];
			r1 = A[j+1 + i*lda] * in[j+1];
			r2 = A[j+2 + i*lda] * in[j+2];
			r3 = A[j+3 + i*lda] * in[j+3];
			out[i] += alpha * (r0 + r1 + r2 + r3);
		}
		for( ; j < n; j+=1 )
		{
			out[i] += alpha * A[j   + i*lda] * in[j];
		}
	}
}

static
inline void DnMulTpY( const int n, const int m, const double* A, const int lda,
                      double* out, const double* in, const double alpha,
                      const double* Y, double beta, const int incIN, const int incY )
{
	register int i, j;
	register double r0, r1, r2, r3, r4, r5, r6, r7, c0, c1, c2, c3, c4;

	for( i = 0; i+3 < m; i+=4 )
	{
		out[i  ] = beta * Y[i];
		out[i+1] = beta * Y[i+1];
		out[i+2] = beta * Y[i+2];
		out[i+3] = beta * Y[i+3];
		for( j = 0; j+3 < n; j+=4 )
		{
			r0 = A[j + (i  )*lda] * in[j];
			r1 = A[j + (i+1)*lda] * in[j];
			r2 = A[j + (i+2)*lda] * in[j];
			r3 = A[j + (i+3)*lda] * in[j];
			
			r4 = A[j+1 + (i )*lda] * in[j+1];
			r5 = A[j+1 + (i+1)*lda] * in[j+1];
			r6 = A[j+1 + (i+2)*lda] * in[j+1];
			r7 = A[j+1 + (i+3)*lda] * in[j+1];
			
			r0 += A[j+2 + (i  )*lda] * in[j+2];
			r1 += A[j+2 + (i+1)*lda] * in[j+2];
			r2 += A[j+2 + (i+2)*lda] * in[j+2];
			r3 += A[j+2 + (i+3)*lda] * in[j+2];
			
			r4 += A[j+3 + (i  )*lda] * in[j+3];
			r5 += A[j+3 + (i+1)*lda] * in[j+3];
			r6 += A[j+3 + (i+2)*lda] * in[j+3];
			r7 += A[j+3 + (i+3)*lda] * in[j+3];
			out[i  ] += alpha * (r0 + r4);
			out[i+1] += alpha * (r1 + r5);
			out[i+2] += alpha * (r2 + r6);
			out[i+3] += alpha * (r3 + r7);
		}
		for( ; j < n; j+=1 )
		{
			out[i  ] += alpha * A[j + (i  )*lda] * in[j];
			out[i+1] += alpha * A[j + (i+1)*lda] * in[j];
			out[i+2] += alpha * A[j + (i+2)*lda] * in[j];
			out[i+3] += alpha * A[j + (i+3)*lda] * in[j];
		}
	}
	for( ; i < m; i+=1 )
	{
		out[i] = beta * Y[i];
		for( j = 0; j+3 < n; j+=4 )
		{
			r0 = A[j   + i*lda] * in[j];
			r1 = A[j+1 + i*lda] * in[j+1];
			r2 = A[j+2 + i*lda] * in[j+2];
			r3 = A[j+3 + i*lda] * in[j+3];
			out[i] += alpha * (r0 + r1 + r2 + r3);
		}
		for( ; j < n; j+=1 )
		{
			out[i] += alpha * A[j   + i*lda] * in[j];
		}
	}
}

///
/// The real routines
///

void cblas_mmx( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                double* out, const double* in, double alpha )
{
	if( alpha == 1.0 )
	{
		if( trans == Trans )
		{
			DnMulT( n, m, A, lda, out, in, 1.0, 1 );
		}else
		{
			DnMulN( n, m, A, lda, out, in, 1.0, 1 );;
		}
	}else
	if( alpha == -1.0 )
	{
		if( trans == Trans )
		{
			DnMulT( n, m, A, lda, out, in, -1.0, 1 );
		}else
		{
			DnMulN( n, m, A, lda, out, in, -1.0, 1 );;
		}
	}else
	{
		if( trans == Trans )
		{
			DnMulT( n, m, A, lda, out, in, alpha, 1 );
		}else
		{
			DnMulN( n, m, A, lda, out, in, alpha, 1 );;
		}
	}
}

void cblas_mmxpy( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                  double* out, const double* in, double alpha, const double* C, double beta )
{
	if( alpha == 1.0 )
	{
		if( trans == Trans )
		{
			DnMulTpY( n, m, A, lda, out, in, 1.0, C, beta, 1, 1 );
		}else
		{
			DnMulNpY( n, m, A, lda, out, in, 1.0, C, beta, 1, 1 );
		}
	}else
	if( alpha == -1.0 )
	{
		if( trans == Trans )
		{
			DnMulTpY( n, m, A, lda, out, in, -1.0, C, beta, 1, 1 );
		}else
		{
			DnMulNpY( n, m, A, lda, out, in, -1.0, C, beta, 1, 1 );
		}
	}else
	{
		if( trans == Trans )
		{
			DnMulTpY( n, m, A, lda, out, in, alpha, C, beta, 1, 1 );
		}else
		{
			DnMulNpY( n, m, A, lda, out, in, alpha, C, beta, 1, 1 );
		}
	}
}

void cblas_mmxm( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                 double* out, int ldout, const double* in, int ldin, double alpha, int nrhs )
{
	int c;
	if( alpha == 1.0 )
	{
		if( trans == Trans )
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulT( n, m, A, lda, &out[c*ldout], &in[c*ldin], 1.0, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulN( n, m, A, lda, &out[c*ldout], &in[c*ldin], 1.0, 1 );
			}
		}
	}else
	{
		if( trans == Trans )
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulT( n, m, A, lda, &out[c*ldout], &in[c*ldin], alpha, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulN( n, m, A, lda, &out[c*ldout], &in[c*ldin], alpha, 1 );
			}
		}
	}
}

void cblas_mmxmpym( enum cspblas_Trans trans, const int n, const int m, const double* A, const int lda,
                    double* out, int ldout,
                    const double* in, int ldin, double alpha,
                    const double* Y, int ldY, double beta, int nrhs )
{
	int c;
	if( alpha == 1.0 )
	{
		if( trans == Trans )
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulTpY( n, m, A, lda, &out[c*ldout], &in[c*ldin], 1.0, &Y[c*ldY], beta, 1, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulNpY( n, m, A, lda, &out[c*ldout], &in[c*ldin], 1.0, &Y[c*ldY], beta, 1, 1 );
			}
		}
	}else
	{
		if( trans == Trans )
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulTpY( n, m, A, lda, &out[c*ldout], &in[c*ldin], alpha, &Y[c*ldY], beta, 1, 1 );
			}
		}else
		{
			for( c = 0; c < nrhs; c++ )
			{
				DnMulNpY( n, m, A, lda, &out[c*ldout], &in[c*ldin], alpha, &Y[c*ldY], beta, 1, 1 );
			}
		}
	}
}
