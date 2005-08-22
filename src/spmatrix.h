// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef SPMATRIX_H
#define SPMATRIX_H

#include "matrix.h"
#include "plot.h"
#include <vector>

extern "C" {

#include "cspblas.h"
#include "umfpack.h"

}

// **************************************************************************//
//                                                                           //
// **********************        SpMatrix Class        **********************//
//                                                                           //
// **************************************************************************//

class SpMatrix
{
	
	protected:
	
		char format;
		int n;
		int m;
		int size;
		int *Ap;
		int *Ai;
		double *Ax;

	public:
		
		SpMatrix( ) : format('R'), n(0), m(0), size(0), Ap(0), Ai(0), Ax(0) { }
		SpMatrix( char F, int n_, int m_, int nz ) { Init( F, n_, m_, nz ); }
		SpMatrix( char F, int n_, int nz ) { Init( F, n_, n_, nz ); }
		SpMatrix( const SpMatrix& M ) /*: BaseMatrix()*/ { Init( M ); }
		virtual ~SpMatrix();
		
		void Init( char F, int n_, int m_, int nz );
		void Init( const SpMatrix& M );
		
		virtual void mmx( enum cspblas_Trans trans, double* out, const double* in, double alpha ) const
		{
			cspblas_mmx( format, trans, n, m, Ap, Ai,Ax, out, in, alpha );
		}
		virtual void mmxpy( enum cspblas_Trans trans, double* out, const double* in, double alpha, const double* Y, double beta ) const
		{
			cspblas_mmxpy( format, trans, n, m, Ap, Ai, Ax, out, in, alpha, Y, beta );
		}
		virtual void mmxm( enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha, int nrhs ) const
		{
			cspblas_mmxm( format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, nrhs );
		}
		virtual void mmxmpym( enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha,
		                      const double* Y, int ldY, double beta, int nrhs ) const
		{
			cspblas_mmxmpym( format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, Y, ldY, beta, nrhs );
		}
		
		__op_mul_vec<SpMatrix,Vector> operator*( Vector& v ) { return __op_mul_vec<SpMatrix,Vector>( *this, v ); }
		__op_mul_vec<SpMatrix,Matrix> operator*( Matrix& v ) { return __op_mul_vec<SpMatrix,Matrix>( *this, v ); }
		__vec_trans<SpMatrix>         operator!( ) { return __vec_trans<SpMatrix>( *this ); }
		
		void AX( double* out, const double* in, double alpha, bool trans ) const;
		void AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const;
		
		int Row() const;
		int Col() const;
		
		virtual void Clear();
		virtual void Clear(char F);
		
		void Check();
		
		// Fill in routines
		void AddC( int n_, int *idx, double *vec );
		void AddR( int n_, int *idx, double *vec );
		
		int NewL( int size_ );
		
		int& WrLi( int l, int e );
		
		double& WrLx( int l, int e );
		
		int GetL( int n_ );
		
		int GetNZ(){ return Ap[n]; }
		int GetN(){ return n; }
		// Computation routines
		void Swap();
		// void Eigval( Vector& re, Vector& im ); // not used anymore
		void StrPlot( GnuPlot& pl );
		void PrintAp(){ for(int i=0; i<n+1; i++) cout<<Ap[i]<<'\t'; cout<<'\n'; }
		void Print();
};

// **************************************************************************//
//                                                                           //
// **********************         SpFact Class         **********************//
//                                                                           //
// **************************************************************************//

class SpFact : public SpMatrix
{
	private:
	
		bool   fact; 
		void   *Numeric; 
		double Control[UMFPACK_CONTROL];
		double Info[UMFPACK_INFO];
		// temporary space
		int    *Wi;
		double *W;
		
	public:
	
		SpFact( char F, int n_, int m_, int nz );
		SpFact( char F, int n_, int nz );
		SpFact( SpMatrix& M );
		SpFact( SpFact& ) : SpMatrix( 'F', 1, 1 )
		{
			cout<<"SpFact::SpFact(SpFact&): not implemented\n";
		}
		~SpFact();

		void New() { fact = false; }
		
		void SetIter( int N ){ Control[UMFPACK_IRSTEP] = N; }
		void Clear( );
		void Clear( char F );
		int  Size() const { return n; }
		void Scale( Vector& x, const Vector& b );
		
		void Solve( double* x, double* b, bool trans = false );
		void Solve( Vector& x, const Vector& b, bool trans = false );
		void Solve( Matrix& x, const Matrix& b, bool trans = false );
		
	private:
		void Fact();
};

// void NColloc::StabJac( StabMat& AB, const Vector& par, const JagMatrix3D& solData );
// this computes its own matrix structures, because it is nowhere else needed: kkSI, eeSI, rrSI, ddSI, etc.
// However the variables will be contained within the NColloc class

class StabMatrix
{
	public:
	
		StabMatrix( int nmat_, int n_, int nz ) : A0( 'R', n_,  nz ), AI(nmat_)
		{
			for( int i=0; i<AI.Size(); i++ ) AI(i).Init( 'R', n_, n_, nz );
			RESID = new double[ nmat_ * n_ + 1 ];
			isINIT = false;
		}
		~StabMatrix() { delete[] RESID; }
		
		int                nmat() { return AI.Size(); }
		SpFact&            getA0() { return A0; }
		Array1D<SpMatrix>& getAI() { return AI; }
		SpMatrix&          getAI( int i ) { return AI(i); }
		
		void Eigval( Vector& wr, Vector& wi );
	
	private:
	
		SpFact             A0;
		Array1D<SpMatrix>  AI;
		
		double* RESID;
		bool    isINIT;
};

// Implementation of SpMatrix

inline void SpMatrix::Init( char F, int n_, int m_, int nz )
{
	if( (F != 'R')&&(F != 'C') ) cout<<"SpMatrix::CONSTRUCTOR: invalid format specification.\n";
	format = F;
	n = 0;
	m = m_;
	size = nz;
	Ap = new int[n_+1];
	Ai = new int[nz];
	Ax = new double[nz];
	for( int i = 0; i < nz; i++ ){ Ai[i] = 0; Ax[i] = 0.0; }
	for( int i = 0; i < n_+1; i++ ){ Ap[i] = 0; }
}

inline void SpMatrix::Init( const SpMatrix& M )
{
	format = M.format;
	n = M.n;
	m = M.m;
	size = M.size;
	
	Ap = new int[M.n+1];
	Ai = new int[size];
	Ax = new double[size];
	for( int i = 0; i < size; i++ ){ Ai[i] = M.Ai[i]; Ax[i] = M.Ax[i]; }
	for( int i = 0; i < n+1; i++ ){ Ap[i] = M.Ap[i]; }
}

inline SpMatrix::~SpMatrix()
{
	delete []Ap;
	delete []Ai;
	delete []Ax;
}

inline void SpMatrix::Clear()
{
	
	for( int i = 0; i < n+1; i++ ){ Ap[i] = 0; }
	n = 0;
	for( int i = 0; i < size; i++ ){ Ai[i] = 0; Ax[i] = 0.0; }
}

inline void SpMatrix::Clear(char F)
{
	
	for( int i = 0; i < n+1; i++ ){ Ap[i] = 0; }
	n = 0;
	format = F;
	for( int i = 0; i < size; i++ ){ Ai[i] = 0; Ax[i] = 0.0; }
}

inline int SpMatrix::Row() const { if( format == 'R' ) return n; else return m; }
inline int SpMatrix::Col() const { if( format == 'R' ) return m; else return n; }

inline int SpMatrix::NewL( int size_ )
{
	n++; Ap[n] = Ap[n-1] + size_;
	return n-1;
}

inline int& SpMatrix::WrLi( int l, int e )
{
#ifdef DEBUG
	if( (l < n)&&(Ap[l] + e < Ap[l+1])&&(e >= 0)&&(l >= 0)&&(Ap[l]+e < size) )
	{
#endif
		return Ai[Ap[l]+e];
#ifdef DEBUG
	}else
	{
		cout<<"WrLi bound "<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n"; throw(1); return Ai[Ap[l]+e];
	}
#endif
}

inline double& SpMatrix::WrLx( int l, int e )
{
#ifdef DEBUG
	if( (l < n)&&(Ap[l] + e < Ap[l+1])&&(e >= 0)&&(l >= 0)&&(Ap[l]+e < size) )
	{
#endif
		return Ax[Ap[l]+e];
#ifdef DEBUG
	}else
	{
		cout<<"WrLx bound "<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n"; throw(1); return Ax[Ap[l]+e];
	}
#endif
}

inline int SpMatrix::GetL( int n_ )
{
#ifdef DEBUG
	if( n_ < n )
	{
#endif
		return Ap[n_+1] - Ap[n_];
#ifdef DEBUG
	}else
	{
		cout<<"SpMatrix::GetL: Error\n"; return -1;
	}
#endif
}

inline void SpMatrix::AX( double* out, const double* in, double alpha, bool trans ) const
{
	mmx( trans ? Trans : NoTrans, out, in, alpha );
}

inline void SpMatrix::AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const
{
	mmxpy( trans ? Trans : NoTrans, out, in, alpha, Y, beta );
}

// End of implementation of SpMatrix


// Implementation of Vector

inline Vector& Vector::operator=( const __op_mul_vec<SpMatrix,Vector> op )
{
	op.op.mmx( NoTrans, this->v, op.vecA.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __op_trans_mul_vec<SpMatrix,Vector> op )
{
	op.op.mmx( Trans, this->v, op.vecA.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_mul_vec<SpMatrix,Vector> op )
{
	op.op.mmx( NoTrans, this->v, op.vecA.v, op.alpha );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_trans_mul_vec<SpMatrix,Vector> op )
{
	op.op.mmx( Trans, this->v, op.vecA.v, op.alpha );
	return *this;
}
inline Vector& Vector::operator=( const __op_mul_vec_plus_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, 1.0, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __op_trans_mul_vec_plus_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, 1.0, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_mul_vec_plus_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, op.alpha, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_trans_mul_vec_plus_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, op.alpha, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __op_mul_vec_plus_scal_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, 1.0, op.vecB.v, op.beta );
	return *this;
}
inline Vector& Vector::operator=( const __op_trans_mul_vec_plus_scal_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, 1.0, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_mul_vec_plus_scal_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, op.alpha, op.vecB.v, op.beta );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_trans_mul_vec_plus_scal_vec<SpMatrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, op.alpha, op.vecB.v, op.beta );
	return *this;
}

// End of implementation of Vector

// Implementation of Matrix

inline Matrix& Matrix::operator=( const __op_mul_vec<SpMatrix,Matrix> op )
{
	op.op.mmxm( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_trans_mul_vec<SpMatrix,Matrix> op )
{
	op.op.mmxm( Trans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_mul_vec<SpMatrix,Matrix> op )
{
	op.op.mmxm( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_trans_mul_vec<SpMatrix,Matrix> op )
{
	op.op.mmxm( Trans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_mul_vec_plus_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_trans_mul_vec_plus_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_mul_vec_plus_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_trans_mul_vec_plus_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_mul_vec_plus_scal_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, op.beta, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_trans_mul_vec_plus_scal_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_mul_vec_plus_scal_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, op.beta, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_trans_mul_vec_plus_scal_vec<SpMatrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, op.beta, op.vecA.c );
	return *this;
}

#endif
