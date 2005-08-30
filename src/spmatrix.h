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
#include "error.h"
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
		
		inline SpMatrix( ) : format('R'), n(0), m(0), size(0), Ap(0), Ai(0), Ax(0) { }
		inline SpMatrix( char F, int n_, int m_, int nz ) { Init( F, n_, m_, nz ); }
		inline SpMatrix( char F, int n_, int nz ) { Init( F, n_, n_, nz ); }
		inline SpMatrix( const SpMatrix& M ) { Init( M ); }
		inline virtual ~SpMatrix();
		
		inline void Init( char F, int n_, int m_, int nz );
		inline void Init( const SpMatrix& M );
		
		inline void mmx( enum cspblas_Trans trans, double* out, const double* in, double alpha ) const
		{
			cspblas_mmx( format, trans, n, m, Ap, Ai,Ax, out, in, alpha );
		}
		inline void mmxpy( enum cspblas_Trans trans, double* out, const double* in, double alpha, const double* Y, double beta ) const
		{
			cspblas_mmxpy( format, trans, n, m, Ap, Ai, Ax, out, in, alpha, Y, beta );
		}
		inline void mmxm( enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha, int nrhs ) const
		{
			cspblas_mmxm( format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, nrhs );
		}
		inline void mmxmpym( enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha,
		                      const double* Y, int ldY, double beta, int nrhs ) const
		{
			cspblas_mmxmpym( format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, Y, ldY, beta, nrhs );
		}
		
		inline __op_mul_vec<SpMatrix,Vector>     operator*( Vector& v )
				{ return __op_mul_vec<SpMatrix,Vector>( __scal_vec_trans<SpMatrix>( *this ), __scal_vec_trans<Vector>( v ) ); }
		inline __op_mul_vec<SpMatrix,Matrix>     operator*( Matrix& v )
				{ return __op_mul_vec<SpMatrix,Matrix>( __scal_vec_trans<SpMatrix>( *this ), __scal_vec_trans<Matrix>( v ) ); }
		inline __op_mul_vec_rng<SpMatrix,Vector> operator*( __scal_vec_trans_rng<Vector> v )
				{ return __op_mul_vec_rng<SpMatrix,Vector>( __scal_vec_trans_rng<SpMatrix>( *this ), v ); }
		inline __op_mul_vec_rng<SpMatrix,Matrix> operator*( __scal_vec_trans_rng<Matrix> v )
				{ return __op_mul_vec_rng<SpMatrix,Matrix>( __scal_vec_trans_rng<SpMatrix>( *this ), v ); }
		inline __scal_vec_trans<SpMatrix>        operator!( )
				{ return __scal_vec_trans<SpMatrix>( *this, 1.0, Trans ); }
		inline __scal_vec_trans_rng<SpMatrix>    operator[ ] ( rng r )
				{ return __scal_vec_trans_rng<SpMatrix>( *this, r ); }
		
		inline void AX( double* out, const double* in, double alpha, bool trans ) const;
		inline void AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const;
		
		inline int Row() const { if( format == 'R' ) return n; else return m; }
		inline int Col() const { if( format == 'R' ) return m; else return n; }
		
		/// clear the matrix
		/// these have to be virtual, because these might be called as SpFact
		virtual void Clear();
		virtual void Clear(char F);
		/// checks the structure of the matrix
		void Check();
		
		// Fill in routines
		/// Creates a new line in the matrix
		inline int NewL( int size_ );
		/// Writes or returns the row or column index 
		/// into line `l' and the `e'-th element
		inline int& WrLi( int l, int e );
		/// Writes into line `l' and the `e'-th element
		inline double& WrLx( int l, int e );
		/// returns the length of the n_ -th line in the matrix
		inline int GetL( int n_ );
		/// returns the nonzero elements in the matrix
		inline int GetNZ(){ return Ap[n]; }
		/// returns the number of lines, e.g. columns or rows in the matrix depending on format
		inline int GetN(){ return n; }
		// Computation routines
		/// transposes the matrix into the other format
		void Swap();
		// these are used for debugging only
		/// plots the structure of the matrix
		void StrPlot( GnuPlot& pl );
		/// prints out Ap
		void PrintAp(){ for(int i=0; i<n+1; i++) std::cout<<Ap[i]<<'\t'; std::cout<<'\n'; }
		/// prints the whole matrix onto the screen
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
		inline SpFact( SpFact& ) : SpMatrix( 'F', 1, 1 )
		{
			std::cout<<"SpFact::SpFact(SpFact&): not implemented\n";
		}
		~SpFact();

		inline void New() { fact = false; }
		
		inline void SetIter( int N ){ Control[UMFPACK_IRSTEP] = N; }
		void Clear( );
		void Clear( char F );
		
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
	if( (F != 'R')&&(F != 'C') ) std::cout<<"SpMatrix::CONSTRUCTOR: invalid format specification.\n";
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
		std::cout<<"WrLi bound "<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n"; PDError(1); return Ai[Ap[l]+e];
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
		std::cout<<"WrLx bound "<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n"; PDError(1); return Ax[Ap[l]+e];
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
		std::cout<<"SpMatrix::GetL: Error\n"; return -1;
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

inline Vector& Vector::operator=( const __op_mul_vec<SpMatrix,Vector> R )
{
	R.op.vec.mmx( R.op.tr, this->v, R.vecA.vec.v, R.op.alpha );
// 	std::cout<<"__op_mul_vec<Matrix,Vector>\n";
	return *this;
}

inline Vector& Vector::operator=( const __op_mul_vec_plus_vec<SpMatrix,Vector> R )
{
	R.op.vec.mmxpy( R.op.tr, this->v, R.vecA.vec.v, R.op.alpha, R.vecB.vec.v, R.vecB.alpha );
// 	std::cout<<"__op_mul_vec_plus_vec<Matrix,Vector>\n";
	return *this;
}

inline Matrix& Matrix::operator=( const __op_mul_vec<SpMatrix,Matrix> R )
{
	R.op.vec.mmxm( R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecA.vec.c );
// 	std::cout<<"__op_mul_vec<Matrix,Matrix>\n";
	return *this;
}

inline Matrix& Matrix::operator=( const __op_mul_vec_plus_vec<SpMatrix,Matrix> R )
{
	R.op.vec.mmxmpym( R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecB.vec.m, R.vecB.vec.r, R.vecB.alpha, R.vecB.vec.c );
// 	std::cout<<"__op_mul_vec_plus_vec<Matrix,Matrix>\n";
	return *this;
}

#endif
