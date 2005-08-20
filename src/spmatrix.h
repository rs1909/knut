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

class SpMatrix : virtual public BaseMatrix
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
		SpMatrix( const SpMatrix& M ) : BaseMatrix() { Init( M ); }
		virtual ~SpMatrix();
		
		void Init( char F, int n_, int m_, int nz );
		void Init( const SpMatrix& M );

		void AX( double* out, const double* in, double alpha, bool trans ) const
		{
			cspblas_mmx( format, trans ? Trans : NoTrans, n, m, Ap, Ai, Ax, out, in, alpha );
		}
		
		void AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const
		{
			cspblas_mmxpy( format, trans ? Trans : NoTrans, n, m, Ap, Ai, Ax, out, in, alpha, Y, beta );
		}
		
		void AX( Vector& X, const Vector& B, double alpha, bool trans ) const
		{
			BaseMatrix::AX( X, B, alpha, trans );
		}
	
		void AX( Matrix& out, const Matrix& in, double alpha, bool trans ) const
		{
			cspblas_mmxm( format, trans ? Trans : NoTrans, n, m, Ap, Ai, Ax, out.m, out.r, in.m, in.r, alpha, in.c );
		}
		
		void AXpY( Matrix& out, const Matrix& in, const Matrix& Y, double alpha, double beta, bool trans ) const
		{
			cspblas_mmxmpym( format, trans ? Trans : NoTrans, n, m, Ap, Ai, Ax, out.m, out.r, in.m, in.r, alpha, Y.m, Y.r, beta, in.c );
		}
		
		void AXpY( Vector& out, const Vector& X, const Vector& Y, double alpha, double beta, bool trans ) const 
		{
			BaseMatrix::AXpY( out, X, Y, alpha, beta, trans );
		}
	
		
		// int Size() const { if( format=='Creturn n; }
		int Row() const { if( format == 'R' ) return n; else return m; }
		int Col() const { if( format == 'R' ) return m; else return n; }
		
		virtual void Clear();
		virtual void Clear(char F);
		
		void Check();
		
		// Fill in routines
		void AddC( int n_, int *idx, double *vec );
		void AddR( int n_, int *idx, double *vec );
		
		int NewL( int size_ )
		{
			n++; Ap[n] = Ap[n-1] + size_;
			return n-1;
		}
		
		int& WrLi( int l, int e )
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
		
		double& WrLx( int l, int e )
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
		
		int GetL( int n_ )
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

class SpFact : public SpMatrix, public BaseFact
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
		SpFact( SpFact& ) : BaseMatrix(), SpMatrix( 'F', 1, 1 ), BaseFact()
		{
			cout<<"SpFact::SpFact(SpFact&): not implemented\n";
		}
		~SpFact();

		void New() { fact = false; }  
		int  Row() const { return SpMatrix::Row(); }
		int  Col() const { return SpMatrix::Col(); }
		void AX( double* X, const double* B, double alpha, bool trans ) const 
		{
			SpMatrix::AX( X, B, alpha, trans );
		}
	
		void AX( Matrix& X, const Matrix& B, double alpha, bool trans ) const 
		{ 
			SpMatrix::AX( X, B, alpha, trans );
		}
	
		void AXpY( double* out, const double* X, double* Y, double alpha, double beta, bool trans ) const 
		{ 
			SpMatrix::AXpY( out, X, Y, alpha, beta, trans );
		}
	
		void AXpY( Matrix& out, const Matrix& X, Matrix& Y, double alpha, double beta, bool trans ) const 
		{ 
			SpMatrix::AXpY( out, X, Y, alpha, beta, trans );
		}
		
		void SetIter( int N ){ Control[UMFPACK_IRSTEP] = N; }
		void Clear( );
		void Clear( char F );
		int  Size() const { return n; }
		void Scale( Vector& x, const Vector& b );
		
		void Solve( double* x, double* b, bool trans = false );
		void Solve( Vector& x, const Vector& b, bool trans = false );
		void Solve( Matrix& x, const Matrix& b, bool trans = false );

		// friend void GenEigval( SpFact& A, SpMatrix& B, Vector& wr, Vector& wi );
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

// void GenEigval( SpFact& A, SpMatrix& B, Vector& wr, Vector& wi );

#endif
