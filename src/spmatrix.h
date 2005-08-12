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
		int size;
		int *Ap;
		int *Ai;
		double *Ax;

	public:
		
		SpMatrix( ) : format('R'), n(0), size(0), Ap(0), Ai(0), Ax(0) { }
		SpMatrix( char F, int n_, int nz ) { Init( F, n_, nz ); }
		SpMatrix( const SpMatrix& M ) : BaseMatrix() { Init( M ); }
		virtual ~SpMatrix();
		
		void Init( char F, int n_, int nz );
		void Init( const SpMatrix& M );

		
		void AX( double* X, const double* B, double alpha, bool trans ) const;
		void AX( Vector& X, const Vector& B, double alpha, bool trans ) const 
		{ 
			BaseMatrix::AX( X, B, alpha, trans ); 
		}
	
		void AX( Matrix& X, const Matrix& B, double alpha, bool trans ) const;
		void AXpY( double* out, const double* X, const double* Y, double alpha, double beta, bool trans ) const;
		void AXpY( Vector& out, const Vector& X, const Vector& Y, double alpha, double beta, bool trans ) const 
		{ 
			BaseMatrix::AXpY( out, X, Y, alpha, beta, trans ); 
		}
	
		void AXpY( Matrix& out, const Matrix& X, const Matrix& Y, double alpha, double beta, bool trans ) const;
    
		int Size() const { return n; }
		int Row() const { return n; }
		int Col() const { return n; }
		
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
  
		// friend void GenEigval( SpFact& A, SpMatrix& B, Vector& wr, Vector& wi );

//		friend class SpFact;
};

// for 1 par continuation
// void BEMSolve(    BaseFact&   A,   Vector&   B, 
//                   Vector&     C,   double&   D,
//                   Vector&     F,   double&   G,
//                   Vector&     X,   double&   Y );
// void BEMWSolve(   BaseFact&   A,   Matrix&   B, 
//                   Matrix&     C,   Matrix&   D,
//                   Vector&     F,   Vector&   G,
//                   Vector&     X,   Vector&   Y );
// void GMBESolve(   BaseFact&   A11,                Vector& A13,
//                   BaseMatrix& A21, BaseFact& A22, Vector& A23,
//                   Vector&     A31, Vector&   A32, double& A33,
//                   Vector&     F1,  Vector&   F2,  double& F3,
//                   Vector&     X1,  Vector&   X2,  double& X3 );
// void GMBEWSolve(  BaseFact&   A11,                Matrix& A13,
//                   BaseMatrix& A21, BaseFact& A22, Matrix& A23,
//                   Matrix&     A31, Matrix&   A32, Matrix& A33,
//                   Vector&     F1,  Vector&   F2,  Vector& F3,
//                   Vector&     X1,  Vector&   X2,  Vector& X3 );

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
			for( int i=0; i<AI.Size(); i++ ) AI(i).Init( 'R', n_, nz );
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
