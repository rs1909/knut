// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef HYPERMATRIX_H
#define HYPERMATRIX_H

#include "matrix.h"
#include "spmatrix.h"

class HyperVector
{
	public:
	
		HyperVector( int i, int j, int k ) : V1(i), V2(j), V3(k) { }
		~HyperVector() { }
		
		Vector& getV1() { return V1; }
		Vector& getV2() { return V2; }
		Vector& getV3() { return V3; }
		
	private:
	
		Vector V1;
		Vector V2;
		Vector V3;
};

class HyperMatrix
{
	public:
	
		// constructing Sparse-Dense type hypermatrix
		HyperMatrix( int i, int j, int k, int nz );
		~HyperMatrix();
		
		SpFact&      getA11() { return *A11; }
		// MatFact&     getA11D() { return static_cast<MatFact&>(*A11); }
		JagVector2D& getA13() { return *A13; }
		Vector&      getA13( int i ) { return (*A13)(i); }
		Matrix&      getA21() { return *A21; }
		MatFact&     getA22() { return *A22; }
		JagVector2D& getA23() { return *A23; }
		Vector&      getA23( int i ) { return (*A23)(i); }
		JagVector2D& getA31() { return *A31; }
		Vector&      getA31( int i ) { return (*A31)(i); }
		JagVector2D& getA32() { return *A32; }
		Vector&      getA32( int i ) { return (*A32)(i); }
		Matrix&      getA33() { return *A33; }
		double&      getA33( int i, int j ) { return (*A33)(i,j); }
		
		inline void GMBE( Vector& X1, Vector& X2, double& X3, const Vector& F1, const Vector& F2, const double& F3 );
		inline void GMBEW( int bord, Vector& X1, Vector& X2, Vector& X3, const Vector& F1, const Vector& F2, const Vector& F3 );

		void AX( HyperVector& X, HyperVector& F, bool tan=true );
		
		void Solve( HyperVector& X, HyperVector& F );
		
		void Solve( HyperVector& X, HyperVector& F, int bord );
		
		void Check( HyperVector& X, HyperVector& F );
		
		void SolveDIRECT( HyperVector& X, HyperVector& F );
		
		void Solve( Vector& x, Vector& f ) { A11->Solve( x, f ); }
		
		void Solve( Vector& x, double& z, const Vector& f, const double& h );
		
		void Solve( int bord, Vector& x, Vector& z, const Vector& f, const Vector& h ); // BEMW
		
		void Solve( Vector& x, Vector& y, double& z, const Vector& f, const Vector& g, const double& h ); // GMBE
		
		void Solve( int bord, Vector& x, Vector& y, Vector& z, const Vector& f, const Vector& g, const Vector& h ); // GMBEW
		
	protected:
	
		// matrix blocks
		
		SpFact*			A11;
		JagVector2D*	A13;
		Matrix*			A21;
		MatFact*			A22;
		JagVector2D*	A23;
		JagVector2D*	A31;
		JagVector2D*	A32;
		Matrix*			A33;
		
	private:
		// temporary storage for factorization	
		
		// BEMW 1.
		JagVector2D*	VV1;
		JagVector2D*	VV1Star;
		JagVector2D*	FF1;
		Vector*			GG1;
		JagVector2D*	XX1;
		Vector*			YY1;
		Vector*			delta1;
		Vector*			delta1Star;
		
		// BEMW 2.
		JagVector2D*	VV2;
		JagVector2D*	VV2Star;
		JagVector2D*	FF2;
		Vector*			GG2;
		JagVector2D*	XX2;
		Vector*			YY2;
		Vector*			delta2;
		Vector*			delta2Star;

	
};
#endif
