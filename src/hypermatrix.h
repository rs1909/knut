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
	
		inline HyperVector( int i, int j, int k ) : V1(i), V2(j), V3(k) { }
		inline ~HyperVector() { }
		
		inline Vector& getV1() { return V1; }
		inline Vector& getV2() { return V2; }
		inline Vector& getV3() { return V3; }
		
	private:
	
		Vector V1;
		Vector V2;
		Vector V3;
};

#define FACT SpFact
class HyperMatrix
{
	public:
	
		// constructing Sparse-Dense type hypermatrix
		HyperMatrix( int i, int j, int k, int nz );
		~HyperMatrix();
		
		inline FACT&        getA11() { return *A11; }
		inline JagVector2D& getA13() { return *A13; }
		inline Vector&      getA13( int i ) { return (*A13)(i); }
		inline Matrix&      getA21() { return *A21; }
		inline MatFact&     getA22() { return *A22; }
		inline JagVector2D& getA23() { return *A23; }
		inline Vector&      getA23( int i ) { return (*A23)(i); }
		inline JagVector2D& getA31() { return *A31; }
		inline Vector&      getA31( int i ) { return (*A31)(i); }
		inline JagVector2D& getA32() { return *A32; }
		inline Vector&      getA32( int i ) { return (*A32)(i); }
		inline Matrix&      getA33() { return *A33; }
		inline double&      getA33( int i, int j ) { return (*A33)(i,j); }
		
		template< class T > inline void __BEM( T& _A, Vector& _b, Vector& _bStar, double& _d, Vector& x, double& y, const Vector& f, const double& g );
		template< class T > inline void __BEMWS
			(
				T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
				JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
				int j, const JagVector2D& V, const JagVector2D& VStar, const Vector& delta, const Vector& deltaStar, bool trans
			);
		template< class T > inline void __BEMWF
			(
				int bord,
				T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix& D,
				JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector& G,
				JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector& deltaStar
			);
		template< class FACT > inline void __BEMW
			(
				int bord,
				FACT&        A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
				JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
				JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector&       deltaStar,
				Vector&      x,  Vector&      y,      const Vector& f,     const Vector& g
			);
		template< class T > inline void GMBE( Vector& X1, Vector& X2, double& X3, const Vector& F1, const Vector& F2, const double& F3 );
		template< class T > inline void GMBEW( int bord, Vector& X1, Vector& X2, Vector& X3, const Vector& F1, const Vector& F2, const Vector& F3 );

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
		
		FACT*				A11;
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

		// BEM temporaries
		Vector         bem_v;
		Vector         bem_vStar;
		Vector         bem_w;
		Vector         bem_f1;
	
		// GMBE temporaries
		Vector         gmbe_xi;
		Vector         gmbe_cc;
		Vector         cm_fr;
	
		// GMBEW temporaries
		JagVector2D    gmbew_xi;
		JagVector2D    gmbew_cc;
		Matrix         gmbew_dd;
		Vector         gmbew_gg;
		Vector         gmbew_gr;
		Vector         gmbew_yy;
		Matrix         gmbew_m_beta;
	
};

#undef FACT

#endif
