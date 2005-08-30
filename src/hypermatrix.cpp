// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "matrix.h"
#include "spmatrix.h"
#include "hypermatrix.h"

#include "plot.h"
#include <cmath>

#define FACT_A11 SpFact
// constructing Sparse-Dense type hypermatrix
HyperMatrix :: HyperMatrix( int i, int j, int k, int nz ) :
	bem_v( i ),
	bem_vStar( i ),
	bem_w( i ),
	bem_f1( i ),
	gmbe_xi( j ),
	gmbe_cc( i ),
	cm_fr( j ),
	gmbew_xi( k, j ),
	gmbew_cc( k, i ),
	gmbew_dd( k, k ),
	gmbew_gg( k ),
	gmbew_gr( k ),
	gmbew_yy( k ),
	gmbew_m_beta( k, k )
{
	// constructing the data members
	if( i != 0 ) A11 = new FACT_A11( 'R', i, nz ); else A11 = 0;
	if( j != 0 ) A22 = new MatFact( j, j ); else A22 = 0;
	if( k != 0 ) A33 = new Matrix( k, k ); else A33 = 0;
	
	if( ( i != 0 )&&( j != 0 ) ) A21 = new Matrix( i, j ); else A21 = 0;
	if( ( i != 0 )&&( k != 0 ) )
	{
		A13 = new JagVector2D( k, i );
		A31 = new JagVector2D( k, i );
		
		VV1 = new JagVector2D( k+1 );
		VV1Star = new JagVector2D( k+1 );
		FF1 = new JagVector2D( k+1 );
		GG1 = new Vector( k+1 );
		XX1 = new JagVector2D( k+1 );
		YY1 = new Vector( k+1 );
		delta1 = new Vector( k+1 );
		delta1Star = new Vector( k+1 );
		
		for( int r=0; r<k+1; r++ )
		{
			(*VV1)(r).Init( i + r );
			(*VV1Star)(r).Init( i + r );
			(*XX1)(r).Init( i + r );
			(*FF1)(r).Init( i + r );
		}
	}else
	{
		A31 = 0;
		A13 = 0;
		
		VV1 = 0;
		VV1Star = 0;
		FF1 = 0;
		GG1 = 0;
		XX1 = 0;
		YY1 = 0;
		delta1 = 0;
		delta1Star = 0;
	}
	if( ( j != 0 )&&( k != 0 ) )
	{
		A23 = new JagVector2D( k, j );
		A32 = new JagVector2D( k, j );
		
		VV2 = new JagVector2D( k+1 );
		VV2Star = new JagVector2D( k+1 );
		FF2 = new JagVector2D( k+1 );
		GG2 = new Vector( k+1 );
		XX2 = new JagVector2D( k+1 );
		YY2 = new Vector( k+1 );
		delta2 = new Vector( k+1 );
		delta2Star = new Vector( k+1 );
					
		for( int r=0; r<k+1; r++ )
		{
			(*VV2)(r).Init( j + r );
			(*VV2Star)(r).Init( j + r );
			(*XX2)(r).Init( j + r );
			(*FF2)(r).Init( j + r );
		}
	}else
	{
		A32 = 0;
		A23 = 0;
		
		VV2 = 0;
		VV2Star = 0;
		FF2 = 0;
		GG2 = 0;
		XX2 = 0;
		YY2 = 0;
		delta2 = 0;
		delta2Star = 0;
	}
}

HyperMatrix :: ~HyperMatrix()
{
	// deleting the data members
	delete A11;
	delete A13;
	delete A21;
	delete A22;
	delete A23;
	delete A31;
	delete A32;
	delete A33;
	
	delete VV1;
	delete VV1Star;
	delete FF1;
	delete GG1;
	delete XX1;
	delete YY1;
	delete delta1;
	delete delta1Star;
	
	delete VV2;
	delete VV2Star;
	delete FF2;
	delete GG2;
	delete XX2;
	delete YY2;
	delete delta2;
	delete delta2Star;
}

// 	BEM temporaries
// 	Vector v( _b.Size() ); == A.Size()                     bem_v
// 	Vector vStar( _b.Size() );                             bem_vStar
// 	Vector w( _b.Size() );                                 bem_w
// 	Vector f1( _b.Size() );                                bem_f1;

// 	GMBE temporaries 
// 	Vector xi( _A32.Size() ); == dim2                      gmbe_xi;
// 	Vector cc( _A31 ); == dim1                             gmbe_cc
// 	Vector fr( F2 ); == dim2                               cm_fr

// GMBEW temporaries
// 	JagVector2D xi( bord, _A32(0).Size() );  -> different  gmbew_xi
// 	JagVector2D cc( _A31 );                  -> different  gmbew_cc
// 	Matrix dd( bord, bord );                               gmbew_dd
// 	Vector gg( bord );                                     gmbew_gg
// 	Vector fr( F2.Size() ); -> same as GMBE                cm_fr
// 	Vector gr( bord );                                     gmbew_gr
// 	Vector yy( bord );                                     gmbew_yy
// 	Matrix m_beta( bord, bord );                           gmbew_m_beta

// BEM
template< class FACT >
inline void HyperMatrix::__BEM( FACT& _A, Vector& _b, Vector& _bStar, double& _d, Vector& x, double& y, const Vector& f, const double& g )
{

	double delta, deltaStar;
	double g1;
	double y1, y2;
	
	// Doolittle
	_A.Solve( bem_vStar, _bStar, true );                              // Step 1.
	deltaStar = _d - bem_vStar*_b;                                    // Step 2. Scalar + ddot
	// Crout
	_A.Solve( bem_v, _b );                                            // Step 3.
	delta = _d - _bStar*bem_v;                                        // Step 4. Scalar + ddot
	// approx Y
	y1 = (g - bem_vStar*f)/deltaStar;                                 // Step 5. Scalar + ddot
	// residuals
	bem_f1 = f;
	bem_f1 -= y1 * _b;                                                // Step 6. daxpy
// 	for( int i=0; i<_b.Size(); i++ ) bem_f1(i) = f(i) - _b(i)*y1;
	g1 = g - _d*y1;                                                   // Step 7. Scalar
	
	// residual corrections
	_A.Solve( bem_w, bem_f1 );                                        // Step 8.
	y2 = (g1 - _bStar*bem_w)/delta;                                   // Step 9. Scalar + ddot
	x = bem_w;
	bem_w -= y2 * bem_v;                                              // Step 10. daxpy
// 	for( int i=0; i<_b.Size(); i++ ) x(i) = bem_w(i) - bem_v(i)*y2;
	y = y1 + y2;                                                      // Step 11. Scalar

}

template< class FACT >
inline void HyperMatrix :: __BEMWS
	(
		FACT&        A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
		JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
		int j, const JagVector2D& V, const JagVector2D& VStar, const Vector& delta, const Vector& deltaStar, bool trans
	)
{

	for( int k = j; k > 0; k-- )
	{
		const int len   = A.Col() + k-1;
		const int aalen = A.Col();
		const int ddlen = k-1;
		const int idx   = k-1;
		
		// Y_k = ( -v*_k 1 )^T . f_k
		Y(k) = F(k)( len );
		for( int r=0; r < len; r++ ) Y(k) -= VStar(k)( r )*F(k)( r ); // ddot with range
		Y(k) /= deltaStar(k);
		
		// residuals
		for( int r=0; r < aalen; r++ ) F(k-1)( r )         = F(k)( r )         - B(idx)( r )*Y(k); // ???
		for( int r=0; r < ddlen; r++ ) F(k-1)( aalen + r ) = F(k)( aalen + r ) - D( r, idx )*Y(k);
		G(k) = F(k)( len ) - D( idx, idx )*Y(k);
	}
	
	A.Solve( X(0), F(0), trans );
	double ypp;
	
	for( int k = 1; k <= j; k++ )
	{
		const int len   = A.Col() + k-1;
		const int aalen = A.Col();
		const int ddlen = k-1;
		const int idx   = k-1;
		
		ypp = G(k);
		for( int r=0; r < aalen; r++ ) ypp -= BStar(idx)( r )*X(k-1)( r );      // ddot with range
		for( int r=0; r < ddlen; r++ ) ypp -= D( idx, r )*X(k-1)( aalen + r );  // ddot with range this won't be vectorized.
		ypp /= delta(k);
		
		for( int r=0; r < len; r++ ) X(k)( r ) = X(k-1)( r ) - V(k)( r )*ypp;   // daxpy with range
		X(k)( len ) = Y(k) + ypp;
	}
}

template< class FACT >
inline void HyperMatrix :: __BEMWF
	(
		int bord,
		FACT&        A,  JagVector2D& B,      JagVector2D&  BStar, Matrix& D,
		JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector& G,
		JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector& deltaStar
	)

{
	for( int k = 1; k < bord+1; k++ )
	{
		const int len   = A.Col() + k-1;
		const int aalen = A.Col();
		const int ddlen = k-1;
		const int idx   = k-1;
		
		// with the ordinary solver
		for( int r=0; r < aalen; r++ ) F(k-1)( r )         = B(idx)( r );  // dcopy with range
		for( int r=0; r < ddlen; r++ ) F(k-1)( aalen + r ) = D( r, idx );  // dcopy with range
		
		__BEMWS<FACT>( A, B, BStar, D, X, Y, F, G, k-1, V, VStar, delta, deltaStar, false );
		
		// VV
		for( int r=0; r < len; r++ ) V(k)( r ) = X(k-1)( r ); // dcopy with range
		// delta
		delta(k) = D(idx,idx);
		for( int r=0; r < aalen; r++ ) delta(k) -= BStar(idx)( r )*V(k)( r );  // ddot with range
		for( int r=0; r < ddlen; r++ ) delta(k) -= D( idx, r )*V(k)( aalen + r );
		
		// with the adjoint solver
		for( int r=0; r < aalen; r++ ) F(k-1)( r )         = BStar(idx)( r );  // dcopy with range
		for( int r=0; r < ddlen; r++ ) F(k-1)( aalen + r ) = D( idx, r );
		
		__BEMWS<FACT>( A, B, BStar, D, X, Y, F, G, k-1, VStar, V, deltaStar, delta, true );
		
		// VV
		for( int r=0; r < len; r++ ) VStar(k)( r ) = X(k-1)( r );            // dcopy with range
		// delta
		deltaStar(k) = D(idx,idx);
		for( int r=0; r < aalen; r++ ) deltaStar(k) -= B(idx)( r )*VStar(k)( r );         // ddot with range
		for( int r=0; r < ddlen; r++ ) deltaStar(k) -= D( r, idx )*VStar(k)( aalen + r );
	}
}

// BEMW
template< class FACT >
void inline HyperMatrix :: __BEMW
	(
		int bord,
		FACT&        A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
		JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
		JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector&       deltaStar,
		Vector&      x,  Vector&      y,      const Vector& f,     const Vector& g
	)
{
	__BEMWF<FACT>( bord, A, B, BStar, D, X, Y, F, G, V, VStar, delta, deltaStar );
	
	for( int i = 0; i < A.Row(); i++ ) F(bord)( i )           = f( i ); // dcopy with range : F(bord)[rng(0,A.Row)] = f;
	for( int i = 0; i < bord; i++ )    F(bord)( A.Row() + i ) = g( i ); // dcopy with range : F(bord)[rng(A.Row,A.Row+bord)] = g[rng(0,bord)];
	
	__BEMWS<FACT>( A, B, BStar, D, X, Y, F, G, bord, V, VStar, delta, deltaStar, false );
	
	for( int i = 0; i < A.Row(); i++ ) x(i) = X(bord)( i );             // dcopy with range : x = X(bord)[rng(0,A.Row)];
	for( int i = 0; i < bord; i++ )    y(i) = X(bord)( A.Row() + i );   // dcopy with range : y = X(bord)[rng(A.Row,A.Row+bord)];
}


template< class FACT >
inline void HyperMatrix :: GMBE( Vector& X1, Vector& X2, double& X3, const Vector& F1, const Vector& F2, const double& F3 )
{
	FACT&       _A11 = *A11;
	Vector&     _A13 = (*A13)(0);
	Matrix&     _A21 = *A21;
	MatFact&    _A22 = *A22;
	Vector&     _A23 = (*A23)(0);
	Vector&     _A31 = (*A31)(0);
	Vector&     _A32 = (*A32)(0);
	double&     _A33 = (*A33)(0);
	
	double dd;
	double gg;
	
	_A22.Solve( gmbe_xi, _A32, true );
	// gmbe_cc = A31 - gmbe_xi*A21
	gmbe_cc = _A31;
	gmbe_cc -= _A21*gmbe_xi; // 1. mpm gmbe_cc writeable, so initialize with _A31 !!!
	// 	_A21.AXpY( gmbe_cc, gmbe_xi, _A31, -1.0, 1.0, false ); 
	dd = _A33 - gmbe_xi*_A23;
	gg = F3 - gmbe_xi*F2;
	
	__BEM<FACT>( _A11, _A13, gmbe_cc, dd, X1, X3, F1, gg );
	
	double gr;
	double yy;
	double beta = 1.0;
	double m_beta = -beta;
	
	gmbe_xi *= beta;
	
	// cm_fr = F2 - A21*X1 - A23*X3
	cm_fr = F2;
	cm_fr -= !_A21*X1;
	cm_fr -= X3*_A23;
// 	_A21.AXpY( cm_fr, X1, _A23, -1.0, -X3, true );
	gr = F3 - _A31*X1 - _A33*X3;
	
	__BEM<MatFact>( _A22, gmbe_xi, _A32, m_beta, X2, yy, cm_fr, gr );
}

template< class FACT >
inline void HyperMatrix :: GMBEW( int bord, Vector& X1, Vector& X2, Vector& X3, const Vector& F1, const Vector& F2, const Vector& F3 )
{
	FACT&        _A11 = *A11;
	JagVector2D& _A13 = *A13;
	Matrix&      _A21 = *A21;
	MatFact&     _A22 = *A22;
	JagVector2D& _A23 = *A23;
	JagVector2D& _A31 = *A31;
	JagVector2D& _A32 = *A32;
	Matrix&      _A33 = *A33;
 
	for( int i=0; i < bord; i++ )
	{
		_A22.Solve( gmbew_xi(i), _A32(i), true );
		// gmbew_cc = A31 - gmbew_xi*A21
		gmbew_cc(i) = _A31(i);
		gmbew_cc(i) -= _A21*gmbew_xi(i);
// 		_A21.AXpY( gmbew_cc(i), gmbew_xi(i), _A31(i), -1.0, 1.0, false ); 
		for( int j=0; j < bord; j++ )
		{
			gmbew_dd( i, j ) = _A33( i, j ) - gmbew_xi(i)*_A23(j);         // gmbew_dd = A33 - gmbew_xi*A23;
		}
		gmbew_gg(i) = F3( i ) - gmbew_xi(i)*F2;                           // gmbew_gg = F3 - gmbew_xi*F2;
	}
	
	__BEMW<FACT>( bord, _A11, _A13, gmbew_cc, gmbew_dd, *XX1, *YY1, *FF1, *GG1, *VV1, *VV1Star, *delta1, *delta1Star, X1, X3, F1, gmbew_gg );
	
	double beta = 1.0;
	
	for( int i=0; i < bord; i++ )
	{
		gmbew_xi(i) *= beta;
	}
	for( int i=0; i < bord; i++ ) gmbew_m_beta(i,i) = -beta;
	
	// cm_fr = F2 - A21*x1 - A23*X3
	cm_fr = F2;
	cm_fr = -1.0*!_A21*X1;
	// _A21.AX( cm_fr, X1, -1.0, true );
	for( int i=0; i < cm_fr.Size(); i++ ) // daxpy with range!
	{
		cm_fr(i) += F2(i);
		for( int j=0; j < bord; j++ ) cm_fr(i) -= _A23(j)(i)*X3(j); // cm_fr
	}
	
	// gmbew_gr = F3 - A31*X1 - A33*X3;
	for( int i=0; i < bord; i++ )
	{
		gmbew_gr(i) = F3(i) - _A31(i)*X1;
		for( int j=0; j < bord; j++ ) gmbew_gr(i) -= _A33(i,j)*X3(j);
	}

	__BEMW<MatFact>( bord, _A22, gmbew_xi, _A32, gmbew_m_beta, *XX2, *YY2, *FF2, *GG2, *VV2, *VV2Star, *delta2, *delta2Star, X2, gmbew_yy, cm_fr, gmbew_gr );

}

// Wrapper functions
void HyperMatrix::Solve( Vector& x, double& z, const Vector& f, const double& h ) // BEM
{
	__BEM<FACT_A11>( *A11, (*A13)(0), (*A31)(0), (*A33)(0,0), x, z, f, h );
}

void HyperMatrix::Solve( int bord, Vector& X1, Vector& X3, const Vector& F1, const Vector& F3 ) // BEMW
{
	__BEMW<FACT_A11>( bord, *A11, *A13, *A31, *A33, *XX1, *YY1, *FF1, *GG1, *VV1, *VV1Star, *delta1, *delta1Star, X1, X3, F1, F3 );
}

void HyperMatrix::Solve( Vector& x, Vector& y, double& z, const Vector& f, const Vector& g, const double& h ) // GMBE
{
	GMBE<FACT_A11>( x, y, z, f, g, h );
}

void HyperMatrix::Solve( int bord, Vector& X1, Vector& X2, Vector& X3, const Vector& F1, const Vector& F2, const Vector& F3 ) // GMBEW
{
	GMBEW<FACT_A11>( bord, X1, X2, X3, F1, F2, F3 );
}

void HyperMatrix::Check( HyperVector& X, HyperVector& F )
{
	// this may be buggy with some compilers
	
	FACT_A11&    _A11 = *A11;
	JagVector2D& _A13 = *A13;
	Matrix&      _A21 = *A21;
	MatFact&     _A22 = *A22;
	JagVector2D& _A23 = *A23;
	JagVector2D& _A31 = *A31;
	JagVector2D& _A32 = *A32;
	Matrix&      _A33 = *A33;
	Vector&      X1   = X.getV1();
	Vector&      X2   = X.getV2();
	Vector&      X3   = X.getV3();
	Vector&      F1   = F.getV1();
	Vector&      F2   = F.getV2();
	Vector&      F3   = F.getV3();
	
	//multiply back...
	Vector R1( X1.Size() );
	Vector R2( X2.Size() );
	Vector R2_t( X2.Size() );
	Vector R3( X3.Size() );
	Vector R3_t( X3.Size() );
	
	if( A11 )
	{
		R1 = _A11 * X1;
// 		_A11.AX( R1, X1 );
		if( A33 )
		{
			for( int i=0; i<X1.Size(); i++ )
			{
				for( int j=0; j<X3.Size(); j++ )
				{
					R1(i) += _A13(j)(i) * X3(j);
				}
				R1(i) -= F1(i);
			}
		}
		std::cout<<"F1: "<<sqrt(F1*F1)<<"R1: "<<sqrt(R1*R1)<<"\n";
	}
	
	if( A22 )
	{
		if( A11 ) R2 = !_A21 * X1; // _A21.AX( R2, X1, 1.0, true );
		R2_t = _A22 * X2; // _A22.AX( R2_t, X2 );
		R2 += R2_t;
		if( A33 )
		{
			for( int i=0; i<X2.Size(); i++ )
			{
				for( int j=0; j<X3.Size(); j++ )
				{
					R2(i) += _A23(j)(i) * X3(j);
				}
				R2(i) -= F2(i);
			}
		}
		std::cout<<"F2: "<<sqrt(F2*F2)<<"R2: "<<sqrt(R2*R2)<<"\n";
	}
	
	for( int i=0; i<F3.Size(); i++)
	{
		R3(i) = - F3(i);
		if( A11 ) R3(i) += _A31(i)*X1;
		if( A22 ) R3(i) += _A32(i)*X2;
	}
	if( A33 )
	{
		R3_t = _A33 * X3; //_A33.AX( R3_t, X3, 1.0, false );
		R3 += R3_t;
		std::cout<<"F3: "<<sqrt(F3*F3)<<"R3: "<<sqrt(R3*R3)<<"\n";
	}
}


void HyperMatrix::Solve( HyperVector& X, HyperVector& F )
{
	if( A11 != 0 )
	{
		if( A22 != 0 )
		{
			if( A33 != 0 )
			{
				if( A33->Col() == 1 )
				{
					Solve( X.getV1(), X.getV2(), X.getV3()(0), F.getV1(), F.getV2(), F.getV3()(0) );
				}else
				{
					Solve( A33->Col(), X.getV1(), X.getV2(), X.getV3(), F.getV1(), F.getV2(), F.getV3() );
				}
			}
			else
			{
				// error
				std::cout<<"HyperMatrix::Solve Error\n";
				PDError(-1);
			}
		}
		else
		{
			if( A33 != 0 )
			{
				// MBEW v BEM
				if( A33->Col() == 1 )
				{
					Solve( X.getV1(), X.getV3()(0), F.getV1(), F.getV3()(0) );
				}
				else{
					Solve( A33->Col(), X.getV1(), X.getV3(), F.getV1(), F.getV3() );
				}
			}
			else
			{
				// simply A11
				A11->Solve( X.getV1(), F.getV1() );
			}
		}
	}
	else
	{
		std::cout<<"HyperMatrix::Solve Error\n";
		PDError(-1);
	}
	
//	Check( X, F );
}

void HyperMatrix::Solve( HyperVector& X, HyperVector& F, int bord )
{
	if( A11 != 0 )
	{
		if( A22 != 0 )
		{
			if( A33 != 0 )
			{
				if( bord == 1 )
				{
					Solve( X.getV1(), X.getV2(), X.getV3()(0), F.getV1(), F.getV2(), F.getV3()(0) );
				}else
				{
					Solve( bord, X.getV1(), X.getV2(), X.getV3(), F.getV1(), F.getV2(), F.getV3() );
				}
			}
			else
			{
				// error
				std::cout<<"HyperMatrix::Solve Error\n";
				PDError(-1);
			}
		}
		else
		{
			if( A33 != 0 )
			{
				// MBEW v BEM
				if( bord == 1 )
				{
					Solve( X.getV1(), X.getV3()(0), F.getV1(), F.getV3()(0) );
				}
				else{
					Solve( bord, X.getV1(), X.getV3(), F.getV1(), F.getV3() );
				}
			}
			else
			{
				// simply A11
				A11->Solve( X.getV1(), F.getV1() );
			}
		}
	}
	else
	{
		std::cout<<"HyperMatrix::Solve Error\n";
		PDError(-1);
	}
	
//	Check( X, F );
}

void HyperMatrix::SolveDIRECT( HyperVector& X, HyperVector& F )
{	
	int dim1, dim2, dim3;
	if( A11 ) dim1 = A11->Col(); else dim1 = 0;
	if( A22 ) dim2 = A22->Col(); else dim2 = 0;
	if( A33 ) dim3 = A33->Col(); else dim3 = 0;
	
	SpFact* A11S = A11;
	
	SpFact  AA( 'R', dim1+dim2+dim3, A11S->GetNZ() + (dim1+dim2+dim3)*(dim2+dim3) + dim1*dim3 + 10 );
	Vector  XX( dim1+dim2+dim3 ), FF( dim1+dim2+dim3 );
	
	for( int i=0; i<dim1; i++ )
	{
		AA.NewL( A11S->GetL(i) + dim3 );
		for( int j = 0; j < A11S->GetL(i); j++ )
		{
			AA.WrLi( i, j ) = A11S->WrLi( i, j );
			AA.WrLx( i, j ) = A11S->WrLx( i, j );
		}
		for( int j = 0; j < dim3; j++ )
		{
			AA.WrLi( i, A11S->GetL(i) + j ) = dim1+dim2+j;
			AA.WrLx( i, A11S->GetL(i) + j ) = getA13(j)(i);
		}
		FF( i ) = F.getV1()(i);
	}
	for( int i=0; i<dim2; i++ )
	{
		AA.NewL( dim1 + dim2 + dim3 );
		for( int j=0;j<dim1;j++ )
		{
			AA.WrLi( dim1 + i, j ) = j;
			AA.WrLx( dim1 + i, j ) = getA21()(j,i);	
		}
		for( int j=0;j<dim2;j++ )
		{
			AA.WrLi( dim1 + i, dim1 + j ) = dim1 + j;
			AA.WrLx( dim1 + i, dim1 + j ) = getA22()(i,j);
		}
		for( int j=0;j<dim3;j++ )
		{
			AA.WrLi( dim1 + i, dim1 + dim2 + j ) = dim1 + dim2 + j;
			AA.WrLx( dim1 + i, dim1 + dim2 + j ) = getA23(j)(i);
		}
		FF( dim1 + i ) = F.getV2()(i);
	}
	for( int i=0; i<dim3; i++ )
	{
		AA.NewL( dim1 + dim2 + dim3 );
		for( int j=0;j<dim1;j++ )
		{
			AA.WrLi( dim1 + dim2 + i, j ) = j;
			AA.WrLx( dim1 + dim2 + i, j ) = getA31(i)(j);
		}
		for( int j=0;j<dim2;j++ )
		{
			AA.WrLi( dim1 + dim2 + i, dim1 + j ) = dim1 + j;
			AA.WrLx( dim1 + dim2 + i, dim1 + j ) = getA32(i)(j);
		}
		for( int j=0;j<dim3;j++ )
		{
			AA.WrLi( dim1 + dim2 + i, dim1 + dim2 + j ) = dim1 + dim2 + j;
			AA.WrLx( dim1 + dim2 + i, dim1 + dim2 + j ) = getA33(i,j);
		}
		FF( dim1 + dim2 + i ) = F.getV3()(i);
	}
	
	AA.Solve( XX, FF );
	
	for( int i=0; i<dim1; i++ )  X.getV1()(i) = XX( i );
	for( int i=0; i<dim2; i++ )  X.getV2()(i) = XX( dim1 + i );
	for( int i=0; i<dim3; i++ )  X.getV3()(i) = XX( dim1 + dim2 + i );
	
//	Check( X, F );
	
}

// currently tan is not implemented: it is assumed to be false
// tan == false -> the last line not multiplied
// tan == true  -> normal multiplication
void HyperMatrix::AX( HyperVector& out, HyperVector& in, bool /*tan*/ )
{
	FACT_A11&    _A11 = *A11;
	JagVector2D& _A13 = *A13;
	Matrix&      _A21 = *A21;
	MatFact&     _A22 = *A22;
	JagVector2D& _A23 = *A23;
	JagVector2D& _A31 = *A31;
	JagVector2D& _A32 = *A32;
	Matrix&      _A33 = *A33;
	Vector&      in1   = in.getV1();
	Vector&      in2   = in.getV2();
	Vector&      in3   = in.getV3();
	Vector&      out1  = out.getV1();
	Vector&      out2  = out.getV2();
	Vector&      out3  = out.getV3();
	
	//multiply back...
	Vector R2_t( in2.Size() );
	Vector R3_t( in3.Size() );
	
	// first line: no A12!
	if( A11 )
	{
		out1 = _A11 * in1; //_A11.AX( out1, in1 );
		if( A33 )
		{
			for( int i=0; i<in1.Size(); i++ )
			{
				for( int j=0; j<in3.Size(); j++ )
				{
					out1(i) += _A13(j)(i) * in3(j);
				}
			}
		}
	}
	// second line
	if( A22 )
	{
		if( A11 ) out2 = !_A21 * in1; // _A21.AX( out2, in1, 1.0, true ); 
		else out2.Clear();
		R2_t = _A22 * in2; // _A22.AX( R2_t, in2 );
		out2 += R2_t;
		if( A33 )
		{
			for( int i=0; i<in2.Size(); i++ )
			{
				for( int j=0; j<in3.Size(); j++ )
				{
					out2(i) += _A23(j)(i) * in3(j);
				}
			}
		}
	}
	// third line
	for( int i=0; i<in3.Size(); i++)
	{
		out3(i) = 0.0;
		if( A11 ) out3(i) += _A31(i)*in1;
		if( A22 ) out3(i) += _A32(i)*in2;
	}
	if( A33 )
	{
		R3_t = _A33 * in3; // _A33.AX( R3_t, in3, 1.0, false );
		out3 += R3_t;
	}
}
