// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#ifndef MATRIX_H
#define MATRIX_H

// from f2c.h
typedef long int integer;
typedef long int logical;
typedef long int ftnlen;
typedef double doublereal;

#define TRUE_ (1)
#define FALSE_ (0)
// end - from f2c.h

#include <iostream>

#ifndef PDDESYS_H
#include "plot.h"
#include "error.h"
extern "C" {
#include "cspblas.h"
}
#endif

using namespace std;


template<class T>
class Array1D{

protected:

  int n;
  T* v;

public:

	inline Array1D() { n = 0; v = 0; }

	inline Array1D( int i )
	{
		n = i;
		v = new T[i];
		Clear();
	}
	
	// specially for JagMatrix2D i.e. Array1D< Vector >, which is indexed as (i)(j)
	inline Array1D( int i, int j )
	{
		n = i;
		v = new T[i];
		for( int r=0; r<i; r++ ) v[r].Init(j);
	}
	
	// specially for JagMatrix3D i.e. Array1D< Matrix >, which is indexed as (i)(j,k)
	inline Array1D( int i, int j, int k )
	{
		n = i;
		v = new T[i+1];
		for( int r=0; r<i; r++ ) v[r].Init(j,k);
	}

	inline Array1D( const Array1D<T> &V )
	{
		n = V.n;
		v = new T[V.n];
		for( int i=0; i<V.n; i++ ) v[i]=V.v[i];
	}
	
	inline virtual ~Array1D() { delete[] v; }
	
	inline void Init( int i )
	{
		delete v;
		n = i;
		v = new T[i];
		Clear();
	}

	inline void Clear( ) { for(int j=0; j<n; j++) v[j]=T(); } // it was T(0), but for built in types it must be the same
	
	inline int Size() const { return n; }
	
	inline Array1D<T>& operator=( const Array1D<T>& V )
	{
	 #ifdef DEBUG
		if( n != V.n )
		{
			std::cout<<"Array1D::operator=(): incompatible sizes\n";
			PDError(-1);
		}
	 #endif
		for( int i=0; i<V.n; i++ ) v[i]=V.v[i];
		return *this;
	}

	inline T& operator()( int i )
	{
	 #ifdef DEBUG
		if( i>=n ){ std::cout<<"Array1D::bound11& "<<n<<", "<<i<<"\n"; PDError(12); }
	 #endif
		return v[i];
	}

	inline T operator()( int i ) const
	{
	 #ifdef DEBUG
		if( i>=n ){ std::cout<<"vec_bound11_\n"; PDError(12); }
	 #endif
		return v[i];
	}
};


template<class T>
class Array2D{

 private:
	T* v;
 protected:

	T* m;
	int r, c;

 public:
	inline Array2D( ) { r = 0; c = 0; v = 0; m = 0; }

	inline Array2D( int _r, int _c ) : r(_r), c(_c)
	{
		v = new T[r*c+2];
		if( (unsigned int)v % (2*sizeof(T)) == 0 ) m = v; else m = v + 1;
		Clear();
	}

	inline Array2D( const Array2D<T>& M ) : r(M.r), c(M.c)
	{
		v = new T[r*c+2];
		if( (unsigned int)v % (2*sizeof(T)) == 0 ) m = v; else m = v + 1;
		for( int i = 0; i < r*c; i++ ) m[i] = M.m[i];
	}

	inline virtual ~Array2D(){ delete []v; }

	inline void Init( int _r, int _c )
	{
		delete v;
		r = _r;
		c = _c;
		v = new T[r*c+2];
		if( (unsigned int)v % (2*sizeof(T)) == 0 ) m = v; else m = v + 1;
		Clear();
	}

	inline void Clear() { for( int i = 0; i < r*c; i++ ) m[i] = T(0); }

	inline Array2D<T>& operator= ( const Array2D<T>& M )
	{
		if( (M.r == r)&&(M.c == c) )
		{
			for( int i = 0; i < r*c; i++ ) m[i] = M.m[i];
		} else
		{
			cout<<"Array2D<T>::operator= : incompatible sizes\n";
		}
		return *this;
	} 

	inline T& operator()( const int i, const int j )
	{
	 #ifdef DEBUG
		if((i>=r)||(j>=c)){ cout<<"bound& "<<r<<", "<<c<<",-"<<i<<", "<<j<<"\n"; PDError(1); }
		if((i<0)||(j<0)){ cout<<"lbound& "<<i<<", "<<j<<"\n"; PDError(1); }
	 #endif
		return m[i + r*j];
	}

	inline T operator()( const int i, const int j ) const
	{
	 #ifdef DEBUG
		if((i>=r)||(j>=c)){ cout<<"bound_\n"; PDError(1); }
		if((i<0)||(j<0)){ cout<<"lbound_\n"; PDError(1); }
	 #endif
		return m[i + r*j];
	}

	inline T& operator()( int i )
	{
	 #ifdef DEBUG
		if( i>=r*c ){ cout<<"bound11&\n"; PDError(1); }
		if((r!=1)&&(c!=1)){ cout<<"H&\n"; PDError(1); }
	 #endif
		return m[i];
	}

	inline T operator()( int i ) const
	{
	 #ifdef DEBUG
		if( i>=r*c ){ cout<<"bound11_\n"; PDError(1); }
		if((r!=1)&&(c!=1)){ cout<<"H&\n"; PDError(1); }
	 #endif
		return m[i];
	}

};

template<class T>
class Array3D{

 protected:

	T* m;
	int d1,d2,d3;

 public:
	inline Array3D( ) { d1 = d2 = d3 = 0; m = 0; }

	inline Array3D( int _d1, int _d2, int _d3 ) : d1(_d1), d2(_d2), d3(_d3)
	{
		m = new T[d1*d2*d3+1];
		Clear();
	}

	inline Array3D( const Array3D<T>& M ) : d1(M.d1), d2(M.d2), d3(M.d3)
	{
		m = new T[d1*d2*d3+1];
		for( int i = 0; i < d1*d2*d3; i++ ) m[i] = M.m[i];
	}

	inline virtual ~Array3D(){ delete []m; }

	inline void Init( int _d1, int _d2, int _d3 )
	{
		delete m;
		d1 = _d1;
		d2 = _d2;
		d3 = _d3;
		m = new T[d1*d2*d3+1];
		Clear();
	}
	
	inline void Clear() { for( int i = 0; i < d1*d2*d3; i++ ) m[i] = 0; }
	
	inline Array3D<T>& operator= ( const Array3D<T>& M )
	{
		if( (M.d1 == d1)&&(M.d2 == d2)&&(M.d3 == d3) )
		{
			for( int i = 0; i < d1*d2*d3; i++ ) m[i] = M.m[i];
		} else 
		{
			cout<<"Array3D<T>::operator= : incompatible sizes\n";
		}
		return *this;
	}
	
	inline T& operator()( const int i, const int j, const int k )
	{
	 #ifdef DEBUG
		if((i>=d1)||(j>=d2)||(k>=d3)) cout<<"bound&\n";
		if((i<0)||(j<0)||(k<0)) cout<<"lbound&\n";
	 #endif
		return m[i + d1*(j + d2*k)];
	}

	inline T operator()( const int i, const int j, const int k ) const
	{
	 #ifdef DEBUG
		if((i>=d1)||(j>=d2)||(k>=d3)) cout<<"bound&\n";
		if((i<0)||(j<0)||(k<0)) cout<<"lbound&\n";
	 #endif
		return m[i + d1*(j + d2*k)];
	}
	
};

#ifndef PDDESYS_H

class SpMatrix;
class Matrix;
class Vector;

template< class MT >           class __scal_vec_trans;
template< class MT, class VT > class __op_mul_vec;
template< class MT, class VT > class __op_mul_vec_plus_vec;

template< class MT >           class __scal_vec_trans_rng;
template< class MT, class VT > class __op_mul_vec_rng;
template< class MT, class VT > class __op_mul_vec_plus_vec_rng;

class rng
{
 public:
	inline rng( ) : i1(0), i2(0), j1(0), j2(0) { }
	inline rng( int i1_, int i2_ ) : i1(i1_), i2(i2_), j1(0), j2(0) { }
	inline rng( int i1_, int i2_, int j1_, int j2_ ) : i1(i1_), i2(i2_), j1(j1_), j2(j2_) { }
	
	const int i1, i2, j1, j2;
};

/// with range
template< class VT > class __scal_vec_trans_rng : public rng
{
 public:
	inline __scal_vec_trans_rng( const VT& v )                                        : rng(0, v.Row(), 0, v.Col()), vec( v ), alpha( 1.0 ), tr( NoTrans ) { }
	inline __scal_vec_trans_rng( const VT& v, double a )                              : rng(0, v.Row(), 0, v.Col()), vec( v ), alpha( a ),   tr( NoTrans ) { }
	inline __scal_vec_trans_rng( const VT& v, double a, enum cspblas_Trans t )        : rng(0, v.Row(), 0, v.Col()), vec( v ), alpha( a ),   tr( t ) { }
	inline __scal_vec_trans_rng( const VT& v, rng r )                                 : rng( r ), vec( v ), alpha( 1.0 ), tr( NoTrans ) { }
	inline __scal_vec_trans_rng( const VT& v, rng r, double a )                       : rng( r ), vec( v ), alpha( a ),   tr( NoTrans ) { }
	inline __scal_vec_trans_rng( const VT& v, rng r, double a, enum cspblas_Trans t ) : rng( r ), vec( v ), alpha( a ),   tr( t ) { }
	inline __scal_vec_trans_rng( __scal_vec_trans<VT> v );
	
	inline __scal_vec_trans_rng< VT >     operator!( ) { return __scal_vec_trans_rng<VT>( alpha, *this, Trans ); }
	inline __scal_vec_trans_rng< VT >     operator-( ) { return __scal_vec_trans_rng<VT>( -alpha, *this, tr ); }
	inline __op_mul_vec_rng< VT, Vector > operator*( __scal_vec_trans_rng<Vector> v );
	inline __op_mul_vec_rng< VT, Matrix > operator*( __scal_vec_trans_rng<Matrix> v );
	inline __op_mul_vec_rng< VT, Vector > operator*( const Vector& v );
	inline __op_mul_vec_rng< VT, Matrix > operator*( const Matrix& v );
	
	const VT& vec;
	const double alpha;
	const enum cspblas_Trans tr;
};

/// with range
template< class MT, class VT > class __op_mul_vec_rng
{
 public:
	inline __op_mul_vec_rng( __scal_vec_trans_rng<MT> m, __scal_vec_trans_rng<VT> v ) : op( m ), vecA( v ) { }
	
	inline __op_mul_vec_plus_vec_rng<MT,VT>   operator+( __scal_vec_trans_rng<VT> v );
	inline __op_mul_vec_plus_vec_rng<MT,VT>   operator+( const VT& v );
	inline __op_mul_vec_plus_vec_rng<MT,VT>   operator+( __scal_vec_trans<VT> v );
	inline __op_mul_vec_plus_vec_rng<MT,VT>   operator-( __scal_vec_trans_rng<VT> v );
	inline __op_mul_vec_plus_vec_rng<MT,VT>   operator-( const VT& v );
	inline __op_mul_vec_plus_vec_rng<MT,VT>   operator-( __scal_vec_trans<VT> v );
	
	const __scal_vec_trans_rng<MT> op;
	const __scal_vec_trans_rng<VT> vecA;
};

/// with range: this is the output
template< class MT, class VT > class __op_mul_vec_plus_vec_rng : public __op_mul_vec_rng<MT,VT>
{
 public:
	inline __op_mul_vec_plus_vec_rng( __op_mul_vec_rng<MT,VT> m, __scal_vec_trans_rng<VT> v ) : __op_mul_vec_rng<MT,VT>( m ), vecB( v ) { }
	
	const __scal_vec_trans_rng<VT> vecB;
};

/// WITHOUT RANGES
template< class VT > class __scal_vec_trans
{
 public:
	inline __scal_vec_trans( const VT& m ) : vec(m), alpha( 1.0 ), tr( NoTrans ) { }
	inline __scal_vec_trans( const VT& m, double a ) : vec( m ), alpha( a ), tr( NoTrans ) { }
	inline __scal_vec_trans( const VT& m, double a, enum cspblas_Trans t ) : vec( m ), alpha( a ), tr(t) { }
	
	inline __scal_vec_trans< VT >                        operator!( ) { return __scal_vec_trans<VT>( alpha, vec, Trans ); }
	inline __scal_vec_trans< VT >                        operator-( ) { return __scal_vec_trans<VT>( -alpha, vec, tr ); }
	inline __op_mul_vec< VT, Vector >                    operator*( const Vector& v );
	inline __op_mul_vec< VT, Matrix >                    operator*( const Matrix& v );
	inline __op_mul_vec_rng< VT, Vector >                operator*( __scal_vec_trans_rng<Vector> v );
	inline __op_mul_vec_rng< VT, Matrix >                operator*( __scal_vec_trans_rng<Matrix> v );
	
	const VT& vec;
	double alpha;
	const enum cspblas_Trans tr;
};

template< class MT, class VT > class __op_mul_vec
{
 public:
	inline __op_mul_vec( __scal_vec_trans<MT> m, __scal_vec_trans<VT> v ) : op( m ), vecA( v ) { }
	
	inline __op_mul_vec_plus_vec<MT,VT>      operator+( __scal_vec_trans<VT> v );
	inline __op_mul_vec_plus_vec<MT,VT>      operator-( __scal_vec_trans<VT> v );
	inline __op_mul_vec_plus_vec_rng<MT,VT>  operator+( __scal_vec_trans_rng<VT> v );
	inline __op_mul_vec_plus_vec_rng<MT,VT>  operator-( __scal_vec_trans_rng<VT> v );
	inline __op_mul_vec_plus_vec<MT,VT>      operator+( const VT& v );
	inline __op_mul_vec_plus_vec<MT,VT>      operator-( const VT& v );

	const __scal_vec_trans<MT> op;
	const __scal_vec_trans<VT> vecA;
};

template< class MT, class VT > class __op_mul_vec_plus_vec : public __op_mul_vec<MT,VT>
{
 public:
	inline __op_mul_vec_plus_vec( __op_mul_vec<MT,VT> m, __scal_vec_trans<VT> v ) : __op_mul_vec<MT,VT>( m ), vecB( v ) { }
	
	const __scal_vec_trans<VT> vecB;
};

#endif // PDDESYS_H

class Vector : public Array1D<double>
{

public:

	inline Vector() { }

	inline Vector( int i ) : Array1D<double>( i ) { }

	inline Vector( const Vector& V ) : Array1D<double>( V ) { }

	inline ~Vector() { }

	inline double *Pointer(){ return v; }
	
	inline void Print() const
	{
		for(int j=0; j<n; j++) std::cout<<v[j]<<'\t'; std::cout<<'\n';
	}

#ifndef PDDESYS_H
	
	inline Vector& operator=( const Vector& V )
	{
// 		std::cout<<"op= ";
	#ifdef DEBUG
		if( n != V.n ){ std::cout<<"Vector::operator=(): incompatible sizes\n"; }
	#endif // DEBUG
// 		for( int i=0; i < n; i++ ) v[i] = V.v[i];
		cblas_dnvcopy( n, v, V.v, 1.0 );
		return *this;
	}

	inline Vector& operator-=( const Vector& V )
	{
// 		std::cout<<"op-= ";
	#ifdef DEBUG
		if( n != V.n ){ std::cout<<"Vector::operator-=(): incompatible sizes\n"; }
	#endif // DEBUG
		for( int i=0; i < n; i++ ) v[i] -= V.v[i];
// 		cblas_dnaxpy( n, v, v, 1.0, V.v, -1.0 );
		return *this;
	}
	
	inline Vector& operator+=( const Vector& V )
	{
// 		std::cout<<"op+= ";
		if( n != V.n ){ std::cout<<"Vector::operator+=(): incompatible sizes\n"; }
		for( int i=0; i < n; i++ ) v[i] += V.v[i];
// 		cblas_dnaxpy( n, v, v, 1.0, V.v, 1.0 );
		return *this;
	}
	
	inline Vector& operator/=( double div )
	{
// 		std::cout<<"op/= ";
// 		for( int i=0; i < n; i++ ) v[i] /= div;
		cblas_dnvcopy( n, v, v, 1.0/div );
		return *this;
	}
	
	inline Vector& operator*=( double mul )
	{
// 		for( int i=0; i < n; i++ ) v[i] *= mul;
// 		std::cout<<"op*= ";
		cblas_dnvcopy( n, v, v, mul );
		return *this;
	}

	
	inline double operator*( const Vector& V ) const
	{
	#ifdef DEBUG
		if( n != V.n ){ std::cout<<"Vector::operator*(): incompatible sizes\n"; }
	#endif // DEBUG
		return cblas_dndot( n, this->v, V.v, 1.0 );
	}
	
	// Matrix operations
	inline Vector& operator=( const __scal_vec_trans<Vector> );
	inline Vector& operator=( const __op_mul_vec<Matrix,Vector> );
	inline Vector& operator=( const __op_mul_vec_plus_vec<Matrix,Vector> );
	inline Vector& operator=( const __op_mul_vec<SpMatrix,Vector> );
	inline Vector& operator=( const __op_mul_vec_plus_vec<SpMatrix,Vector> );
	// with ranges
	inline Vector& operator=( const __scal_vec_trans_rng<Vector> );
	inline Vector& operator=( const __op_mul_vec_rng<Matrix,Vector> );
	inline Vector& operator=( const __op_mul_vec_plus_vec_rng<Matrix,Vector> );
	inline Vector& operator=( const __op_mul_vec_rng<SpMatrix,Vector> );
	inline Vector& operator=( const __op_mul_vec_plus_vec_rng<SpMatrix,Vector> );

	inline __scal_vec_trans<Vector>     operator-( ) { return __scal_vec_trans<Vector>( *this, -1.0 ); }
	inline __scal_vec_trans_rng<Vector> operator[ ] ( rng r ) { return __scal_vec_trans_rng<Vector>( *this, r ); }
	
#endif // PDDESYS_H
	
	friend class SpMatrix;
	friend class SpFact;
	friend class Matrix;
	friend class MatFact;

};

class Matrix : public Array2D<double>
{
 public:

	inline Matrix() { }

	inline Matrix( int i, int j ) : Array2D<double>( i, j ) { }

	inline Matrix( const Matrix& M ) : Array2D<double>( M ) { }

	inline virtual ~Matrix() { }

	inline int Row() const { return r; }
	inline int Col() const { return c; }
	inline int Size() const { if((r==1)||(c==1)){return r*c;}else{cout<<"Hs\n"; return 0;} }

#ifndef PDDESYS_H
	inline void mmx( enum cspblas_Trans trans, double* out, const double* in, double alpha ) const
	{
		cblas_mmx( trans, r, c,  m, r, out, in, alpha );
	}
	inline void mmxpy( enum cspblas_Trans trans, double* out, const double* in, double alpha, const double* Y, double beta ) const
	{
		cblas_mmxpy( trans, r, c, m, r, out, in, alpha, Y, beta );
	}
	inline void mmxm( enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha, int nrhs ) const
	{
		cblas_mmxm( trans, r, c, m, r, out, ldout, in, ldin, alpha, nrhs );
	}
	inline void mmxmpym( enum cspblas_Trans trans, double* out, int ldout, 
	              const double* in, int ldin, double alpha, const double* Y, int ldY, double beta, int nrhs ) const
	{
		cblas_mmxmpym( trans, r, c, m, r, out, ldout, in, ldin, alpha, Y, ldY, beta, nrhs );
	}
	
	inline Matrix& operator=( const __scal_vec_trans<Matrix> );
	inline Matrix& operator=( const __op_mul_vec<Matrix,Matrix> );
	inline Matrix& operator=( const __op_mul_vec_plus_vec<Matrix,Matrix> );
	inline Matrix& operator=( const __op_mul_vec<SpMatrix,Matrix> );
	inline Matrix& operator=( const __op_mul_vec_plus_vec<SpMatrix,Matrix> );
		// with ranges
	inline Matrix& operator=( const __scal_vec_trans_rng<Matrix> );
	inline Matrix& operator=( const __op_mul_vec_rng<Matrix,Matrix> );
	inline Matrix& operator=( const __op_mul_vec_plus_vec_rng<Matrix,Matrix> );
	inline Matrix& operator=( const __op_mul_vec_rng<SpMatrix,Matrix> );
	inline Matrix& operator=( const __op_mul_vec_plus_vec_rng<SpMatrix,Matrix> );
	
	inline void AX( double* out, const double* in, double alpha, bool trans ) const;
	inline void AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const;
	
	void Eigval( Vector& re, Vector& im );
	void Eigval( Vector& re, Vector& im, Matrix& lev, Matrix& rev );
	void StrPlot( GnuPlot& pl );
	
	/* operators */
	inline Matrix& operator*=( double mul )
	{
		for( int i=0; i < r*c; i++ ) m[i] *= mul;
		return *this;
	}
	
	inline __op_mul_vec<Matrix,Vector>     operator*( Vector& v )
		{ return __op_mul_vec<Matrix,Vector>( __scal_vec_trans<Matrix>( *this ), __scal_vec_trans<Vector>( v ) ); }
	inline __op_mul_vec<Matrix,Matrix>     operator*( Matrix& v )
		{ return __op_mul_vec<Matrix,Matrix>( __scal_vec_trans<Matrix>( *this ), __scal_vec_trans<Matrix>( v ) ); }
	inline __op_mul_vec_rng<Matrix,Vector> operator*( __scal_vec_trans_rng<Vector> v )
		{ return __op_mul_vec_rng<Matrix,Vector>( __scal_vec_trans_rng<Matrix>( *this ), v ); }
	inline __op_mul_vec_rng<Matrix,Matrix> operator*( __scal_vec_trans_rng<Matrix> v )
		{ return __op_mul_vec_rng<Matrix,Matrix>( __scal_vec_trans_rng<Matrix>( *this ), v ); }
	inline __scal_vec_trans<Matrix>        operator!( )
		{ return __scal_vec_trans<Matrix>( *this, 1.0, Trans ); }
	inline __scal_vec_trans_rng<Matrix>    operator[ ] ( rng r )
		{ return __scal_vec_trans_rng<Matrix>( *this, r ); }
	
#endif // PDDESYS_H
	
	inline void Print()
	{
		double sum=0.0;
		for( int i = 0; i < r; i++ )
		{
			for( int j = 0; j < c; j++ )
			{
				std::cout<<(*this)(i,j)<<'\t';
				sum += (*this)(i,j);
			}
			std::cout<<'\n';
		}
		cout<<"SUM: "<<sum<<'\n';
	}
	
	friend class Vector;
	friend class MatFact;
	friend class SpMatrix;
	friend class SpFact;
};

#ifndef PDDESYS_H

class MatFact : public Matrix
{

		bool     fact; // == true if factorized
		// for factorizig
		double*  mf;
		integer* ipiv;
		double*  work;
		integer* iwork;
		integer  info;
		// for solving
		double rcond;
		double ferr;  // in case of Vector these are just double
		double berr;

	public:

		inline MatFact( int _r, int _c ) : Matrix( _r, _c )
		{
			fact = false;
			mf = new double[r*c+1];
			ipiv = new integer[this->r+1];
			work = new double[4*(this->r)];
			iwork = new integer[this->r];
		}
		
		inline MatFact( const Matrix& M ) : Matrix( M )
		{
			fact = false;
			mf = new double[r*c+1];
			ipiv = new integer[this->r+1];
			work = new double[4*(this->r)];
			iwork = new integer[this->r];
		}
		
		inline ~MatFact(){ delete[] iwork; delete[] work; delete[] ipiv; delete[] mf; }

		inline MatFact& operator=( MatFact& M )
		{
			Matrix::operator=( M );
			cout<<"Copying MatFact is not yet implemented, though the Matrix part will be copied normally\n";
			return *this;
		}
		
		inline void New() { fact = false; }
		
		inline int  Row() const { return Matrix::Row(); }
		inline int  Col() const { return Matrix::Col(); }
		
		void Fact();
		void Solve( Vector& X, const Vector& B, bool TRANS=false );
		void Solve( Matrix& X, const Matrix& B, bool TRANS=false );
};

typedef Array1D< Matrix > JagMatrix3D;
typedef Array1D< Vector > JagVector2D;

/// specialized versions of the Clear function
template< > inline void Array1D< Array1D<int> >::Clear( ) { for( int i=0; i<n; i++ ) v[i].Clear(); }
template< > inline void Array1D< Array1D<double> >::Clear( ) { for( int i=0; i<n; i++ ) v[i].Clear(); }
template< > inline void Array1D< Vector >::Clear( ) { for( int i=0; i<n; i++ ) v[i].Clear(); }
template< > inline void Array1D< Array2D<int> >::Clear( ) { for( int i=0; i<n; i++ ) v[i].Clear(); }
template< > inline void Array1D< Array2D<double> >::Clear( ) { for( int i=0; i<n; i++ ) v[i].Clear(); }
template< > inline void Array1D< Matrix >::Clear( ) { for( int i=0; i<n; i++ ) v[i].Clear(); }

/// Member functions and Operators for __scal_vec_trans_rng
template<> inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng( const Vector& v )
	: rng(0, v.Size() ), vec( v ), alpha( 1.0 ), tr( NoTrans ) { }

template<> inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng( const Vector& v, double a )
	: rng(0, v.Size() ), vec( v ), alpha( a ),   tr( NoTrans ) { }

template<> inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng( const Vector& v, double a, enum cspblas_Trans t )
	: rng(0, v.Size() ), vec( v ), alpha( a ),   tr( t ) { }

template< class VT > inline __scal_vec_trans_rng<VT>::__scal_vec_trans_rng( __scal_vec_trans<VT> v )
	: rng( 0, v.vec.Row(), 0, v.vec.Col() ), vec( v.vec ), alpha( v.alpha ), tr( v.tr ) { }

template< >          inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng( __scal_vec_trans<Vector> v )
	: rng( 0, v.vec.Size() ), vec( v.vec ), alpha( v.alpha ), tr( v.tr ) { }

template< class VT > inline __op_mul_vec_rng< VT, Vector > __scal_vec_trans_rng<VT>::operator*( __scal_vec_trans_rng<Vector> v )
	{ return __op_mul_vec_rng< VT, Vector >( *this, v ); }

template< class VT > inline __op_mul_vec_rng< VT, Matrix > __scal_vec_trans_rng<VT>::operator*( __scal_vec_trans_rng<Matrix> v )
	{ return __op_mul_vec_rng< VT, Matrix >( *this, v ); }

template< class VT > inline __op_mul_vec_rng< VT, Vector > __scal_vec_trans_rng<VT>::operator*( const Vector& v )
	{ return __op_mul_vec_rng< VT, Vector >( *this, __scal_vec_trans_rng<Vector>( v ) ); }

template< class VT > inline __op_mul_vec_rng< VT, Matrix > __scal_vec_trans_rng<VT>::operator*( const Matrix& v )
	{ return __op_mul_vec_rng< VT, Matrix >( *this, __scal_vec_trans_rng<Matrix>( v ) ); }


/// Member functions and Operators for __op_mul_vec_rng
template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>   __op_mul_vec_rng<MT,VT>::operator+( __scal_vec_trans_rng<VT> v )
	{ return __op_mul_vec_plus_vec_rng<MT,VT>( *this, v ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>   __op_mul_vec_rng<MT,VT>::operator+( const VT& v )
		{ return __op_mul_vec_plus_vec_rng<MT,VT>( *this, __scal_vec_trans_rng<VT>( v ) ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>   __op_mul_vec_rng<MT,VT>::operator+( __scal_vec_trans<VT> v )
	{ return __op_mul_vec_plus_vec_rng<MT,VT>( *this, __scal_vec_trans_rng<VT>( v ) ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>   __op_mul_vec_rng<MT,VT>::operator-( __scal_vec_trans_rng<VT> v )
	{ return __op_mul_vec_plus_vec_rng<MT,VT>( *this, __scal_vec_trans_rng<VT>( v, -v.alpha ) ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>   __op_mul_vec_rng<MT,VT>::operator-( const VT& v )
	{ return __op_mul_vec_plus_vec_rng<MT,VT>( *this, __scal_vec_trans_rng<VT>( v, -1.0 ) ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>   __op_mul_vec_rng<MT,VT>::operator-( __scal_vec_trans<VT> v )
	{ return __op_mul_vec_plus_vec_rng<MT,VT>( *this, __scal_vec_trans_rng<VT>( v.vec, -v.alpha, v.tr ) ); }


/// Member functions and Operators for __scal_vec_trans
template< class VT > inline __op_mul_vec< VT, Vector >     __scal_vec_trans<VT>::operator*( const Vector& v )
	{ return __op_mul_vec< VT, Vector >( *this, __scal_vec_trans<Vector>( v ) ); }

template< class VT > inline __op_mul_vec< VT, Matrix >     __scal_vec_trans<VT>::operator*( const Matrix& v )
	{ return __op_mul_vec< VT, Matrix >( *this, __scal_vec_trans<Matrix>( v ) ); }

template< class VT > inline __op_mul_vec_rng< VT, Vector > __scal_vec_trans<VT>::operator*( __scal_vec_trans_rng<Vector> v )
	{ return __op_mul_vec_rng< VT, Vector >( __scal_vec_trans_rng<VT>( *this ), v ); }

template< class VT > inline __op_mul_vec_rng< VT, Matrix > __scal_vec_trans<VT>::operator*( __scal_vec_trans_rng<Matrix> v )
	{ return __op_mul_vec_rng< VT, Matrix >( __scal_vec_trans_rng<VT>( *this ), v ); }


/// Member functions and Operators for __op_mul_vec
template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>  __op_mul_vec<MT,VT>::operator+( __scal_vec_trans_rng<VT> v )
	{ return __op_mul_vec_plus_vec_rng<MT,VT>( __op_mul_vec_rng<MT,VT>( this->op, this->vecA ), v ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT,VT>  __op_mul_vec<MT,VT>::operator-( __scal_vec_trans_rng<VT> v )
	{ return __op_mul_vec_plus_vec_rng<MT,VT>( __op_mul_vec_rng<MT,VT>( this->op, this->vecA ), __scal_vec_trans_rng<VT>( v.vec, -v.alpha, v.tr ) ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT,VT>      __op_mul_vec<MT,VT>::operator+( __scal_vec_trans<VT> v )
	{ return __op_mul_vec_plus_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT,VT>      __op_mul_vec<MT,VT>::operator-( __scal_vec_trans<VT> v )
	{ return __op_mul_vec_plus_vec<MT,VT>( *this, __scal_vec_trans<VT>( v.vec, -v.alpha, v.tr ) ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT,VT>      __op_mul_vec<MT,VT>::operator+( const VT& v )
	{ return __op_mul_vec_plus_vec<MT,VT>( *this, __scal_vec_trans<VT>( v ) ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT,VT>      __op_mul_vec<MT,VT>::operator-( const VT& v )
	{ return __op_mul_vec_plus_vec<MT,VT>( *this, __scal_vec_trans<VT>( v, -1.0 ) ); }


/// The global operators
inline __scal_vec_trans < Vector >                      operator*( double a, const Vector& v )
	{ return __scal_vec_trans<Vector>( v, a ); }

inline __scal_vec_trans < Matrix >                      operator*( double a, const Matrix& v )
	{ return __scal_vec_trans<Matrix>( v, a ); }

inline __scal_vec_trans < SpMatrix >                    operator*( double a, const SpMatrix& v )
	{ return __scal_vec_trans<SpMatrix>( v, a ); }

template< class VT > inline __scal_vec_trans < VT  >    operator*( double a, __scal_vec_trans< VT > v )
	{ return __scal_vec_trans<VT>( v.vec, a * v.alpha, v.tr ); }

template< class VT > inline __scal_vec_trans_rng < VT > operator*( double a, __scal_vec_trans_rng< VT > v )
	{ return __scal_vec_trans_rng<VT>( v.vec, v, a * v.alpha, v.tr ); }

inline __scal_vec_trans < Vector >                      operator-( const Vector& v )
	{ return __scal_vec_trans<Vector>( v, -1.0 ); }

inline __scal_vec_trans < Matrix >                      operator-( const Matrix& v )
	{ return __scal_vec_trans<Matrix>( v, -1.0 ); }

inline __scal_vec_trans < SpMatrix >                    operator-( const SpMatrix& v )
	{ return __scal_vec_trans<SpMatrix>( v, -1.0 ); }

// Implementation of Matrix

inline void Matrix::AX( double* out, const double* in, double alpha, bool trans ) const
{
	mmx( trans ? Trans : NoTrans, out, in, alpha );
}

inline void Matrix::AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const
{
	mmxpy( trans ? Trans : NoTrans, out, in, alpha, Y, beta );
}

// End of implementation of Matrix

// Implementation of Vector

inline Vector& Vector::operator=( const __scal_vec_trans<Vector> )
{
	std::cout<<"__scal_vec\n"; return *this;
}

inline Vector& Vector::operator=( const __op_mul_vec<Matrix,Vector> R )
{
	R.op.vec.mmx( R.op.tr, this->v, R.vecA.vec.v, R.op.alpha );
// 	if( R.vecA.alpha != 1.0 || R.vecA.tr != NoTrans ) std::cout<<"Dmmx err";
// 	std::cout<<"__op_mul_vec<Matrix,Vector>\n";
	return *this;
}

inline Vector& Vector::operator=( const __op_mul_vec_plus_vec<Matrix,Vector> R )
{
	R.op.vec.mmxpy( R.op.tr, this->v, R.vecA.vec.v, R.op.alpha, R.vecB.vec.v, R.vecB.alpha );
// 	if( R.vecA.alpha != 1.0 || R.vecA.tr != NoTrans || R.vecB.tr != NoTrans ) std::cout<<"Dmmxpy err";
// 	std::cout<<"__op_mul_vec_plus_vec<Matrix,Vector>\n";
	return *this;
}

inline Vector& Vector::operator=( const __scal_vec_trans_rng<Vector> )
{
	std::cout<<"__scal_vec_rng\n"; return *this;
	return *this;
}

inline Vector& Vector::operator=( const __op_mul_vec_rng<Matrix,Vector> R )
{
	cblas_mmx( R.op.tr, R.op.i2 - R.op.i1, R.op.j2 - R.op.j1, &R.op.vec.m[ R.op.i1 + R.op.vec.r * R.op.j1 ],
	           R.op.vec.r, this->v, &R.vecA.vec.v[ R.vecA.i1 ], R.op.alpha );

// 	if( R.vecA.alpha != 1.0 || R.vecA.tr != NoTrans ) std::cout<<"Dmmx err";
// 	std::cout<<"__op_mul_vec_rngMatrix,Vector>\n";
	return *this;
}

inline Vector& Vector::operator=( const __op_mul_vec_plus_vec_rng<Matrix,Vector> R )
{
	cblas_mmxpy( R.op.tr, R.op.i2 - R.op.i1, R.op.j2 - R.op.j1, &R.op.vec.m[ R.op.i1 + R.op.vec.r * R.op.j1 ], R.op.vec.r,
	             this->v, &R.vecA.vec.v[ R.vecA.i1 ], R.op.alpha, &R.vecB.vec.v[ R.vecB.i1 ], R.vecB.alpha );
// 	std::cout<<"__op_mul_vec_plus_vec_rng<Matrix,Vector>\n";
	return *this;
}

// End of implementation of Vector

// Implementation of Matrix

inline Matrix& Matrix::operator=( const __scal_vec_trans<Matrix> )
{
	std::cout<<"__scal_vec\n"; return *this;
	return *this;
}

inline Matrix& Matrix::operator=( const __op_mul_vec<Matrix,Matrix> R )
{
	R.op.vec.mmxm( R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecA.vec.c );
// 	if( R.vecA.alpha != 1.0 || R.vecA.tr != NoTrans ) std::cout<<"Dmmxm err";
// 	std::cout<<"__op_mul_vec<Matrix,Matrix>\n";
	return *this;
}

inline Matrix& Matrix::operator=( const __op_mul_vec_plus_vec<Matrix,Matrix> R )
{
	R.op.vec.mmxmpym( R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecB.vec.m, R.vecB.vec.r, R.vecB.alpha, R.vecA.vec.c );
// 	if( R.vecA.alpha != 1.0 || R.vecA.tr != NoTrans || R.vecB.tr != NoTrans ) std::cout<<"Dmmxmpym err";
	std::cout<<"__op_mul_vec_plus_vec<Matrix,Matrix>\n";
	return *this;
}

inline Matrix& Matrix::operator=( const __scal_vec_trans_rng<Matrix> )
{
	std::cout<<"__scal_vec_rng\n"; return *this;
	return *this;
}

inline Matrix& Matrix::operator=( const __op_mul_vec_rng<Matrix,Matrix> R )
{
	cblas_mmxm( R.op.tr, R.op.i2 - R.op.i1, R.op.j2 - R.op.j1, &R.op.vec.m[ R.op.i1 + R.op.vec.r * R.op.j1 ], R.op.vec.r,
	           this->m, this->r, &R.vecA.vec.m[ R.vecA.i1 + R.vecA.vec.r * R.vecA.j1 ], R.vecA.vec.r, R.op.alpha, R.vecA.j2 - R.vecA.j1 );
// 	std::cout<<"__op_mul_vec_rngMatrix,Matrix>\n";
	return *this;
}

inline Matrix& Matrix::operator=( const __op_mul_vec_plus_vec_rng<Matrix,Matrix> R )
{
	cblas_mmxmpym( R.op.tr, R.op.i2 - R.op.i1, R.op.j2 - R.op.j1, &R.op.vec.m[ R.op.i1 + R.op.vec.r * R.op.j1 ], R.op.vec.r,
	              this->m, this->r, // out
	              &R.vecA.vec.m[ R.vecA.i1 + R.vecA.vec.r * R.vecA.j1 ], R.vecA.vec.r, R.op.alpha, // in1
	              &R.vecB.vec.m[ R.vecB.i1 + R.vecB.vec.r * R.vecB.j1 ], R.vecB.vec.r, R.vecB.alpha, // in2
	              R.vecA.j2 - R.vecA.j1 ); // nrhs
// 	std::cout<<"__op_mul_vec_plus_vec_rng<Matrix,Matrix>\n";
	return *this;
}

// End of implementation of Matrix

#endif // PDDESYS_H

#endif
