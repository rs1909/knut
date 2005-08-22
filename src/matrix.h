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

	Array1D() { n = 0; v = 0; }

	Array1D( int i )
	{
		n = i;
		v = new T[i];
		Clear();
	}

	// specially for JagMatrix2D i.e. Array1D< Vector >, which is indexed as (i)(j)
	Array1D( int i, int j ) 
	{
		n = i;
		v = new T[i];
		for( int r=0; r<i; r++ ) v[r].Init(j);
	}

	// specially for JagMatrix3D i.e. Array1D< Matrix >, which is indexed as (i)(j,k)
	Array1D( int i, int j, int k )
	{
		n = i;
		v = new T[i+1];
		for( int r=0; r<i; r++ ) v[r].Init(j,k);
	}

	Array1D( const Array1D<T> &V )
	{
		n = V.n; 
		v = new T[V.n];
		for( int i=0; i<V.n; i++ ) v[i]=V.v[i]; 
	}

  ~Array1D()
	{
		delete[] v; 
	}
	
	void Init( int i )
	{
		delete v;
		n = i; 
		v = new T[i];
		Clear(); 
	}

	void Clear()
	{
		for(int j=0; j<n; j++) v[j]=T(); // it was T(0), but for built in types it must be the same
	}

	int Size() const { return n; }
	
  	Array1D<T>& operator=( const Array1D<T>& V )
	{
#ifdef DEBUG
		if( n != V.n )
		{
			std::cout<<"Array1D::operator=(): incompatible sizes\n";
			throw(-1);
		}
#endif
		for( int i=0; i < n; i++ ) v[i] = V.v[i];
		return *this;
	}

	T& operator()( int i )
	{

#ifdef DEBUG

		if( i>=n ){ std::cout<<"Array1D::bound11& "<<n<<", "<<i<<"\n"; throw(12); }

#endif

		return v[i];
	}

	T operator()( int i ) const 
	{

#ifdef DEBUG

		if( i>=n ){ std::cout<<"vec_bound11_\n"; throw(12); }

#endif

		return v[i]; 
	}

};


template<class T>
class Array2D{

 protected:

	T* m;
	int r,c;

 public:
	Array2D( ) { r = 0; c = 0; m = 0; }

	Array2D( int _r, int _c ) 
		: r(_r), c(_c) 
	{
		m = new T[r*c+1]; 
		Clear(); 
	}

	Array2D( const Array2D<T>& M ) 
		: r(M.r), c(M.c)
	{
		m = new T[r*c+1];
		for( int i = 0; i < r*c; i++ ) m[i] = M.m[i]; 
	}

	virtual ~Array2D(){ delete []m; }

	void Init( int _r, int _c )
	{
		delete m;
		r = _r; 
		c = _c;
		m = new T[r*c+1]; 
		Clear(); 
	}

	void Clear()
	{
		for( int i = 0; i < r*c; i++ ) m[i] = T(0); 
	}


	Array2D<T>& operator= ( const Array2D<T>& M )
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

	T& operator()( const int i, const int j )
	{

#ifdef DEBUG 

		if((i>=r)||(j>=c)){ cout<<"bound& "<<r<<", "<<c<<",-"<<i<<", "<<j<<"\n"; throw(1); } 
		if((i<0)||(j<0)){ cout<<"lbound& "<<i<<", "<<j<<"\n"; throw(1); }

#endif
		return m[i + r*j]; 
	}

	T operator()( const int i, const int j ) const 
	{

#ifdef DEBUG 

		if((i>=r)||(j>=c)){ cout<<"bound_\n"; throw(1); }
		if((i<0)||(j<0)){ cout<<"lbound_\n"; throw(1); }

#endif

		return m[i + r*j]; 
	}

	T& operator()( int i )
	{

#ifdef DEBUG

		if( i>=r*c ){ cout<<"bound11&\n"; throw(1); }
		if((r!=1)&&(c!=1)){ cout<<"H&\n"; throw(1); }

#endif

		return m[i];
	}

	T operator()( int i ) const
	{

#ifdef DEBUG

		if( i>=r*c ){ cout<<"bound11_\n"; throw(1); }
		if((r!=1)&&(c!=1)){ cout<<"H&\n"; throw(1); }

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
	Array3D( ) { d1 = d2 = d3 = 0; m = 0; }

	Array3D( int _d1, int _d2, int _d3 ) 
		: d1(_d1), d2(_d2), d3(_d3) 
	{
		m = new T[d1*d2*d3+1]; 
		Clear(); 
	}

	Array3D( const Array3D<T>& M ) 
		: d1(M.d1), d2(M.d2), d3(M.d3)
	{
		m = new T[d1*d2*d3+1];
		for( int i = 0; i < d1*d2*d3; i++ ) m[i] = M.m[i]; 
	}

	virtual ~Array3D(){ delete []m; }

	void Init( int _d1, int _d2, int _d3 ) 
	{ 
		delete m;
		d1 = _d1;
		d2 = _d2;
		d3 = _d3;
		m = new T[d1*d2*d3+1]; 
		Clear(); 
	}
	void Clear()
	{
		for( int i = 0; i < d1*d2*d3; i++ ) m[i] = 0; 
	}

	Array3D<T>& operator= ( const Array3D<T>& M )
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

	T& operator()( const int i, const int j, const int k )
	{

#ifdef DEBUG 

		if((i>=d1)||(j>=d2)||(k>=d3)) cout<<"bound&\n"; 
		if((i<0)||(j<0)||(k<0)) cout<<"lbound&\n";

#endif
		return m[i + d1*(j + d2*k)];
	}

	T operator()( const int i, const int j, const int k ) const 
	{

#ifdef DEBUG 

		if((i>=d1)||(j>=d2)||(k>=d3)) cout<<"bound&\n"; 
		if((i<0)||(j<0)||(k<0)) cout<<"lbound&\n";

#endif
		return m[i + d1*(j + d2*k)];
	}
	
};

// pre declaration

class SpMatrix;
class Matrix;
class Vector;

template< class MT > class __vec_trans;
template< class MT > class __scal_vec;
template< class MT > class __scal_vec_trans;
template< class MT, class VT > class __op_mul_vec;
template< class MT, class VT > class __op_trans_mul_vec;
template< class MT, class VT > class __scal_op_mul_vec;
template< class MT, class VT > class __scal_op_trans_mul_vec;
template< class MT, class VT > class __op_mul_vec_plus_vec;
template< class MT, class VT > class __op_trans_mul_vec_plus_vec;
template< class MT, class VT > class __scal_op_mul_vec_plus_vec;
template< class MT, class VT > class __scal_op_trans_mul_vec_plus_vec;
template< class MT, class VT > class __op_mul_vec_plus_scal_vec;
template< class MT, class VT > class __op_trans_mul_vec_plus_scal_vec;
template< class MT, class VT > class __scal_op_mul_vec_plus_scal_vec;
template< class MT, class VT > class __scal_op_trans_mul_vec_plus_scal_vec;

// First level

template< class VT > class __vec_trans
{
 public:
	__vec_trans( VT& m ) : vec( m ) { }
	
	__scal_vec_trans<VT>          operator*( double d );
	__op_trans_mul_vec<VT,Vector> operator*( Vector& v );
	__op_trans_mul_vec<VT,Matrix> operator*( Matrix& v );
	
	VT& vec;
};

template< class VT > class __scal_vec
{
 public:
	__scal_vec( double a, VT& v ) : alpha( a ), vec( v ) { }
	
	__scal_op_mul_vec<VT,Vector> operator*( Vector& v );
	__scal_op_mul_vec<VT,Matrix> operator*( Matrix& v );

	double alpha;
	VT& vec;
};

template< class VT > class __scal_vec_trans
{
 public:
	__scal_vec_trans( double a, __vec_trans<VT> m ) : alpha( a ), vec( m.vec ) { }
	
	__scal_op_trans_mul_vec<VT,Vector> operator*( Vector& v );
	__scal_op_trans_mul_vec<VT,Matrix> operator*( Matrix& v );
	
	double alpha;
	VT& vec;
};

// Second level
// Multiplied by VT or __scal_vec

template< class MT, class VT > class __op_mul_vec
{
 public:
	__op_mul_vec( MT& m, VT& v ) : op( m ), vecA( v ) { }
	
	__op_mul_vec_plus_vec<MT,VT>      operator+( VT& v );
	__op_mul_vec_plus_scal_vec<MT,VT> operator+( __scal_vec<VT> v );
	
	MT& op;
	VT& vecA;
};

template< class MT, class VT > class __op_trans_mul_vec
{
 public:
	__op_trans_mul_vec( __vec_trans<MT> m, VT& v ) : op( m.vec ), vecA( v ) { }
	
	__op_trans_mul_vec_plus_vec<MT,VT>      operator+( VT& v );
	__op_trans_mul_vec_plus_scal_vec<MT,VT> operator+( __scal_vec<VT> v );
	
	MT& op;
	VT& vecA;
};

template< class MT, class VT > class __scal_op_mul_vec
{
 public:
	__scal_op_mul_vec( __scal_vec<MT> m, VT& v ) : alpha( m.alpha ), op( m.vec ), vecA( v ) { }

	__scal_op_mul_vec_plus_vec<MT,VT>      operator+( VT& v );
	__scal_op_mul_vec_plus_scal_vec<MT,VT> operator+( __scal_vec<VT> v );

	double alpha;
	MT& op;
	VT& vecA;
};

template< class MT, class VT > class __scal_op_trans_mul_vec
{
 public:
	__scal_op_trans_mul_vec( __scal_vec_trans<MT> m, VT& v ) : alpha( m.alpha ), op( m.vec ), vecA( v ) { }
	
	__scal_op_trans_mul_vec_plus_vec<MT,VT>      operator+( VT& v );
	__scal_op_trans_mul_vec_plus_scal_vec<MT,VT> operator+( __scal_vec<VT> v );
	
	double alpha;
	MT& op;
	VT& vecA;
};

// Third level
// No operation. These are the arguments of the respective operator=

template< class MT, class VT > class __op_mul_vec_plus_vec : public __op_mul_vec<MT,VT>
{
 public:
	__op_mul_vec_plus_vec( __op_mul_vec<MT,VT> m, VT& v ) : __op_mul_vec<MT,VT>( m ), vecB( v ) { }

	VT& vecB;
};

template< class MT, class VT > class __op_trans_mul_vec_plus_vec : public __op_trans_mul_vec<MT,VT>
{
 public:
	__op_trans_mul_vec_plus_vec( __op_trans_mul_vec<MT,VT> m, VT& v ) : __op_trans_mul_vec<MT,VT>( m ), vecB( v ) { }
	
	VT& vecB;
};

template< class MT, class VT > class __scal_op_mul_vec_plus_vec : public __scal_op_mul_vec<MT,VT>
{
 public:
	__scal_op_mul_vec_plus_vec( __scal_op_mul_vec<MT,VT> m, VT& v ) : __scal_op_mul_vec<MT,VT>( m ), vecB( v ) { }
	
	VT& vecB;
};

template< class MT, class VT > class __scal_op_trans_mul_vec_plus_vec : public __scal_op_trans_mul_vec<MT,VT>
{
 public:
	__scal_op_trans_mul_vec_plus_vec( __scal_op_trans_mul_vec<MT,VT> m, VT& v ) : __scal_op_trans_mul_vec<MT,VT>( m ), vecB( v ) { }
	
	VT& vecB;
};

template< class MT, class VT > class __op_mul_vec_plus_scal_vec : public __op_mul_vec<MT,VT>
{
 public:
	__op_mul_vec_plus_scal_vec( __op_mul_vec<MT,VT> m, __scal_vec<VT>& v ) : __op_mul_vec<MT,VT>( m ), beta( v.alpha ), vecB( v.vec ) { }
	
	double beta;
	VT& vecB;
};

template< class MT, class VT > class __op_trans_mul_vec_plus_scal_vec : public __op_trans_mul_vec<MT,VT>
{
 public:
	__op_trans_mul_vec_plus_scal_vec( __op_trans_mul_vec<MT,VT> m, __scal_vec<VT>& v ) : __op_trans_mul_vec<MT,VT>( m ), beta( v.alpha ), vecB( v.vec ) { }
	
	double beta;
	VT& vecB;
};

template< class MT, class VT > class __scal_op_mul_vec_plus_scal_vec : public __scal_op_mul_vec<MT,VT>
{
 public:
	__scal_op_mul_vec_plus_scal_vec( __scal_op_mul_vec<MT,VT> m, __scal_vec<VT>& v ) : __scal_op_mul_vec<MT,VT>( m ), beta( v.alpha ), vecB( v.vec ) { }
	
	double beta;
	VT& vecB;
};

template< class MT, class VT > class __scal_op_trans_mul_vec_plus_scal_vec : public __scal_op_trans_mul_vec<MT,VT>
{
 public:
	__scal_op_trans_mul_vec_plus_scal_vec( __scal_op_trans_mul_vec<MT,VT> m, __scal_vec<VT>& v ) : __scal_op_trans_mul_vec<MT,VT>( m ), beta( v.alpha ), vecB( v.vec ) { }
	
	double beta;
	VT& vecB;
};

// the global operators
inline __scal_vec<Vector>         operator*( double a, Vector& v ) { return __scal_vec<Vector>( a, v ); }
inline __scal_vec<Matrix>         operator*( double a, Matrix& v ) { return __scal_vec<Matrix>( a, v ); }
inline __scal_vec<SpMatrix>       operator*( double a, SpMatrix& v ) { return __scal_vec<SpMatrix>( a, v ); }
inline __scal_vec_trans<Vector>   operator*( double a, __vec_trans<Vector> v ) { return __scal_vec_trans<Vector>( a, v ); }
inline __scal_vec_trans<Matrix>   operator*( double a, __vec_trans<Matrix> v ) { return __scal_vec_trans<Matrix>( a, v ); }
inline __scal_vec_trans<SpMatrix> operator*( double a, __vec_trans<SpMatrix> v ) { return __scal_vec_trans<SpMatrix>( a, v ); }

template< class VT > inline __scal_vec_trans<VT> __vec_trans<VT>::operator*( double d )
{ return __scal_vec_trans<VT>( d, *this ); }

template< class VT > inline __op_trans_mul_vec<VT,Vector> __vec_trans<VT>::operator*( Vector& v ) 
{ return __op_trans_mul_vec<VT,Vector>( *this, v ); }

template< class VT > inline __op_trans_mul_vec<VT,Matrix> __vec_trans<VT>::operator*( Matrix& v ) 
{ return __op_trans_mul_vec<VT,Matrix>( *this, v ); }

template< class VT > inline __scal_op_mul_vec<VT,Vector> __scal_vec<VT>::operator*( Vector& v ) 
{ return __scal_op_mul_vec<VT,Vector>( *this, v ); }

template< class VT > inline __scal_op_mul_vec<VT,Matrix> __scal_vec<VT>::operator*( Matrix& v ) 
{ return __scal_op_mul_vec<VT,Matrix>( *this, v ); }

template< class VT > inline __scal_op_trans_mul_vec<VT,Vector> __scal_vec_trans<VT>::operator*( Vector& v ) 
{ return __scal_op_trans_mul_vec<VT,Vector>( *this, v ); }

template< class VT > inline __scal_op_trans_mul_vec<VT,Matrix> __scal_vec_trans<VT>::operator*( Matrix& v ) 
{ return __scal_op_trans_mul_vec<VT,Matrix>( *this, v ); }

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT,VT> __op_mul_vec<MT,VT>::operator+( VT& v )
{ return __op_mul_vec_plus_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __op_mul_vec_plus_scal_vec<MT,VT> __op_mul_vec<MT,VT>::operator+( __scal_vec<VT> v ) 
{ return __op_mul_vec_plus_scal_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __op_trans_mul_vec_plus_vec<MT,VT> __op_trans_mul_vec<MT,VT>::operator+( VT& v ) 
{ return __op_trans_mul_vec_plus_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __op_trans_mul_vec_plus_scal_vec<MT,VT> __op_trans_mul_vec<MT,VT>::operator+( __scal_vec<VT> v ) 
{ return __op_trans_mul_vec_plus_scal_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __scal_op_mul_vec_plus_vec<MT,VT> __scal_op_mul_vec<MT,VT>::operator+( VT& v )
{ return __scal_op_mul_vec_plus_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __scal_op_mul_vec_plus_scal_vec<MT,VT> __scal_op_mul_vec<MT,VT>::operator+( __scal_vec<VT> v )
{ return __scal_op_mul_vec_plus_scal_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __scal_op_trans_mul_vec_plus_vec<MT,VT> __scal_op_trans_mul_vec<MT,VT>::operator+( VT& v ) 
{ return __scal_op_trans_mul_vec_plus_vec<MT,VT>( *this, v ); }

template< class MT, class VT > inline __scal_op_trans_mul_vec_plus_scal_vec<MT,VT> __scal_op_trans_mul_vec<MT,VT>::operator+( __scal_vec<VT> v ) 
{ return __scal_op_trans_mul_vec_plus_scal_vec<MT,VT>( *this, v ); }

class Vector : public Array1D<double> 
{

public:

	Vector() { }

	Vector( int i ) : Array1D<double>( i ) { }

	Vector( const Vector& V ) : Array1D<double>( V ) { }

	~Vector() { }

	double *Pointer(){ return v; }
	
  	Vector& operator+=( const Vector& V )
	{
		if( n != V.n ){ std::cout<<"Vector::operator+=(): incompatible sizes\n"; }
		for( int i=0; i < n; i++ ) v[i] += V.v[i];
		return *this;
	}

	Vector& operator-=( const Vector& V )
	{
		if( n != V.n ){ std::cout<<"Vector::operator-=(): incompatible sizes\n"; }
		for( int i=0; i < n; i++ ) v[i] -= V.v[i];
		return *this;
  }

	Vector& operator/=( double div )
	{
		for( int i=0; i < n; i++ ) v[i] /= div;
		return *this;
	}

	Vector& operator*=( double mul )
	{
		for( int i=0; i < n; i++ ) v[i] *= mul;
		return *this;
	}

	double operator*( const Vector& V ) const
	{
		register double sum=0.0;
		if( n != V.n ){ std::cout<<"Vector::operator*(): incompatible sizes\n"; }
		for( int i=0; i < n; i++ ) sum += v[i]*V.v[i];
		return sum;
	}

	void Print() const
	{
		for(int j=0; j<n; j++) std::cout<<v[j]<<'\t'; std::cout<<'\n'; 
	}

	// Matrix operations
	
	Vector& operator=( const __scal_vec<Vector> );
	Vector& operator=( const __op_mul_vec<Matrix,Vector> );
	Vector& operator=( const __op_trans_mul_vec<Matrix,Vector> );
	Vector& operator=( const __scal_op_mul_vec<Matrix,Vector> );
	Vector& operator=( const __scal_op_trans_mul_vec<Matrix,Vector> );
	Vector& operator=( const __op_mul_vec_plus_vec<Matrix,Vector> );
	Vector& operator=( const __op_trans_mul_vec_plus_vec<Matrix,Vector> );
	Vector& operator=( const __scal_op_mul_vec_plus_vec<Matrix,Vector> );
	Vector& operator=( const __scal_op_trans_mul_vec_plus_vec<Matrix,Vector> );
	Vector& operator=( const __op_mul_vec_plus_scal_vec<Matrix,Vector> );
	Vector& operator=( const __op_trans_mul_vec_plus_scal_vec<Matrix,Vector> );
	Vector& operator=( const __scal_op_mul_vec_plus_scal_vec<Matrix,Vector> );
	Vector& operator=( const __scal_op_trans_mul_vec_plus_scal_vec<Matrix,Vector> );

	// Sparse matrix operations
	
	Vector& operator=( const __op_mul_vec<SpMatrix,Vector> );
	Vector& operator=( const __op_trans_mul_vec<SpMatrix,Vector> );
	Vector& operator=( const __scal_op_mul_vec<SpMatrix,Vector> );
	Vector& operator=( const __scal_op_trans_mul_vec<SpMatrix,Vector> );
	Vector& operator=( const __op_mul_vec_plus_vec<SpMatrix,Vector> );
	Vector& operator=( const __op_trans_mul_vec_plus_vec<SpMatrix,Vector> );
	Vector& operator=( const __scal_op_mul_vec_plus_vec<SpMatrix,Vector> );
	Vector& operator=( const __scal_op_trans_mul_vec_plus_vec<SpMatrix,Vector> );
	Vector& operator=( const __op_mul_vec_plus_scal_vec<SpMatrix,Vector> );
	Vector& operator=( const __op_trans_mul_vec_plus_scal_vec<SpMatrix,Vector> );
	Vector& operator=( const __scal_op_mul_vec_plus_scal_vec<SpMatrix,Vector> );
	Vector& operator=( const __scal_op_trans_mul_vec_plus_scal_vec<SpMatrix,Vector> );

	friend class SpMatrix;
	friend class SpFact;
	friend class Matrix;
	friend class MatFact;

};

class Matrix : public Array2D<double>
{
 public:

	Matrix() { }

	Matrix( int i, int j ) : Array2D<double>( i, j ) { }

	Matrix( const Matrix& M ) : Array2D<double>( M ) { }

	virtual ~Matrix() { }

	int Row() const { return r; }
	int Col() const { return c; }
	int Size() const { if((r==1)||(c==1)){return r*c;}else{cout<<"Hs\n"; return 0;} }
	void Copy( int i, int j, Matrix& mat );
	
	void CopyCol( int k, Vector& C )
	{
		if((r==C.Size())&&(k<c))
		{
			for( int i=0; i<r; i++ ) m[i+k*r] = C(i);
		}else
		{
			cout<<"CopyCol:err"; 
		}
	}
	
	void CopyRow( int k, Vector& C )
	{
		if((c==C.Size())&&(k<r))
		{
			for( int i=0; i<c; i++ ) m[k+i*r] = C(i);
		}else
		{
			cout<<"CopyRow:err";
		}
	}

#ifndef PDDESYS_H
	void mmx( enum cspblas_Trans trans, double* out, const double* in, double alpha ) const
	{
		cblas_mmx( trans, r, c,  m, r, out, in, alpha );
	}
	void mmxpy( enum cspblas_Trans trans, double* out, const double* in, double alpha, const double* Y, double beta ) const
	{
		cblas_mmxpy( trans, r, c, m, r, out, in, alpha, Y, beta );
	}
	void mmxm( enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha, int nrhs ) const
	{
		cblas_mmxm( trans, r, c, m, r, out, ldout, in, ldin, alpha, nrhs );
	}
	void mmxmpym( enum cspblas_Trans trans, double* out, int ldout, 
	              const double* in, int ldin, double alpha, const double* Y, int ldY, double beta, int nrhs ) const
	{
		cblas_mmxmpym( trans, r, c, m, r, out, ldout, in, ldin, alpha, Y, ldY, beta, nrhs );
	}

	void AX( double* out, const double* in, double alpha, bool trans ) const;
	void AX( Matrix& out, const Matrix& in, double alpha, bool trans ) const;
	void AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const;
	void AXpY( Matrix& out, const Matrix& in, const Matrix& Y, double alpha, double beta, bool trans ) const;
	
	// with Vectors
	void AX( Vector& X, const Vector& B, double alpha, bool trans ) const;
	void AXpY( Vector& out, const Vector& X, const Vector& Y, double alpha, double beta, bool trans ) const;
	
	// other variants
	void AX( Vector& X, const Vector& B ) const { this->AX( X, B, 1.0, false ); }

	void Eigval( Vector& re, Vector& im );
	void Eigval( Vector& re, Vector& im, Matrix& lev, Matrix& rev );
	void StrPlot( GnuPlot& pl );
#endif

	/* operators */
	Matrix& operator*=( double mul )
	{
		for( int i=0; i < r*c; i++ ) m[i] *= mul;
		return *this;
	}
	
	__op_mul_vec<Matrix,Vector> operator*( Vector& v ) { return __op_mul_vec<Matrix,Vector>( *this, v ); }
	__op_mul_vec<Matrix,Matrix> operator*( Matrix& v ) { return __op_mul_vec<Matrix,Matrix>( *this, v ); }
	__vec_trans<Matrix>         operator!( ) { return __vec_trans<Matrix>( *this ); }
	
	// Sparse matrix operations
	
	Matrix& operator=( const __scal_vec<Matrix> );
	Matrix& operator=( const __op_mul_vec<Matrix,Matrix> );
	Matrix& operator=( const __op_trans_mul_vec<Matrix,Matrix> );
	Matrix& operator=( const __scal_op_mul_vec<Matrix,Matrix> );
	Matrix& operator=( const __scal_op_trans_mul_vec<Matrix,Matrix> );
	Matrix& operator=( const __op_mul_vec_plus_vec<Matrix,Matrix> );
	Matrix& operator=( const __op_trans_mul_vec_plus_vec<Matrix,Matrix> );
	Matrix& operator=( const __scal_op_mul_vec_plus_vec<Matrix,Matrix> );
	Matrix& operator=( const __scal_op_trans_mul_vec_plus_vec<Matrix,Matrix> );
	Matrix& operator=( const __op_mul_vec_plus_scal_vec<Matrix,Matrix> );
	Matrix& operator=( const __op_trans_mul_vec_plus_scal_vec<Matrix,Matrix> );
	Matrix& operator=( const __scal_op_mul_vec_plus_scal_vec<Matrix,Matrix> );
	Matrix& operator=( const __scal_op_trans_mul_vec_plus_scal_vec<Matrix,Matrix> );
	
	// Sparse matrix operations
	
	Matrix& operator=( const __op_mul_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __op_trans_mul_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __scal_op_mul_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __scal_op_trans_mul_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __op_mul_vec_plus_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __op_trans_mul_vec_plus_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __scal_op_mul_vec_plus_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __scal_op_trans_mul_vec_plus_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __op_mul_vec_plus_scal_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __op_trans_mul_vec_plus_scal_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __scal_op_mul_vec_plus_scal_vec<SpMatrix,Matrix> );
	Matrix& operator=( const __scal_op_trans_mul_vec_plus_scal_vec<SpMatrix,Matrix> );
	
	void Print()
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

		MatFact( int _r, int _c ) : Matrix( _r, _c )
		{
			fact = false;
			mf = new double[r*c+1];
			ipiv = new integer[this->r+1];
			work = new double[4*(this->r)];
			iwork = new integer[this->r];
		}
		
		MatFact( const Matrix& M ) : Matrix( M )
		{
			fact = false;
			mf = new double[r*c+1];
			ipiv = new integer[this->r+1];
			work = new double[4*(this->r)];
			iwork = new integer[this->r];
		}
		
		~MatFact(){ delete[] iwork; delete[] work; delete[] ipiv; delete[] mf; }

		MatFact& operator=( MatFact& M )
		{
			Matrix::operator=( M );
			cout<<"Copying MatFact is not yet implemented, though the Matrix part will be copied normally\n";
			return *this;
		}
		
		void New() { fact = false; }
		
		int  Row() const { return Matrix::Row(); }
		int  Col() const { return Matrix::Col(); }
		
		void Fact();
		void Solve( Vector& X, const Vector& B, bool TRANS=false );
		void Solve( Matrix& X, const Matrix& B, bool TRANS=false );
};

typedef Array1D< Matrix > JagMatrix3D;
typedef Array1D< Vector > JagVector2D;

// Implementation of Matrix

inline void Matrix::AX( double* out, const double* in, double alpha, bool trans ) const
{
	mmx( trans ? Trans : NoTrans, out, in, alpha );
}

inline void Matrix::AX( Matrix& out, const Matrix& in, double alpha, bool trans ) const
{
	mmxm( trans ? Trans : NoTrans, out.m, out.r, in.m, in.r, alpha, in.c );
}

inline void Matrix::AXpY( double* out, const double* in, const double* Y, double alpha, double beta, bool trans ) const
{
	mmxpy( trans ? Trans : NoTrans, out, in, alpha, Y, beta );
}

inline void Matrix::AXpY( Matrix& out, const Matrix& in, const Matrix& Y, double alpha, double beta, bool trans ) const
{
	mmxmpym( trans ? Trans : NoTrans, out.m, out.r, in.m, in.r, alpha, Y.m, Y.r, beta, Y.c );
}

// with Vectors
inline void Matrix::AX( Vector& X, const Vector& B, double alpha, bool trans ) const
{ 
	if( !trans )
	{
		if( (X.Size() != Row())||(B.Size() != Col()) )
		{
			cout<<"Matrix::AX_V:bad dimensions\n";
			throw(12); return;
		}
	}else
	{
		if( (X.Size() != Col())||(B.Size() != Row()) ) 
		{
			cout<<"Matrix::AX_V:bad dimensions\n"; 
			throw(12); return;
		}
	}
	this->AX( X.v, B.v, alpha, trans ); 
}

inline void Matrix::AXpY( Vector& out, const Vector& X, const Vector& Y, double alpha, double beta, bool trans ) const
{

	if( !trans )
	{ 
		if( (out.Size() != Row())||(X.Size() != Col())||(Y.Size() != Row()) ) 
		{ 
			cout<<"BaseMatrix::AXpY_V:bad dimensions\n"; return; 
		} 
	}else
	{ 
		if( (out.Size() != Col())||(X.Size() != Row())||(Y.Size() != Col()) ) 
		{
			cout<<"BaseMatrix::AXpY_V:bad dimensions\n"; return; 
		}
	}
	this->AXpY( out.v, X.v, Y.v, alpha, beta, trans );
}

// End of implementation of Matrix

// Implementation of Vector

inline Vector& Vector::operator=( const __scal_vec<Vector> ) { std::cout<<"__scal_vec\n"; return *this; }

inline Vector& Vector::operator=( const __op_mul_vec<Matrix,Vector> op )
{
	op.op.mmx( NoTrans, this->v, op.vecA.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __op_trans_mul_vec<Matrix,Vector> op )
{
	op.op.mmx( Trans, this->v, op.vecA.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_mul_vec<Matrix,Vector> op )
{
	op.op.mmx( NoTrans, this->v, op.vecA.v, op.alpha );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_trans_mul_vec<Matrix,Vector> op )
{
	op.op.mmx( Trans, this->v, op.vecA.v, op.alpha );
	return *this;
}
inline Vector& Vector::operator=( const __op_mul_vec_plus_vec<Matrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, 1.0, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __op_trans_mul_vec_plus_vec<Matrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, 1.0, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_mul_vec_plus_vec<Matrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, op.alpha, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_trans_mul_vec_plus_vec<Matrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, op.alpha, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __op_mul_vec_plus_scal_vec<Matrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, 1.0, op.vecB.v, op.beta );
	return *this;
}
inline Vector& Vector::operator=( const __op_trans_mul_vec_plus_scal_vec<Matrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, 1.0, op.vecB.v, 1.0 );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_mul_vec_plus_scal_vec<Matrix,Vector> op )
{
	op.op.mmxpy( NoTrans, this->v, op.vecA.v, op.alpha, op.vecB.v, op.beta );
	return *this;
}
inline Vector& Vector::operator=( const __scal_op_trans_mul_vec_plus_scal_vec<Matrix,Vector> op )
{
	op.op.mmxpy( Trans, this->v, op.vecA.v, op.alpha, op.vecB.v, op.beta );
	return *this;
}

// End of implementation of Vector

// Implementation of Matrix

inline Matrix& Matrix::operator=( const __scal_vec<Matrix> ) { std::cout<<"__scal_vec\n"; return *this; }

inline Matrix& Matrix::operator=( const __op_mul_vec<Matrix,Matrix> op )
{
	op.op.mmxm( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_trans_mul_vec<Matrix,Matrix> op )
{
	op.op.mmxm( Trans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_mul_vec<Matrix,Matrix> op )
{
	op.op.mmxm( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_trans_mul_vec<Matrix,Matrix> op )
{
	op.op.mmxm( Trans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_mul_vec_plus_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_trans_mul_vec_plus_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_mul_vec_plus_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_trans_mul_vec_plus_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_mul_vec_plus_scal_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, op.beta, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __op_trans_mul_vec_plus_scal_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, 1.0, op.vecB.m, op.vecB.r, 1.0, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_mul_vec_plus_scal_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( NoTrans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, op.beta, op.vecA.c );
	return *this;
}
inline Matrix& Matrix::operator=( const __scal_op_trans_mul_vec_plus_scal_vec<Matrix,Matrix> op )
{
	op.op.mmxmpym( Trans, this->m, this->r, op.vecA.m, op.vecA.r, op.alpha, op.vecB.m, op.vecB.r, op.beta, op.vecA.c );
	return *this;
}

// End of implementation of Matrix


#endif // PDDESYS_H

#endif
