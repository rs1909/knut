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

	friend class BaseMatrix;
	friend class SpMatrix;
	friend class SpFact;
	friend class Matrix;
	friend class MatFact;

};

class Matrix;

// abstract matrix class
class BaseMatrix 
{

	public:
		virtual ~BaseMatrix() {}

		virtual int  Row() const = 0;
		virtual int  Col() const = 0;
		virtual void AX( double* X, const double* B, double alpha, bool trans ) const = 0;
		virtual void AX( Matrix& X, const Matrix& B, double alpha, bool trans ) const = 0;
		virtual void AXpY( double* out, const double* X, const double* Y, double alpha, double beta, bool trans ) const = 0;
		virtual void AXpY( Matrix& out, const Matrix& X, const Matrix& Y, double alpha, double beta, bool trans ) const = 0;
  
		// with Vectors
		void AX( Vector& X, const Vector& B, double alpha, bool trans ) const
		{ 
			if( !trans )
			{ 
				if( (X.Size() != Row())||(B.Size() != Col()) )
				{
					cout<<"BaseMatrix::AX_V:bad dimensions\n";
					throw(12); return;
				}
			}else
			{ 
				if( (X.Size() != Col())||(B.Size() != Row()) ) 
				{ 
					cout<<"BaseMatrix::AX_V:bad dimensions\n"; 
					throw(12); return;
				} 
			}
			this->AX( X.v, B.v, alpha, trans ); 
		}
		
		void AXpY( Vector& out, const Vector& X, const Vector& Y, double alpha, double beta, bool trans ) const
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
		
		// other variants
		void AX( Vector& X, const Vector& B ) const { this->AX( X, B, 1.0, false ); }
		void AX( Matrix& X, const Matrix& B ) const { this->AX( X, B, 1.0, false ); }
		void AtX( Vector& X, const Vector& B ) const { this->AX( X, B, 1.0, true ); }
		void AtX( Matrix& X, const Matrix& B ) const { this->AX( X, B, 1.0, true ); }
		void mAX( Vector& X, const Vector& B ) const { this->AX( X, B, -1.0, false ); }
		void mAX( Matrix& X, const Matrix& B ) const { this->AX( X, B, -1.0, false ); }
		void mAtX( Vector& X, const Vector& B ) const { this->AX( X, B, -1.0, true ); }
		void mAtX( Matrix& X, const Matrix& B ) const { this->AX( X, B, -1.0, true ); }
};

#ifndef PDDESYS_H

class BaseFact : virtual public BaseMatrix
{
	public:
		virtual ~BaseFact() {}
		
		virtual void New() = 0;
		virtual void Solve( Vector& X, const Vector& B, bool TRANS=false ) = 0;
		virtual void Solve( Matrix& X, const Matrix& B, bool TRANS=false ) = 0;
};
#endif // PDDESYS_H

class Matrix : public Array2D<double>
#ifndef PDDESYS_H
	, virtual public BaseMatrix
#endif // PDDESYS_H
{
 public:

	Matrix() { }

	Matrix( int i, int j ) : Array2D<double>( i, j ) { }

	Matrix( const Matrix& M ) : 
#ifndef PDDESYS_H
	BaseMatrix(),
#endif // PDDESYS_H
	Array2D<double>( M ) { }

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
	
  friend class MatFact;
  friend class SpMatrix;
  friend class SpFact;
};

#ifndef PDDESYS_H

class MatFact : public Matrix, public BaseFact
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
	
		void AX( double* X, const double* B, double alpha, bool trans ) const 
		{
			Matrix::AX( X, B, alpha, trans ); 
		}
	
		void AX( Vector& X, const Vector& B, double alpha, bool trans ) const  
		{ 
			BaseMatrix::AX( X, B, alpha, trans ); 
		}
	
		void AX( Matrix& X, const Matrix& B, double alpha, bool trans ) const 
		{ 
			Matrix::AX( X, B, alpha, trans ); 
		}
	
		void AXpY( double* out, const double* X, const double* Y, double alpha, double beta, bool trans ) const 
		{ 
			Matrix::AXpY( out, X, Y, alpha, beta, trans ); 
		}
	
		void AXpY( Vector& out, const Vector& X, const Vector& Y, double alpha, double beta, bool trans ) const 
		{ 
			BaseMatrix::AXpY( out, X, Y, alpha, beta, trans ); 
		}
	
		void AXpY( Matrix& out, const Matrix& X, const Matrix& Y, double alpha, double beta, bool trans ) const 
		{ 
			Matrix::AXpY( out, X, Y, alpha, beta, trans ); 
		}
  
		//   void Clear() { Matrix::Clear(); fact = false; }
		void Fact();
		void Solve( Vector& X, const Vector& B, bool TRANS=false );
		void Solve( Matrix& X, const Matrix& B, bool TRANS=false );
};

typedef Array1D< Matrix > JagMatrix3D;
typedef Array1D< Vector > JagVector2D;

#endif // PDDESYS_H

#endif
