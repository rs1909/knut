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

#include <iostream>
#include <iterator>

#ifndef PDDESYS_H
#include "plot.h"
#include "pderror.h"

extern "C"
{

#include "cspblas.h"
#include "cblas.h"
#include "laarpack.h"

}
#endif

#ifndef P_ASSERT
#  define P_ASSERT(cond) do{ if(!(cond)) { std::cout<<#cond; std::cout.flush(); abort(); } }while(0)
#endif

#ifndef P_ASSERT_X
#  define P_ASSERT_X(cond,msg) do{ if(!(cond)) { std::cout<<#cond<<msg; std::cout.flush(); abort(); } }while(0)
#endif

template<class T> class Array2D;
template<class T> class Array3D;

template<class T>
class Array1D
{

  protected:

    int n;
    T  *v;
  private:
    const bool destructable;
  public:

    inline Array1D() : n(0), v(0), destructable(true)
    { }

    inline Array1D(bool b) : n(0), v(0), destructable(false)
    { }

    inline Array1D(int i) : n(i), v(new T[i]), destructable(true)
    {
      Clear();
    }

    // specially for JagMatrix2D i.e. Array1D< Vector >, which is indexed as (i)(j)
    inline Array1D(int i, int j) : n(i), v(new T[i]), destructable(true)
    {
      for (int r = 0; r < i; r++) v[r].Init(j);
    }

    // specially for JagMatrix3D i.e. Array1D< Matrix >, which is indexed as (i)(j,k)
    inline Array1D(int i, int j, int k) : n(i), v(new T[i]), destructable(true)
    {
      for (int r = 0; r < i; r++) v[r].Init(j, k);
    }

    inline Array1D(const Array1D<T>& V_) : n(V_.n), v(new T[V_.n]), destructable(true)
    {
      *this = V_;
    }

    inline Array1D(T *data, int i) : n(i), v(data), destructable(false)
    { }

    inline Array1D(const Array2D<T>& v, int i);

    inline virtual ~Array1D()
    {
      if (destructable) delete[] v;
    }

    inline void Init(int i)
    {
      if (destructable)
      {
        delete[] v;
        n = i;
        v = new T[i];
        Clear();
      } else throw(-1);
    }

    inline void Init(const Array1D<T>& V_)
    {
      if (destructable)
      {
        delete[] v;
        n = V_.n;
        v = new T[V_.n];
        *this = V_;
      } else throw(-1);
    }

    inline void Init(T *data, int i)
    {
      if (!destructable)
      {
        n = i;
        v = data;
      } else throw(-1);
    }

    inline void Clear()
    {
      for (int j = 0; j < n; j++) v[j] = T();
    } // it was T(0), but for built in types it must be the same

    inline T *Pointer()
    {
      return v;
    }

    inline T *Pointer(const int i)
    {
#ifdef DEBUG
      P_ASSERT_X(i < n && i >= 0, "Array1D::bound11&");
#endif
      return &v[i];
    }

    inline int Size() const
    {
      return n;
    }

    inline Array1D<T>& operator=(const Array1D<T>& V)
    {
#ifdef DEBUG
      P_ASSERT_X(n == V.n, "Array1D::operator=(): incompatible sizes\n");
#endif
      for (int i = 0; i < V.n; i++) v[i] = V.v[i];
      return *this;
    }

    inline T& operator()(int i)
    {
#ifdef DEBUG
      P_ASSERT_X(i < n && i >= 0, "Array1D::bound11&");
#endif
      return v[i];
    }

    inline T operator()(int i) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < n && i >= 0 , "vec_bound11_\n");
#endif
      return v[i];
    }
    class iterator
    {
      private:
        T *pt;
      public:
        iterator() : pt(0) {}
        iterator( const Array1D<T>::iterator& it ) : pt(it.pt) {}
        iterator& operator= (const iterator& it) { pt = it.pt; return *this; }
        iterator  operator+ ( int i ) const { iterator it; it.pt = pt+i; return it; }
        iterator  operator- ( int i ) const { iterator it; it.pt = pt-i; return it; }
        int       operator- (const iterator& it) const { return pt - it.pt; }
        bool      operator!= (const iterator& it) const { return pt != it.pt; }
        bool      operator== (const iterator& it) const { return pt == it.pt; }
        const T&  operator*() const { return *pt; }
        T&  operator*() { return *pt; }
        iterator& operator++() { ++pt; return *this; }
        iterator& operator--() { --pt; return *this; }
        iterator& operator++(int) { ++pt; return *this; }
        iterator& operator--(int) { --pt; return *this; }
        bool      operator< (const iterator& it) const { return pt < it.pt; }
        friend class Array1D<T>;
        typedef std::bidirectional_iterator_tag iterator_category;
        typedef T                        value_type;
        typedef std::ptrdiff_t           difference_type;
        typedef T*                       pointer;
        typedef T&                       reference;
    };
    iterator begin() { iterator it; it.pt = v; return it; }
    iterator end() { iterator it; it.pt = v+n; return it; }
};


template<class T>
class Array2D
{

  protected:

    T* m;
    int r, c;
  private:
    const bool destructable;
  public:
    inline Array2D() : destructable(true)
    {
      r = 0;
      c = 0;
      m = 0;
    }

    inline Array2D(int _r, int _c) : r(_r), c(_c), destructable(true)
    {
      m = new T[r*c+1];
      Clear();
    }

    inline Array2D(const Array2D<T>& M) : r(M.r), c(M.c), destructable(true)
    {
      m = new T[r*c+1];
      for (int i = 0; i < r*c; i++) m[i] = M.m[i];
    }

    inline Array2D(const Array3D<T>& v, int i);

    inline virtual ~Array2D()
    {
      if (destructable) delete[] m;
    }
    inline int dim1() const { return r; }
    inline int dim2() const { return c; }

    inline void Init(int _r, int _c)
    {
      P_ASSERT_X(destructable, "Array2D<T>::Init : trying to resize a non-destructable a");
      delete[] m;
      r = _r;
      c = _c;
      m = new T[r*c+1];
      Clear();
    }

    inline T *Pointer()
    {
      return m;
    }

    inline T *Pointer(const int i, const int j)
    {
#ifdef DEBUG
      P_ASSERT_X(i < r && j < c, "bound& ");
      P_ASSERT_X(i >= 0 && j >= 0, "lbound& ");
#endif
      return &m[i + r*j];
    }

    inline void Clear()
    {
      for (int i = 0; i < r*c; i++) m[i] = T(0);
    }

    inline Array2D<T>& operator= (const Array2D<T>& M)
    {
#ifdef DEBUG
      P_ASSERT_X(M.r == r && M.c == c, "Array2D<T>::operator= : incompatible sizes");
#endif
      for (int i = 0; i < r*c; i++) m[i] = M.m[i];
      return *this;
    }

    inline T& operator()(const int i, const int j)
    {
#ifdef DEBUG
      P_ASSERT_X(i < r && j < c, "bound& ");
      P_ASSERT_X(i >= 0 && j >= 0, "lbound& ");
#endif
      return m[i + r*j];
    }

    inline T operator()(const int i, const int j) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < r && j < c, "bound_ ");
      P_ASSERT_X(i >= 0 && j >= 0, "lbound_ ");
#endif
      return m[i + r*j];
    }

    inline T& operator()(const int i)
    {
#ifdef DEBUG
      P_ASSERT_X(i < r*c, "bound11&\n");
      P_ASSERT_X(r == 1 || c == 1, "H&\n");
#endif
      return m[i];
    }

    inline T operator()(const int i) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < r*c, "bound11_\n");
      P_ASSERT_X(r == 1 || c == 1, "H&\n");
#endif
      return m[i];
    }
    friend class Array1D<T>;
};

template<class T>
class Array3D
{

  protected:

    T* m;
    int d1, d2, d3;

  public:
    inline Array3D()
    {
      d1 = d2 = d3 = 0;
      m = 0;
    }

    inline Array3D(int _d1, int _d2, int _d3) : d1(_d1), d2(_d2), d3(_d3)
    {
      m = new T[d1*d2*d3+1];
      Clear();
    }

    inline Array3D(const Array3D<T>& M) : d1(M.d1), d2(M.d2), d3(M.d3)
    {
      m = new T[d1*d2*d3+1];
      for (int i = 0; i < d1*d2*d3; i++) m[i] = M.m[i];
    }

    inline virtual ~Array3D()
    {
      delete[] m;
    }

    inline int dim1() const { return d1; }
    inline int dim2() const { return d2; }
    inline int dim3() const { return d3; }
    inline void Init(int _d1, int _d2, int _d3)
    {
      delete[] m;
      d1 = _d1;
      d2 = _d2;
      d3 = _d3;
      m = new T[d1*d2*d3+1];
      Clear();
    }

    inline void Clear()
    {
      for (int i = 0; i < d1*d2*d3; i++) m[i] = 0;
    }

    inline Array3D<T>& operator= (const Array3D<T>& M)
    {
#ifdef DEBUG
      P_ASSERT_X((M.d1 == d1) && (M.d2 == d2) && (M.d3 == d3), "Array3D<T>::operator= : incompatible sizes\n");
#endif
      for (int i = 0; i < d1*d2*d3; i++) m[i] = M.m[i];
      return *this;
    }

    inline T& operator()(const int i, const int j, const int k)
    {
#ifdef DEBUG
      P_ASSERT_X((i < d1) && (j < d2) && (k < d3), "bound&\n");
      P_ASSERT_X((i >= 0) && (j >= 0) && (k >= 0), "lbound&\n");
#endif
      return m[i + d1*(j + d2*k)];
    }

    inline T operator()(const int i, const int j, const int k) const
    {
#ifdef DEBUG
      P_ASSERT_X((i < d1) && (j < d2) && (k < d3), "bound\n");
      P_ASSERT_X((i >= 0) && (j >= 0) && (k >= 0), "lbound\n");
#endif
      return m[i + d1*(j + d2*k)];
    }
    friend class Array2D<T>;
};

template<class T> inline Array1D<T>::Array1D(const Array2D<T>& vec, int i) : destructable(false)
{
#ifdef DEBUG
  P_ASSERT_X((i < vec.c), "bound\n");
#endif
  v = &vec.m[vec.r*i];
  n = vec.r;
}

template<class T> inline Array2D<T>::Array2D(const Array3D<T>& vec, int i) : destructable(false)
{
#ifdef DEBUG
  P_ASSERT_X((i < vec.d3), "bound\n");
#endif
  m = &vec.m[vec.d1*vec.d2*i];
  r = vec.d1;
  c = vec.d2;
}

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
    inline rng() : i1(0), i2(0), j1(0), j2(0)
    { }
    inline rng(int i1_, int i2_) : i1(i1_), i2(i2_), j1(0), j2(0)
    { }
    inline rng(int i1_, int i2_, int j1_, int j2_) : i1(i1_), i2(i2_), j1(j1_), j2(j2_)
    { }

    const int i1, i2, j1, j2;
};

/// with range
template< class VT > class __scal_vec_trans_rng : public rng
{
  public:
    inline __scal_vec_trans_rng(const VT& v)                                        : rng(0, v.Row(), 0, v.Col()), vec(v), alpha(1.0), tr(NoTrans)
    { }
    inline __scal_vec_trans_rng(const VT& v, double a)                              : rng(0, v.Row(), 0, v.Col()), vec(v), alpha(a),   tr(NoTrans)
    { }
    inline __scal_vec_trans_rng(const VT& v, double a, enum cspblas_Trans t)        : rng(0, v.Row(), 0, v.Col()), vec(v), alpha(a),   tr(t)
    { }
    inline __scal_vec_trans_rng(const VT& v, rng r)                                 : rng(r), vec(v), alpha(1.0), tr(NoTrans)
    { }
    inline __scal_vec_trans_rng(const VT& v, rng r, double a)                       : rng(r), vec(v), alpha(a),   tr(NoTrans)
    { }
    inline __scal_vec_trans_rng(const VT& v, rng r, double a, enum cspblas_Trans t) : rng(r), vec(v), alpha(a),   tr(t)
    { }
    inline __scal_vec_trans_rng(__scal_vec_trans<VT> v);

    inline __scal_vec_trans_rng< VT >     operator!()
    {
      return __scal_vec_trans_rng<VT>(alpha, *this, Trans);
    }
    inline __scal_vec_trans_rng< VT >     operator-()
    {
      return __scal_vec_trans_rng<VT>(-alpha, *this, tr);
    }
    inline __op_mul_vec_rng< VT, Vector > operator*(__scal_vec_trans_rng<Vector> v);
    inline __op_mul_vec_rng< VT, Matrix > operator*(__scal_vec_trans_rng<Matrix> v);
    inline __op_mul_vec_rng< VT, Vector > operator*(const Vector& v);
    inline __op_mul_vec_rng< VT, Matrix > operator*(const Matrix& v);

    const VT& vec;
    const double alpha;
    const enum cspblas_Trans tr;
};

/// with range
template< class MT, class VT > class __op_mul_vec_rng
{
  public:
    inline __op_mul_vec_rng(__scal_vec_trans_rng<MT> m, __scal_vec_trans_rng<VT> v) : op(m), vecA(v)
    { }

    inline __op_mul_vec_plus_vec_rng<MT, VT>   operator+(__scal_vec_trans_rng<VT> v);
    inline __op_mul_vec_plus_vec_rng<MT, VT>   operator+(const VT& v);
    inline __op_mul_vec_plus_vec_rng<MT, VT>   operator+(__scal_vec_trans<VT> v);
    inline __op_mul_vec_plus_vec_rng<MT, VT>   operator-(__scal_vec_trans_rng<VT> v);
    inline __op_mul_vec_plus_vec_rng<MT, VT>   operator-(const VT& v);
    inline __op_mul_vec_plus_vec_rng<MT, VT>   operator-(__scal_vec_trans<VT> v);

    const __scal_vec_trans_rng<MT> op;
    const __scal_vec_trans_rng<VT> vecA;
};

/// with range: this is the output
template< class MT, class VT > class __op_mul_vec_plus_vec_rng : public __op_mul_vec_rng<MT, VT>
{
  public:
    inline __op_mul_vec_plus_vec_rng(__op_mul_vec_rng<MT, VT> m, __scal_vec_trans_rng<VT> v) : __op_mul_vec_rng<MT, VT>(m), vecB(v)
    { }

    const __scal_vec_trans_rng<VT> vecB;
};

/// WITHOUT RANGES
template< class VT > class __scal_vec_trans
{
  public:
    inline __scal_vec_trans(const VT& m) : vec(m), alpha(1.0), tr(NoTrans)
    { }
    inline __scal_vec_trans(const VT& m, double a) : vec(m), alpha(a), tr(NoTrans)
    { }
    inline __scal_vec_trans(const VT& m, double a, enum cspblas_Trans t) : vec(m), alpha(a), tr(t)
    { }

    inline __scal_vec_trans< VT >                        operator!()
    {
      return __scal_vec_trans<VT>(alpha, vec, Trans);
    }
    inline __scal_vec_trans< VT >                        operator-()
    {
      return __scal_vec_trans<VT>(-alpha, vec, tr);
    }
    inline __op_mul_vec< VT, Vector >                    operator*(const Vector& v);
    inline __op_mul_vec< VT, Matrix >                    operator*(const Matrix& v);
    inline __op_mul_vec_rng< VT, Vector >                operator*(__scal_vec_trans_rng<Vector> v);
    inline __op_mul_vec_rng< VT, Matrix >                operator*(__scal_vec_trans_rng<Matrix> v);

    const VT& vec;
    double alpha;
    const enum cspblas_Trans tr;
};

template< class MT, class VT > class __op_mul_vec
{
  public:
    inline __op_mul_vec(__scal_vec_trans<MT> m, __scal_vec_trans<VT> v) : op(m), vecA(v)
    { }

    inline __op_mul_vec_plus_vec<MT, VT>      operator+(__scal_vec_trans<VT> v);
    inline __op_mul_vec_plus_vec<MT, VT>      operator-(__scal_vec_trans<VT> v);
    inline __op_mul_vec_plus_vec_rng<MT, VT>  operator+(__scal_vec_trans_rng<VT> v);
    inline __op_mul_vec_plus_vec_rng<MT, VT>  operator-(__scal_vec_trans_rng<VT> v);
    inline __op_mul_vec_plus_vec<MT, VT>      operator+(const VT& v);
    inline __op_mul_vec_plus_vec<MT, VT>      operator-(const VT& v);

    const __scal_vec_trans<MT> op;
    const __scal_vec_trans<VT> vecA;
};

template< class MT, class VT > class __op_mul_vec_plus_vec : public __op_mul_vec<MT, VT>
{
  public:
    inline __op_mul_vec_plus_vec(__op_mul_vec<MT, VT> m, __scal_vec_trans<VT> v) : __op_mul_vec<MT, VT>(m), vecB(v)
    { }

    const __scal_vec_trans<VT> vecB;
};

#endif // PDDESYS_H

class Vector : public Array1D<double>
{

  public:

    inline Vector()
    { }

    inline Vector(int i) : Array1D<double>(i)
    { }

    inline Vector(bool b) : Array1D<double>(b)
    { }

    inline Vector(const Vector& V) : Array1D<double>(V)
    { }

    inline Vector(const Array2D<double>& m, int i) : Array1D<double>(m,i)
    { }

    inline virtual ~Vector()
    { }

    inline void Print() const
    {
      for (int j = 0; j < n; j++) std::cout << v[j] << '\t';
      std::cout << '\n';
    }

#ifndef PDDESYS_H

    inline void Rand()
    {
      static int idist = 2;
      static int iseed[4] =
        {
          1, 3, 5, 7
        };
      int N = static_cast<int>(n);
      pdde_dlarnv(&idist, iseed, &N, v);
    }

    inline Vector& operator= (const Vector& V);
    inline Vector& operator+=(const Vector& V);
    inline Vector& operator-=(const Vector& V);
    inline Vector& operator/=(double div);
    inline Vector& operator*=(double mul);
    inline double  operator*(const Vector& V) const;

    // Matrix operations
    inline Vector& operator= (const __scal_vec_trans<Vector>);
    inline Vector& operator+=(const __scal_vec_trans<Vector>);
    inline Vector& operator-=(const __scal_vec_trans<Vector>);
    inline Vector& operator= (const __op_mul_vec<Matrix, Vector>);
    inline Vector& operator+=(const __op_mul_vec<Matrix, Vector>);
    inline Vector& operator-=(const __op_mul_vec<Matrix, Vector>);
    inline Vector& operator+=(const __op_mul_vec_plus_vec<Matrix, Vector>);
    inline Vector& operator-=(const __op_mul_vec_plus_vec<Matrix, Vector>);
    // for SpMatrix
    inline Vector& operator= (const __op_mul_vec<SpMatrix, Vector>);
    // obsolote
    inline Vector& operator= (const __op_mul_vec_plus_vec<SpMatrix, Vector>);
    /// with ranges
    inline Vector& operator=(const __scal_vec_trans_rng<Vector>);
    inline Vector& operator=(const __op_mul_vec_rng<Matrix, Vector>);
    // obsolote
    inline Vector& operator=(const __op_mul_vec_plus_vec_rng<Matrix, Vector>);
    // for SpMatrix
    inline Vector& operator=(const __op_mul_vec_rng<SpMatrix, Vector>);
    // obsolote
    inline Vector& operator=(const __op_mul_vec_plus_vec_rng<SpMatrix, Vector>);

    inline __scal_vec_trans<Vector>     operator-() const
    {
      return __scal_vec_trans<Vector>(*this, -1.0);
    }
    inline __scal_vec_trans_rng<Vector> operator[ ](const rng r) const
    {
      return __scal_vec_trans_rng<Vector>(*this, r);
    }

#endif // PDDESYS_H

    friend class SpMatrix;
    friend class SpFact;
    friend class Matrix;
    friend class MatFact;

};

class Matrix : public Array2D<double>
{
  public:

    inline Matrix()
    { }

    inline Matrix(int i, int j) : Array2D<double>(i, j)
    { }

    inline Matrix(const Matrix& M) : Array2D<double>(M)
    { }

    inline Matrix(const Array3D<double>& m, int i) : Array2D<double>(m,i)
    { }

    inline virtual ~Matrix()
    { }

    inline int Row() const
    {
      return r;
    }
    inline int Col() const
    {
      return c;
    }
    inline int Size() const
    {
      if ((r == 1) || (c == 1))
      {
        return r*c;
      }
      else
      {
        std::cout << "Hs\n";
        return 0;
      }
    }

#ifndef PDDESYS_H

    inline Matrix& operator= (const __scal_vec_trans<Matrix>);
    inline Matrix& operator= (const __op_mul_vec<Matrix, Matrix>);
    inline Matrix& operator+=(const __op_mul_vec<Matrix, Matrix>);
    inline Matrix& operator-=(const __op_mul_vec<Matrix, Matrix>);
    // obsolote
    inline Matrix& operator=(const __op_mul_vec_plus_vec<Matrix, Matrix>);
    // with SpMatrix
    inline Matrix& operator=(const __op_mul_vec<SpMatrix, Matrix>);
    // obsolote
    inline Matrix& operator=(const __op_mul_vec_plus_vec<SpMatrix, Matrix>);
    /// with ranges
    inline Matrix& operator=(const __scal_vec_trans_rng<Matrix>);
    inline Matrix& operator=(const __op_mul_vec_rng<Matrix, Matrix>);
    // obsolote
    inline Matrix& operator=(const __op_mul_vec_plus_vec_rng<Matrix, Matrix>);
    // with SpMatrix
    inline Matrix& operator=(const __op_mul_vec_rng<SpMatrix, Matrix>);
    // obsolote
    inline Matrix& operator=(const __op_mul_vec_plus_vec_rng<SpMatrix, Matrix>);

    void Eigval(Vector& re, Vector& im);
    void Eigval(Vector& re, Vector& im, Matrix& lev, Matrix& rev);
    void StrPlot(GnuPlot& pl);

    /* operators */

    inline __op_mul_vec<Matrix, Vector>     operator*(const Vector& v) const
    {
      return __op_mul_vec<Matrix, Vector>(__scal_vec_trans<Matrix>(*this), __scal_vec_trans<Vector>(v));
    }
    inline __op_mul_vec<Matrix, Matrix>     operator*(const Matrix& v) const
    {
      return __op_mul_vec<Matrix, Matrix>(__scal_vec_trans<Matrix>(*this), __scal_vec_trans<Matrix>(v));
    }
    inline __op_mul_vec_rng<Matrix, Vector> operator*(const __scal_vec_trans_rng<Vector> v) const
    {
      return __op_mul_vec_rng<Matrix, Vector>(__scal_vec_trans_rng<Matrix>(*this), v);
    }
    inline __op_mul_vec_rng<Matrix, Matrix> operator*(const __scal_vec_trans_rng<Matrix> v) const
    {
      return __op_mul_vec_rng<Matrix, Matrix>(__scal_vec_trans_rng<Matrix>(*this), v);
    }
    inline __scal_vec_trans<Matrix>        operator!() const
    {
      return __scal_vec_trans<Matrix>(*this, 1.0, Trans);
    }
    inline __scal_vec_trans_rng<Matrix>    operator[ ](const rng r_) const
    {
      return __scal_vec_trans_rng<Matrix>(*this, r_);
    }

#endif // PDDESYS_H

    inline void Print()
    {
      double sum = 0.0;
      for (int i = 0; i < r; i++)
      {
        for (int j = 0; j < c; j++)
        {
          std::cout << (*this)(i, j) << '\t';
          sum += (*this)(i, j);
        }
        std::cout << '\n';
      }
      std::cout << "SUM: " << sum << '\n';
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
    int* ipiv;
    double*  work;
    int* iwork;
    int  info;
    // for solving
    double rcond;
    double ferr;  // in case of Vector these are just double
    double berr;

  public:

    inline MatFact(int _r, int _c) : Matrix(_r, _c)
    {
      fact = false;
      mf = new double[r*c+1];
      ipiv = new int[this->r+1];
      work = new double[4*(this->r)];
      iwork = new int[this->r];
    }

    inline MatFact(const Matrix& M) : Matrix(M)
    {
      fact = false;
      mf = new double[r*c+1];
      ipiv = new int[this->r+1];
      work = new double[4*(this->r)];
      iwork = new int[this->r];
    }

    inline virtual ~MatFact()
    {
      delete[] iwork;
      delete[] work;
      delete[] ipiv;
      delete[] mf;
    }

    inline MatFact& operator=(MatFact& M)
    {
      Matrix::operator=(M);
      std::cout << "Copying MatFact is not yet implemented, though the Matrix part will be copied normally\n";
      return *this;
    }

    inline void New()
    {
      fact = false;
    }

    inline int  Row() const
    {
      return Matrix::Row();
    }
    inline int  Col() const
    {
      return Matrix::Col();
    }

    void Fact();
    void Solve(Vector& X, const Vector& B, bool TRANS = false);
    void Solve(Matrix& X, const Matrix& B, bool TRANS = false);
};

typedef Array1D< Matrix > JagMatrix3D;
typedef Array1D< Vector > JagVector2D;

/// specialized versions of the Clear function
template< > inline void Array1D< Array1D<int> >::Clear()
{
  for (int i = 0; i < n; i++) v[i].Clear();
}
template< > inline void Array1D< Array1D<double> >::Clear()
{
  for (int i = 0; i < n; i++) v[i].Clear();
}
template< > inline void Array1D< Vector >::Clear()
{
  for (int i = 0; i < n; i++) v[i].Clear();
}
template< > inline void Array1D< Array2D<int> >::Clear()
{
  for (int i = 0; i < n; i++) v[i].Clear();
}
template< > inline void Array1D< Array2D<double> >::Clear()
{
  for (int i = 0; i < n; i++) v[i].Clear();
}
template< > inline void Array1D< Matrix >::Clear()
{
  for (int i = 0; i < n; i++) v[i].Clear();
}

template< > inline Array1D< Vector >::Array1D(const Array1D< Vector >& V_) : destructable(true)
{
  n = V_.n;
  v = new Vector[V_.n];
  for (int i = 0; i < V_.n; i++) v[i].Init(V_.v[i]);
}

/// Member functions and Operators for __scal_vec_trans_rng
template<> inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng(const Vector& v)
    : rng(0, v.Size()), vec(v), alpha(1.0), tr(NoTrans)
{ }

template<> inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng(const Vector& v, double a)
    : rng(0, v.Size()), vec(v), alpha(a),   tr(NoTrans)
{ }

template<> inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng(const Vector& v, double a, enum cspblas_Trans t)
    : rng(0, v.Size()), vec(v), alpha(a),   tr(t)
{ }

template< class VT > inline __scal_vec_trans_rng<VT>::__scal_vec_trans_rng(__scal_vec_trans<VT> v)
    : rng(0, v.vec.Row(), 0, v.vec.Col()), vec(v.vec), alpha(v.alpha), tr(v.tr)
{ }

template< >          inline __scal_vec_trans_rng<Vector>::__scal_vec_trans_rng(__scal_vec_trans<Vector> v)
    : rng(0, v.vec.Size()), vec(v.vec), alpha(v.alpha), tr(v.tr)
{ }

template< class VT > inline __op_mul_vec_rng< VT, Vector > __scal_vec_trans_rng<VT>::operator*(__scal_vec_trans_rng<Vector> v)
{
  return __op_mul_vec_rng< VT, Vector >(*this, v);
}

template< class VT > inline __op_mul_vec_rng< VT, Matrix > __scal_vec_trans_rng<VT>::operator*(__scal_vec_trans_rng<Matrix> v)
{
  return __op_mul_vec_rng< VT, Matrix >(*this, v);
}

template< class VT > inline __op_mul_vec_rng< VT, Vector > __scal_vec_trans_rng<VT>::operator*(const Vector& v)
{
  return __op_mul_vec_rng< VT, Vector >(*this, __scal_vec_trans_rng<Vector>(v));
}

template< class VT > inline __op_mul_vec_rng< VT, Matrix > __scal_vec_trans_rng<VT>::operator*(const Matrix& v)
{
  return __op_mul_vec_rng< VT, Matrix >(*this, __scal_vec_trans_rng<Matrix>(v));
}


/// Member functions and Operators for __op_mul_vec_rng
template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>   __op_mul_vec_rng<MT, VT>::operator+(__scal_vec_trans_rng<VT> v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(*this, v);
}

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>   __op_mul_vec_rng<MT, VT>::operator+(const VT& v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(*this, __scal_vec_trans_rng<VT>(v));
}

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>   __op_mul_vec_rng<MT, VT>::operator+(__scal_vec_trans<VT> v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(*this, __scal_vec_trans_rng<VT>(v));
}

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>   __op_mul_vec_rng<MT, VT>::operator-(__scal_vec_trans_rng<VT> v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(*this, __scal_vec_trans_rng<VT>(v, -v.alpha));
}

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>   __op_mul_vec_rng<MT, VT>::operator-(const VT& v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(*this, __scal_vec_trans_rng<VT>(v, -1.0));
}

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>   __op_mul_vec_rng<MT, VT>::operator-(__scal_vec_trans<VT> v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(*this, __scal_vec_trans_rng<VT>(v.vec, -v.alpha, v.tr));
}


/// Member functions and Operators for __scal_vec_trans
template< class VT > inline __op_mul_vec< VT, Vector >     __scal_vec_trans<VT>::operator*(const Vector& v)
{
  return __op_mul_vec< VT, Vector >(*this, __scal_vec_trans<Vector>(v));
}

template< class VT > inline __op_mul_vec< VT, Matrix >     __scal_vec_trans<VT>::operator*(const Matrix& v)
{
  return __op_mul_vec< VT, Matrix >(*this, __scal_vec_trans<Matrix>(v));
}

template< class VT > inline __op_mul_vec_rng< VT, Vector > __scal_vec_trans<VT>::operator*(__scal_vec_trans_rng<Vector> v)
{
  return __op_mul_vec_rng< VT, Vector >(__scal_vec_trans_rng<VT>(*this), v);
}

template< class VT > inline __op_mul_vec_rng< VT, Matrix > __scal_vec_trans<VT>::operator*(__scal_vec_trans_rng<Matrix> v)
{
  return __op_mul_vec_rng< VT, Matrix >(__scal_vec_trans_rng<VT>(*this), v);
}


/// Member functions and Operators for __op_mul_vec
template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>  __op_mul_vec<MT, VT>::operator+(__scal_vec_trans_rng<VT> v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(__op_mul_vec_rng<MT, VT>(this->op, this->vecA), v);
}

template< class MT, class VT > inline __op_mul_vec_plus_vec_rng<MT, VT>  __op_mul_vec<MT, VT>::operator-(__scal_vec_trans_rng<VT> v)
{
  return __op_mul_vec_plus_vec_rng<MT, VT>(__op_mul_vec_rng<MT, VT>(this->op, this->vecA), __scal_vec_trans_rng<VT>(v.vec, -v.alpha, v.tr));
}

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT, VT>      __op_mul_vec<MT, VT>::operator+(__scal_vec_trans<VT> v)
{
  return __op_mul_vec_plus_vec<MT, VT>(*this, v);
}

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT, VT>      __op_mul_vec<MT, VT>::operator-(__scal_vec_trans<VT> v)
{
  return __op_mul_vec_plus_vec<MT, VT>(*this, __scal_vec_trans<VT>(v.vec, -v.alpha, v.tr));
}

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT, VT>      __op_mul_vec<MT, VT>::operator+(const VT& v)
{
  return __op_mul_vec_plus_vec<MT, VT>(*this, __scal_vec_trans<VT>(v));
}

template< class MT, class VT > inline __op_mul_vec_plus_vec<MT, VT>      __op_mul_vec<MT, VT>::operator-(const VT& v)
{
  return __op_mul_vec_plus_vec<MT, VT>(*this, __scal_vec_trans<VT>(v, -1.0));
}


/// The global operators
inline __scal_vec_trans < Vector >                      operator*(double a, const Vector& v)
{
  return __scal_vec_trans<Vector>(v, a);
}

inline __scal_vec_trans < Matrix >                      operator*(double a, const Matrix& v)
{
  return __scal_vec_trans<Matrix>(v, a);
}

inline __scal_vec_trans < SpMatrix >                    operator*(double a, const SpMatrix& v)
{
  return __scal_vec_trans<SpMatrix>(v, a);
}

template< class VT > inline __scal_vec_trans < VT  >    operator*(double a, __scal_vec_trans< VT > v)
{
  return __scal_vec_trans<VT>(v.vec, a * v.alpha, v.tr);
}

template< class VT > inline __scal_vec_trans_rng < VT > operator*(double a, __scal_vec_trans_rng< VT > v)
{
  return __scal_vec_trans_rng<VT>(v.vec, v, a * v.alpha, v.tr);
}

inline __scal_vec_trans < Vector >                      operator-(const Vector& v)
{
  return __scal_vec_trans<Vector>(v, -1.0);
}

inline __scal_vec_trans < Matrix >                      operator-(const Matrix& v)
{
  return __scal_vec_trans<Matrix>(v, -1.0);
}

inline __scal_vec_trans < SpMatrix >                    operator-(const SpMatrix& v)
{
  return __scal_vec_trans<SpMatrix>(v, -1.0);
}

// Implementation of Vector

inline Vector& Vector::operator=(const Vector& V)
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "Vector::operator=(): incompatible sizes\n");
#endif // DEBUG
  cblas_dcopy(n, V.v, 1, v, 1);
  return *this;
}

inline Vector& Vector::operator+=(const Vector& V)
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "Vector::operator+=(): incompatible sizes\n");
#endif // DEBUG
  cblas_daxpy(n, 1.0, V.v, 1, v, 1);
  return *this;
}

inline Vector& Vector::operator-=(const Vector& V)
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "Vector::operator-=(): incompatible sizes\n");
#endif // DEBUG
  cblas_daxpy(n, -1.0, V.v, 1, v, 1);
  return *this;
}

inline Vector& Vector::operator/=(double div)
{
  cblas_dscal(n, 1.0 / div, v, 1);
  return *this;
}

inline Vector& Vector::operator*=(double mul)
{
  cblas_dscal(n, mul, v, 1);
  return *this;
}

inline double Vector::operator*(const Vector& V) const
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "Vector::operator*(): incompatible sizes\n");
#endif // DEBUG
  return cblas_ddot(n, V.v, 1, v, 1);
}

// With the intermediate classes

// scal_vec_trans

inline Vector& Vector::operator=(const __scal_vec_trans<Vector> R)
{
  cblas_dcopy(n, R.vec.v, 1, v, 1);
  cblas_dscal(n, R.alpha, v, 1);
//  std::cout<<" = __scal_vec\n";
  return *this;
}

inline Vector& Vector::operator+=(const __scal_vec_trans<Vector> R)
{
  cblas_daxpy(n, R.alpha, R.vec.v, 1, v, 1);
//  std::cout<<" += __scal_vec\n";
  return *this;
}

inline Vector& Vector::operator-=(const __scal_vec_trans<Vector> R)
{
  cblas_daxpy(n, -R.alpha, R.vec.v, 1, v, 1);
//  std::cout<<"-= __scal_vec\n";
  return *this;
}

// op_mul_vec

inline Vector& Vector::operator=(const __op_mul_vec<Matrix, Vector> R)
{
  cblas_dgemv(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, R.op.vec.r, R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 0.0, this->v, 1);
//  std::cout<<"__op_mul_vec<Matrix,Vector>\n";
  return *this;
}

inline Vector& Vector::operator+=(const __op_mul_vec<Matrix, Vector> R)
{
  cblas_dgemv(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, R.op.vec.r, R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
//  std::cout<<"+= __op_mul_vec<Matrix,Vector>\n";
  return *this;
}

inline Vector& Vector::operator-=(const __op_mul_vec<Matrix, Vector> R)
{
  cblas_dgemv(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, R.op.vec.r, R.op.vec.c, -R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
//  std::cout<<"-= __op_mul_vec<Matrix,Vector>\n";
  return *this;
}

// op_mul_vec_plus_vec, but only += and -=

inline Vector& Vector::operator+=(const __op_mul_vec_plus_vec<Matrix, Vector> R)
{
  cblas_daxpy(n, R.vecB.alpha, R.vecB.vec.v, 1, v, 1);
  cblas_dgemv(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, R.op.vec.r, R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
  std::cout << "+=__op_mul_vec_plus_vec<Matrix,Vector>\n";
  return *this;
}

inline Vector& Vector::operator-=(const __op_mul_vec_plus_vec<Matrix, Vector> R)
{
  cblas_daxpy(n, -R.vecB.alpha, R.vecB.vec.v, 1, v, 1);
  cblas_dgemv(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, R.op.vec.r, R.op.vec.c, -R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
  std::cout << "-=__op_mul_vec_plus_vec<Matrix,Vector>\n";
  return *this;
}

/// WITH RANGES

inline Vector& Vector::operator=(const __scal_vec_trans_rng<Vector>)
{
  P_ASSERT_X(false, "__scal_vec_rng\n");
  return *this;
}

inline Vector& Vector::operator=(const __op_mul_vec_rng<Matrix, Vector>)
{
  // cblas_dgemv( ... )
  P_ASSERT_X(false, "__op_mul_vec_rngMatrix,Vector>\n");
  return *this;
}

inline Vector& Vector::operator=(const __op_mul_vec_plus_vec_rng<Matrix, Vector>)
{
  P_ASSERT_X(false, "__op_mul_vec_plus_vec_rng<Matrix,Vector>\n");
  return *this;
}

// End of implementation of Vector

// Implementation of Matrix

inline Matrix& Matrix::operator=(const __scal_vec_trans<Matrix>)
{
  P_ASSERT_X(false, "__scal_vec\n");
  return *this;
}

inline Matrix& Matrix::operator=(const __op_mul_vec<Matrix, Matrix> R)
{
  cblas_dgemm(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, CblasNoTrans,
              this->r, this->c, R.op.tr == Trans ? R.op.vec.r : R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.m, R.vecA.vec.r, 0.0, this->m, this->r);
  P_ASSERT_X(false, "= __op_mul_vec<Matrix,Matrix>\n");
  return *this;
}

inline Matrix& Matrix::operator+=(const __op_mul_vec<Matrix, Matrix> R)
{
  cblas_dgemm(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, CblasNoTrans,
              this->r, this->c, R.op.tr == Trans ? R.op.vec.r : R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.m, R.vecA.vec.r, 1.0, this->m, this->r);
  P_ASSERT_X(false, "+= __op_mul_vec<Matrix,Matrix>\n");
  return *this;
}

inline Matrix& Matrix::operator-=(const __op_mul_vec<Matrix, Matrix> R)
{
  cblas_dgemm(CblasColMajor, R.op.tr == Trans ? CblasTrans : CblasNoTrans, CblasNoTrans,
              this->r, this->c, R.op.tr == Trans ? R.op.vec.r : R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.m, R.vecA.vec.r, -1.0, this->m, this->r);
  P_ASSERT_X(false, "-= __op_mul_vec<Matrix,Matrix>\n");
  return *this;
}

// obsolote : op_mul_vec_plus_vec

inline Matrix& Matrix::operator=(const __op_mul_vec_plus_vec<Matrix, Matrix>)
{
  P_ASSERT_X(false, "__op_mul_vec_plus_vec<Matrix,Matrix>\n");
  return *this;
}

/// WITH RANGES

inline Matrix& Matrix::operator=(const __scal_vec_trans_rng<Matrix>)
{
  P_ASSERT_X(false, "__scal_vec_rng\n");
  return *this;
}

inline Matrix& Matrix::operator=(const __op_mul_vec_rng<Matrix, Matrix>)
{
  // cblas_dgemm( ... );
  P_ASSERT_X(false, "__op_mul_vec_rngMatrix,Matrix>\n");
  return *this;
}

// obsolote
inline Matrix& Matrix::operator=(const __op_mul_vec_plus_vec_rng<Matrix, Matrix>)
{
  P_ASSERT_X(false, "__op_mul_vec_plus_vec_rng<Matrix,Matrix>\n");
  return *this;
}

// End of implementation of Matrix

#endif // PDDESYS_H

#endif
