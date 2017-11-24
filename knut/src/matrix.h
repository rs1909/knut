// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#ifndef MATRIX_H
#define MATRIX_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef KNUTSYS_H
#include <iostream>
#include <iterator>
#include <cstdlib>

#include "plot.h"
#include "knerror.h"
#include "cspblas.h"

extern "C"
{

#include "blas.h"
#include "laarpack.h"

}
#endif

#ifndef KNUTSYS_H

#ifndef P_ASSERT
#  define P_ASSERT(cond) do{ if(!(cond)) { std::cout<<#cond; std::cout.flush(); abort(); } }while(0)
#endif

#ifndef P_ASSERT_X
#  define P_ASSERT_X(cond,msg) do{ if(!(cond)) { std::cout<<#cond<<msg; std::cout.flush(); abort(); } }while(0)
#endif

#else

#ifndef P_ASSERT
#  define P_ASSERT(cond) do{ }while(0)
#endif

#ifndef P_ASSERT_X
#  define P_ASSERT_X(cond,msg) do{ }while(0)
#endif

#ifndef P_ASSERT_X11
#  define P_ASSERT_X11(cond, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) do{ }while(0)
#endif

#endif

template<class T> class KNArray2D;
template<class T> class KNArray3D;

template<class T>
class KNArray1D
{

  protected:

    size_t n;
    T  *v;
  private:
    const bool destructable;
  public:

    inline KNArray1D() : n(0), v(nullptr), destructable(true)
    { }

    inline KNArray1D(bool /*b*/) : n(0), v(nullptr), destructable(false)
    { }

    inline KNArray1D(size_t i) : n(i), v(new T[i]), destructable(true)
    {
      clear();
    }

    // specially for JagMatrix2D i.e. KNArray1D< KNVector >, which is indexed as (i)(j)
    inline KNArray1D(size_t i, size_t j) : n(i), v(new T[i]), destructable(true)
    {
      for (size_t r = 0; r < i; r++) v[r].init(j);
    }

    // specially for JagMatrix3D i.e. KNArray1D< KNMatrix >, which is indexed as (i)(j,k)
    inline KNArray1D(size_t i, size_t j, size_t k) : n(i), v(new T[i]), destructable(true)
    {
      for (size_t r = 0; r < i; r++) v[r].init(j, k);
    }

    inline KNArray1D(const KNArray1D<T>& V_) : n(V_.n), v(new T[V_.n]), destructable(true)
    {
      *this = V_;
    }

    inline KNArray1D(T *data, size_t i) : n(i), v(data), destructable(false)
    { }

    inline KNArray1D(const KNArray2D<T>& v, size_t i);

    inline virtual ~KNArray1D()
    {
      if (destructable) delete[] v;
    }

    inline void init(size_t i)
    {
#ifdef DEBUG
      P_ASSERT_X(destructable, "KNArray1D::init(int)");
#endif
      delete[] v;
      n = i;
      v = new T[n];
      clear();
    }

    inline void init(const KNArray1D<T>& V_)
    {
#ifdef DEBUG
      P_ASSERT_X(destructable, "KNArray1D::init(const KNArray1D<T>&)");
#endif
      delete[] v;
      n = V_.n;
      v = new T[n];
      *this = V_;
    }

    inline void init(T *data, size_t i)
    {
#ifdef DEBUG
      P_ASSERT_X((!destructable), "KNArray1D::init(T *, int)");
#endif
      n = i;
      v = data;
    }

    inline void clear()
    {
      for (size_t j = 0; j < n; j++) v[j] = T();
    } // it was T(0), but for built in types it must be the same

    inline T *pointer()
    {
      return v;
    }

    inline T *pointer(const size_t i)
    {
#ifdef DEBUG
      P_ASSERT_X(i < n, "KNArray1D::bound11&");
#endif
      return &v[i];
    }

    inline const T *pointer(const size_t i) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < n, "KNArray1D::bound11&");
#endif
      return &v[i];
    }

    inline size_t size() const
    {
      return n;
    }

    inline KNArray1D<T>& operator=(const KNArray1D<T>& V)
    {
#ifdef DEBUG
      P_ASSERT_X(n == V.n, "KNArray1D::operator=(): incompatible sizes\n");
#endif
      for (size_t i = 0; i < V.n; i++) v[i] = V.v[i];
      return *this;
    }

    inline T& operator()(size_t i)
    {
#ifdef DEBUG
      P_ASSERT_X(i < n, "KNArray1D::bound11&");
#endif
      return v[i];
    }

    inline T operator()(size_t i) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < n, "vec_bound11_\n");
#endif
      return v[i];
    }
#ifndef KNUTSYS_H
    class iterator
    {
      private:
        T* pt;
      public:
        iterator() : pt(nullptr) {}
        iterator( const KNArray1D<T>::iterator& it ) : pt(it.pt) {}
        iterator& operator= (const iterator& it) { pt = it.pt; return *this; }
        iterator operator+ ( int i ) const { iterator it(*this); it.pt += i; return it; }
        iterator operator- ( int i ) const { iterator it(*this); it.pt -= i; return it; }
        size_t    operator- (const iterator& it) const { return pt - it.pt; }
        bool      operator!= (const iterator& it) const { return pt != it.pt; }
        bool      operator== (const iterator& it) const { return pt == it.pt; }
        const T&  operator*() const { return *pt; }
        T&  operator*() { return *pt; }
        iterator& operator++() { ++pt; return *this; }
        iterator& operator--() { --pt; return *this; }
        iterator& operator++(int) { ++pt; return *this; }
        iterator& operator--(int) { --pt; return *this; }
        bool      operator< (const iterator& it) const { return pt < it.pt; }
        friend class KNArray1D<T>;
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T *;
        using reference = T &;
    };
    iterator begin() { iterator it; it.pt = 0; return it; }
    iterator end() { iterator it; it.pt = n; return it; }
#endif
};


template<class T>
class KNArray2D
{

  protected:

    T* m;
    size_t r, c;
  private:
    const bool destructable;
  public:
    inline KNArray2D() : destructable(true)
    {
      r = 0;
      c = 0;
      m = nullptr;
    }

    inline KNArray2D(size_t _r, size_t _c) : r(_r), c(_c), destructable(true)
    {
      m = new T[r*c+1];
      clear();
    }

    inline KNArray2D(const KNArray2D<T>& M) : r(M.r), c(M.c), destructable(true)
    {
      m = new T[r*c+1];
      for (size_t i = 0; i < r*c; i++) m[i] = M.m[i];
    }

    inline KNArray2D(const KNArray3D<T>& v, size_t i);

    inline virtual ~KNArray2D()
    {
      if (destructable) delete[] m;
    }
    inline size_t dim1() const { return r; }
    inline size_t dim2() const { return c; }

    inline void init(size_t _r, size_t _c)
    {
      P_ASSERT_X(destructable, "KNArray2D<T>::Init : trying to resize a non-destructable a");
      delete[] m;
      r = _r;
      c = _c;
      m = new T[r*c+1];
      clear();
    }

    inline T *pointer()
    {
      return m;
    }

    inline T *pointer(const size_t i, const size_t j)
    {
#ifdef DEBUG
      P_ASSERT_X(i < r && j < c, "bound& ");
//      P_ASSERT_X(i >= 0 && j >= 0, "lbound& ");
#endif
      return &m[i + r*j];
    }

    inline const T *pointer(const size_t i, const size_t j) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < r && j < c, "bound& ");
//      P_ASSERT_X(i >= 0 && j >= 0, "lbound& ");
#endif
      return &m[i + r*j];
    }

    inline void clear()
    {
      for (size_t i = 0; i < r*c; i++) m[i] = T(0);
    }

    inline KNArray2D<T>& operator= (const KNArray2D<T>& M)
    {
#ifdef DEBUG
      P_ASSERT_X(M.r == r && M.c == c, "KNArray2D<T>::operator= : incompatible sizes");
#endif
      for (size_t i = 0; i < r*c; i++) m[i] = M.m[i];
      return *this;
    }

    inline T& operator()(const size_t i, const size_t j)
    {
#ifdef DEBUG
      P_ASSERT_X(i < r && j < c, "bound& ");
//      P_ASSERT_X(i >= 0 && j >= 0, "lbound& ");
#endif
      return m[i + r*j];
    }

    inline T operator()(const size_t i, const size_t j) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < r && j < c, "bound_ ");
//      P_ASSERT_X(i >= 0 && j >= 0, "lbound_ ");
#endif
      return m[i + r*j];
    }

    inline T& operator()(const size_t i)
    {
#ifdef DEBUG
      P_ASSERT_X(i < r*c, "bound11&\n");
      P_ASSERT_X(r == 1 || c == 1, "H&\n");
#endif
      return m[i];
    }

    inline T operator()(const size_t i) const
    {
#ifdef DEBUG
      P_ASSERT_X(i < r*c, "bound11_\n");
      P_ASSERT_X(r == 1 || c == 1, "H&\n");
#endif
      return m[i];
    }
    friend class KNArray1D<T>;
};

template<class T>
class KNArray3D
{

  protected:

    T* m;
    size_t d1, d2, d3;

  public:
    inline KNArray3D()
    {
      d1 = d2 = d3 = 0;
      m = 0;
    }

    inline KNArray3D(size_t _d1, size_t _d2, size_t _d3) : d1(_d1), d2(_d2), d3(_d3)
    {
      if (d1*d2*d3 > 0)
      {
        m = new T[d1*d2*d3+1];
        clear();
      } else m = nullptr;
    }

    inline KNArray3D(const KNArray3D<T>& M) : d1(M.d1), d2(M.d2), d3(M.d3)
    {
      if (d1*d2*d3 > 0)
      {
        m = new T[d1*d2*d3+1];
        for (size_t i = 0; i < d1*d2*d3; i++) m[i] = M.m[i];
      } else m = nullptr;
    }

    inline virtual ~KNArray3D()
    {
      delete[] m;
    }

    inline size_t dim1() const { return d1; }
    inline size_t dim2() const { return d2; }
    inline size_t dim3() const { return d3; }
    inline void init(size_t _d1, size_t _d2, size_t _d3)
    {
      delete[] m;
      d1 = _d1;
      d2 = _d2;
      d3 = _d3;
      m = new T[d1*d2*d3+1];
      clear();
    }

    inline void clear()
    {
      for (size_t i = 0; i < d1*d2*d3; i++) m[i] = 0;
    }
    inline T *pointer(const size_t i, const size_t j, const size_t k)
    {
#ifdef DEBUG
      P_ASSERT_X11((i < d1) && (j < d2) && (k < d3), d1, "|",i , " ", d2, "|", j, " ", d3, "|", k);
//      P_ASSERT_X((i >= 0) && (j >= 0) && (k >= 0), "lbound\n");
#endif
      return &m[i + d1*(j + d2*k)];
    }

    inline KNArray3D<T>& operator= (const KNArray3D<T>& M)
    {
#ifdef DEBUG
      P_ASSERT_X((M.d1 == d1) && (M.d2 == d2) && (M.d3 == d3), "KNArray3D<T>::operator= : incompatible sizes\n");
#endif
      for (size_t i = 0; i < d1*d2*d3; i++) m[i] = M.m[i];
      return *this;
    }

    inline T& operator()(const size_t i, const size_t j, const size_t k)
    {
#ifdef DEBUG
      P_ASSERT_X((i < d1) && (j < d2) && (k < d3), "bound&\n");
//      P_ASSERT_X((i >= 0) && (j >= 0) && (k >= 0), "lbound&\n");
#endif
      return m[i + d1*(j + d2*k)];
    }

    inline T operator()(const size_t i, const size_t j, const size_t k) const
    {
#ifdef DEBUG
      P_ASSERT_X11((i < d1) && (j < d2) && (k < d3), d1, "|",i , " ", d2, "|", j, " ", d3, "|", k);
//      P_ASSERT_X((i >= 0) && (j >= 0) && (k >= 0), "lbound\n");
#endif
      return m[i + d1*(j + d2*k)];
    }
    friend class KNArray2D<T>;
};

template<class T> inline KNArray1D<T>::KNArray1D(const KNArray2D<T>& vec, size_t i) : destructable(false)
{
#ifdef DEBUG
  P_ASSERT_X(i < vec.c, "bound\n");
#endif
  v = &vec.m[vec.r*i];
  n = vec.r;
}

template<class T> inline KNArray2D<T>::KNArray2D(const KNArray3D<T>& vec, size_t i) : destructable(false)
{
#ifdef DEBUG
  P_ASSERT_X(i < vec.d3, "bound\n");
#endif
  m = &vec.m[vec.d1*vec.d2*i];
  r = vec.d1;
  c = vec.d2;
}

#ifndef KNUTSYS_H

class KNSparseMatrix;
class KNMatrix;
class KNVector;

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
    inline rng(size_t i1_, size_t i2_) : i1(i1_), i2(i2_), j1(0), j2(0)
    { }
    inline rng(size_t i1_, size_t i2_, size_t j1_, size_t j2_) : i1(i1_), i2(i2_), j1(j1_), j2(j2_)
    { }

    const size_t i1, i2, j1, j2;
};

/// with range
template< class VT > class __scal_vec_trans_rng : public rng
{
  public:
    inline __scal_vec_trans_rng(const VT& v)                                        : rng(0, v.row(), 0, v.col()), vec(v), alpha(1.0), tr(NoTrans)
    { }
    inline __scal_vec_trans_rng(const VT& v, double a)                              : rng(0, v.row(), 0, v.col()), vec(v), alpha(a),   tr(NoTrans)
    { }
    inline __scal_vec_trans_rng(const VT& v, double a, enum cspblas_Trans t)        : rng(0, v.row(), 0, v.col()), vec(v), alpha(a),   tr(t)
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
    inline __op_mul_vec_rng< VT, KNVector > operator*(__scal_vec_trans_rng<KNVector> v);
    inline __op_mul_vec_rng< VT, KNMatrix > operator*(__scal_vec_trans_rng<KNMatrix> v);
    inline __op_mul_vec_rng< VT, KNVector > operator*(const KNVector& v);
    inline __op_mul_vec_rng< VT, KNMatrix > operator*(const KNMatrix& v);

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
    inline __op_mul_vec< VT, KNVector >                    operator*(const KNVector& v);
    inline __op_mul_vec< VT, KNMatrix >                    operator*(const KNMatrix& v);
    inline __op_mul_vec_rng< VT, KNVector >                operator*(__scal_vec_trans_rng<KNVector> v);
    inline __op_mul_vec_rng< VT, KNMatrix >                operator*(__scal_vec_trans_rng<KNMatrix> v);

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

#endif // KNUTSYS_H

class KNVector : public KNArray1D<double>
{

  public:

    inline KNVector()
    { }

    inline KNVector(size_t i) : KNArray1D<double>(i)
    { }

    inline KNVector(bool b) : KNArray1D<double>(b)
    { }

    inline KNVector(const KNVector& V) : KNArray1D<double>(V)
    { }

    inline KNVector(const KNArray2D<double>& m, size_t i) : KNArray1D<double>(m,i)
    { }

    inline ~KNVector() override
    { }

#ifndef KNUTSYS_H

    inline void print(std::ostream& os) const
    {
      for (size_t j = 0; j < n; j++) os << v[j] << '\t';
      os << '\n';
    }

    inline void random()
    {
      static blasint idist = 2;
      static blasint iseed[4] =
        {
          1, 3, 5, 7
        };
      blasint N = static_cast<blasint>(n);
      knut_dlarnv(&idist, iseed, &N, v);
    }

    inline KNVector& operator= (const KNVector& V);
    inline KNVector& operator+=(const KNVector& V);
    inline KNVector& operator-=(const KNVector& V);
    inline KNVector& operator/=(double div);
    inline KNVector& operator*=(double mul);
    inline double  operator*(const KNVector& V) const;

    // KNMatrix operations
    inline KNVector& operator= (const __scal_vec_trans<KNVector>);
    inline KNVector& operator+=(const __scal_vec_trans<KNVector>);
    inline KNVector& operator-=(const __scal_vec_trans<KNVector>);
    inline KNVector& operator= (const __op_mul_vec<KNMatrix, KNVector>);
    inline KNVector& operator+=(const __op_mul_vec<KNMatrix, KNVector>);
    inline KNVector& operator-=(const __op_mul_vec<KNMatrix, KNVector>);
    inline KNVector& operator+=(const __op_mul_vec_plus_vec<KNMatrix, KNVector>);
    inline KNVector& operator-=(const __op_mul_vec_plus_vec<KNMatrix, KNVector>);
    // for KNSparseMatrix
    inline KNVector& operator= (const __op_mul_vec<KNSparseMatrix, KNVector>);
    // obsolote
    inline KNVector& operator= (const __op_mul_vec_plus_vec<KNSparseMatrix, KNVector>);
    /// with ranges
    inline KNVector& operator=(const __scal_vec_trans_rng<KNVector>);
    inline KNVector& operator=(const __op_mul_vec_rng<KNMatrix, KNVector>);
    // obsolote
    inline KNVector& operator=(const __op_mul_vec_plus_vec_rng<KNMatrix, KNVector>);
    // for KNSparseMatrix
    inline KNVector& operator=(const __op_mul_vec_rng<KNSparseMatrix, KNVector>);
    // obsolote
    inline KNVector& operator=(const __op_mul_vec_plus_vec_rng<KNSparseMatrix, KNVector>);

    inline __scal_vec_trans<KNVector>     operator-() const
    {
      return __scal_vec_trans<KNVector>(*this, -1.0);
    }
    inline __scal_vec_trans_rng<KNVector> operator[ ](const rng r) const
    {
      return __scal_vec_trans_rng<KNVector>(*this, r);
    }

#endif // KNUTSYS_H

    friend class KNSparseMatrix;
    friend class KNLuSparseMatrix;
    friend class KNMatrix;
    friend class KNLuMatrix;

};

class KNMatrix : public KNArray2D<double>
{
  public:

    inline KNMatrix()
    { }

    inline KNMatrix(size_t i, size_t j) : KNArray2D<double>(i, j)
    { }

    inline KNMatrix(const KNMatrix& M) : KNArray2D<double>(M)
    { }

    inline KNMatrix(const KNArray3D<double>& m, size_t i) : KNArray2D<double>(m,i)
    { }

    inline ~KNMatrix() override
    { }

    inline size_t row() const
    {
      return r;
    }
    inline size_t col() const
    {
      return c;
    }
    inline size_t size() const
    {
#ifdef DEBUG
      P_ASSERT_X((r == 1) || (c == 1), "KNMatrix::size(): not a single row or column.\n");
#endif // DEBUG
      return r*c;
    }

#ifndef KNUTSYS_H

    inline KNMatrix& operator= (const __scal_vec_trans<KNMatrix>);
    inline KNMatrix& operator= (const __op_mul_vec<KNMatrix, KNMatrix>);
    inline KNMatrix& operator+=(const __op_mul_vec<KNMatrix, KNMatrix>);
    inline KNMatrix& operator-=(const __op_mul_vec<KNMatrix, KNMatrix>);
    // obsolote
    inline KNMatrix& operator=(const __op_mul_vec_plus_vec<KNMatrix, KNMatrix>);
    // with KNSparseMatrix
    inline KNMatrix& operator=(const __op_mul_vec<KNSparseMatrix, KNMatrix>);
    // obsolote
    inline KNMatrix& operator=(const __op_mul_vec_plus_vec<KNSparseMatrix, KNMatrix>);
    /// with ranges
    inline KNMatrix& operator=(const __scal_vec_trans_rng<KNMatrix>);
    inline KNMatrix& operator=(const __op_mul_vec_rng<KNMatrix, KNMatrix>);
    // obsolote
    inline KNMatrix& operator=(const __op_mul_vec_plus_vec_rng<KNMatrix, KNMatrix>);
    // with KNSparseMatrix
    inline KNMatrix& operator=(const __op_mul_vec_rng<KNSparseMatrix, KNMatrix>);
    // obsolote
    inline KNMatrix& operator=(const __op_mul_vec_plus_vec_rng<KNSparseMatrix, KNMatrix>);

    void eigenvalues(KNVector& re, KNVector& im);
    void eigenvalues(KNVector& re, KNVector& im, KNMatrix& lev, KNMatrix& rev);
    void sparsityPlot(GnuPlot& pl);

    /* operators */

    inline __op_mul_vec<KNMatrix, KNVector>     operator*(const KNVector& v) const
    {
      return __op_mul_vec<KNMatrix, KNVector>(__scal_vec_trans<KNMatrix>(*this), __scal_vec_trans<KNVector>(v));
    }
    inline __op_mul_vec<KNMatrix, KNMatrix>     operator*(const KNMatrix& v) const
    {
      return __op_mul_vec<KNMatrix, KNMatrix>(__scal_vec_trans<KNMatrix>(*this), __scal_vec_trans<KNMatrix>(v));
    }
    inline __op_mul_vec_rng<KNMatrix, KNVector> operator*(const __scal_vec_trans_rng<KNVector> v) const
    {
      return __op_mul_vec_rng<KNMatrix, KNVector>(__scal_vec_trans_rng<KNMatrix>(*this), v);
    }
    inline __op_mul_vec_rng<KNMatrix, KNMatrix> operator*(const __scal_vec_trans_rng<KNMatrix> v) const
    {
      return __op_mul_vec_rng<KNMatrix, KNMatrix>(__scal_vec_trans_rng<KNMatrix>(*this), v);
    }
    inline __scal_vec_trans<KNMatrix>        operator!() const
    {
      return __scal_vec_trans<KNMatrix>(*this, 1.0, Trans);
    }
    inline __scal_vec_trans_rng<KNMatrix>    operator[ ](const rng r_) const
    {
      return __scal_vec_trans_rng<KNMatrix>(*this, r_);
    }

    inline void print(std::ostream& os)
    {
      double sum = 0.0;
      for (size_t i = 0; i < r; i++)
      {
        for (size_t j = 0; j < c; j++)
        {
          os << (*this)(i, j) << '\t';
          sum += (*this)(i, j);
        }
        os << '\n';
      }
      os << "SUM: " << sum << '\n';
    }

#endif // KNUTSYS_H

    friend class KNVector;
    friend class KNLuMatrix;
    friend class KNSparseMatrix;
    friend class KNLuSparseMatrix;
};

#ifndef KNUTSYS_H

class KNLuMatrix : public KNMatrix
{

    bool     fact; // == true if factorized
    // for factorizig
    double*  mf;
    blasint* ipiv;      // --
    double*  work;
    blasint* iwork;     // --
    blasint  info;      // --
    // for solving
    double rcond;
    double ferr;  // in case of KNVector these are just double
    double berr;

  public:

    inline KNLuMatrix(size_t _r, size_t _c) : KNMatrix(_r, _c)
    {
      fact = false;
      mf = new double[r*c+1];
      ipiv = new blasint[this->r+1];
      work = new double[4*(this->r)];
      iwork = new blasint[this->r];
    }

    inline KNLuMatrix(const KNMatrix& M) : KNMatrix(M)
    {
      fact = false;
      mf = new double[r*c+1];
      ipiv = new blasint[this->r+1];
      work = new double[4*(this->r)];
      iwork = new blasint[this->r];
    }

    inline ~KNLuMatrix() override
    {
      delete[] iwork;
      delete[] work;
      delete[] ipiv;
      delete[] mf;
    }

    inline KNLuMatrix& operator=(KNLuMatrix& M)
    {
      KNMatrix::operator=(M);
      std::cout << "Copying KNLuMatrix is not yet implemented, though the KNMatrix part will be copied normally\n";
      return *this;
    }

    inline void modified()
    {
      fact = false;
    }

    inline size_t  row() const
    {
      return KNMatrix::row();
    }
    inline size_t  col() const
    {
      return KNMatrix::col();
    }

    void luFactorize();
    void solve(KNVector& X, const KNVector& B, bool TRANS = false);
    void solve(KNMatrix& X, const KNMatrix& B, bool TRANS = false);
};

using JagMatrix3D = KNArray1D<KNMatrix>;
using JagVector2D = KNArray1D<KNVector>;

/// specialized versions of the Clear function
template< > inline void KNArray1D< KNArray1D<int> >::clear()
{
  for (size_t i = 0; i < n; i++) v[i].clear();
}
template< > inline void KNArray1D< KNArray1D<double> >::clear()
{
  for (size_t i = 0; i < n; i++) v[i].clear();
}
template< > inline void KNArray1D< KNVector >::clear()
{
  for (size_t i = 0; i < n; i++) v[i].clear();
}
template< > inline void KNArray1D< KNArray2D<int> >::clear()
{
  for (size_t i = 0; i < n; i++) v[i].clear();
}
template< > inline void KNArray1D< KNArray2D<double> >::clear()
{
  for (size_t i = 0; i < n; i++) v[i].clear();
}
template< > inline void KNArray1D< KNMatrix >::clear()
{
  for (size_t i = 0; i < n; i++) v[i].clear();
}

template< > inline KNArray1D< KNVector >::KNArray1D(const KNArray1D< KNVector >& V_) : destructable(true)
{
  n = V_.n;
  v = new KNVector[V_.n];
  for (size_t i = 0; i < V_.n; i++) v[i].init(V_.v[i]);
}

/// Member functions and Operators for __scal_vec_trans_rng
template<> inline __scal_vec_trans_rng<KNVector>::__scal_vec_trans_rng(const KNVector& v)
    : rng(0, v.size()), vec(v), alpha(1.0), tr(NoTrans)
{ }

template<> inline __scal_vec_trans_rng<KNVector>::__scal_vec_trans_rng(const KNVector& v, double a)
    : rng(0, v.size()), vec(v), alpha(a),   tr(NoTrans)
{ }

template<> inline __scal_vec_trans_rng<KNVector>::__scal_vec_trans_rng(const KNVector& v, double a, enum cspblas_Trans t)
    : rng(0, v.size()), vec(v), alpha(a),   tr(t)
{ }

template< class VT > inline __scal_vec_trans_rng<VT>::__scal_vec_trans_rng(__scal_vec_trans<VT> v)
    : rng(0, v.vec.row(), 0, v.vec.col()), vec(v.vec), alpha(v.alpha), tr(v.tr)
{ }

template< >          inline __scal_vec_trans_rng<KNVector>::__scal_vec_trans_rng(__scal_vec_trans<KNVector> v)
    : rng(0, v.vec.size()), vec(v.vec), alpha(v.alpha), tr(v.tr)
{ }

template< class VT > inline __op_mul_vec_rng< VT, KNVector > __scal_vec_trans_rng<VT>::operator*(__scal_vec_trans_rng<KNVector> v)
{
  return __op_mul_vec_rng< VT, KNVector >(*this, v);
}

template< class VT > inline __op_mul_vec_rng< VT, KNMatrix > __scal_vec_trans_rng<VT>::operator*(__scal_vec_trans_rng<KNMatrix> v)
{
  return __op_mul_vec_rng< VT, KNMatrix >(*this, v);
}

template< class VT > inline __op_mul_vec_rng< VT, KNVector > __scal_vec_trans_rng<VT>::operator*(const KNVector& v)
{
  return __op_mul_vec_rng< VT, KNVector >(*this, __scal_vec_trans_rng<KNVector>(v));
}

template< class VT > inline __op_mul_vec_rng< VT, KNMatrix > __scal_vec_trans_rng<VT>::operator*(const KNMatrix& v)
{
  return __op_mul_vec_rng< VT, KNMatrix >(*this, __scal_vec_trans_rng<KNMatrix>(v));
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
template< class VT > inline __op_mul_vec< VT, KNVector >     __scal_vec_trans<VT>::operator*(const KNVector& v)
{
  return __op_mul_vec< VT, KNVector >(*this, __scal_vec_trans<KNVector>(v));
}

template< class VT > inline __op_mul_vec< VT, KNMatrix >     __scal_vec_trans<VT>::operator*(const KNMatrix& v)
{
  return __op_mul_vec< VT, KNMatrix >(*this, __scal_vec_trans<KNMatrix>(v));
}

template< class VT > inline __op_mul_vec_rng< VT, KNVector > __scal_vec_trans<VT>::operator*(__scal_vec_trans_rng<KNVector> v)
{
  return __op_mul_vec_rng< VT, KNVector >(__scal_vec_trans_rng<VT>(*this), v);
}

template< class VT > inline __op_mul_vec_rng< VT, KNMatrix > __scal_vec_trans<VT>::operator*(__scal_vec_trans_rng<KNMatrix> v)
{
  return __op_mul_vec_rng< VT, KNMatrix >(__scal_vec_trans_rng<VT>(*this), v);
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
inline __scal_vec_trans < KNVector >                      operator*(double a, const KNVector& v)
{
  return __scal_vec_trans<KNVector>(v, a);
}

inline __scal_vec_trans < KNMatrix >                      operator*(double a, const KNMatrix& v)
{
  return __scal_vec_trans<KNMatrix>(v, a);
}

inline __scal_vec_trans < KNSparseMatrix >                    operator*(double a, const KNSparseMatrix& v)
{
  return __scal_vec_trans<KNSparseMatrix>(v, a);
}

template< class VT > inline __scal_vec_trans < VT  >    operator*(double a, __scal_vec_trans< VT > v)
{
  return __scal_vec_trans<VT>(v.vec, a * v.alpha, v.tr);
}

template< class VT > inline __scal_vec_trans_rng < VT > operator*(double a, __scal_vec_trans_rng< VT > v)
{
  return __scal_vec_trans_rng<VT>(v.vec, v, a * v.alpha, v.tr);
}

inline __scal_vec_trans < KNVector >                      operator-(const KNVector& v)
{
  return __scal_vec_trans<KNVector>(v, -1.0);
}

inline __scal_vec_trans < KNMatrix >                      operator-(const KNMatrix& v)
{
  return __scal_vec_trans<KNMatrix>(v, -1.0);
}

inline __scal_vec_trans < KNSparseMatrix >                    operator-(const KNSparseMatrix& v)
{
  return __scal_vec_trans<KNSparseMatrix>(v, -1.0);
}

// Implementation of KNVector

inline KNVector& KNVector::operator=(const KNVector& V)
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "KNVector::operator=(): incompatible sizes\n");
#endif // DEBUG
  BLAS_dcopy(n, V.v, 1, v, 1);
  return *this;
}

inline KNVector& KNVector::operator+=(const KNVector& V)
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "KNVector::operator+=(): incompatible sizes\n");
#endif // DEBUG
  BLAS_daxpy(n, 1.0, V.v, 1, v, 1);
  return *this;
}

inline KNVector& KNVector::operator-=(const KNVector& V)
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "KNVector::operator-=(): incompatible sizes\n");
#endif // DEBUG
  BLAS_daxpy(n, -1.0, V.v, 1, v, 1);
  return *this;
}

inline KNVector& KNVector::operator/=(double div)
{
  BLAS_dscal(n, 1.0 / div, v, 1);
  return *this;
}

inline KNVector& KNVector::operator*=(double mul)
{
  BLAS_dscal(n, mul, v, 1);
  return *this;
}

inline double KNVector::operator*(const KNVector& V) const
{
#ifdef DEBUG
  P_ASSERT_X(n == V.n, "KNVector::operator*(): incompatible sizes\n");
#endif // DEBUG
  return BLAS_ddot(n, V.v, 1, v, 1);
}

// With the intermediate classes

// scal_vec_trans

inline KNVector& KNVector::operator=(const __scal_vec_trans<KNVector> R)
{
  BLAS_dcopy(n, R.vec.v, 1, v, 1);
  BLAS_dscal(n, R.alpha, v, 1);
//  std::cout<<" = __scal_vec\n";
  return *this;
}

inline KNVector& KNVector::operator+=(const __scal_vec_trans<KNVector> R)
{
  BLAS_daxpy(n, R.alpha, R.vec.v, 1, v, 1);
//  std::cout<<" += __scal_vec\n";
  return *this;
}

inline KNVector& KNVector::operator-=(const __scal_vec_trans<KNVector> R)
{
  BLAS_daxpy(n, -R.alpha, R.vec.v, 1, v, 1);
//  std::cout<<"-= __scal_vec\n";
  return *this;
}

// op_mul_vec

inline KNVector& KNVector::operator=(const __op_mul_vec<KNMatrix, KNVector> R)
{
  BLAS_dgemv(R.op.tr == Trans ? 'T' : 'N', R.op.vec.r, R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 0.0, this->v, 1);
//  std::cout<<"__op_mul_vec<KNMatrix,KNVector>\n";
  return *this;
}

inline KNVector& KNVector::operator+=(const __op_mul_vec<KNMatrix, KNVector> R)
{
  BLAS_dgemv(R.op.tr == Trans ? 'T' : 'N', R.op.vec.r, R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
//  std::cout<<"+= __op_mul_vec<KNMatrix,KNVector>\n";
  return *this;
}

inline KNVector& KNVector::operator-=(const __op_mul_vec<KNMatrix, KNVector> R)
{
  BLAS_dgemv(R.op.tr == Trans ? 'T' : 'N', R.op.vec.r, R.op.vec.c, -R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
//  std::cout<<"-= __op_mul_vec<KNMatrix,KNVector>\n";
  return *this;
}

// op_mul_vec_plus_vec, but only += and -=

inline KNVector& KNVector::operator+=(const __op_mul_vec_plus_vec<KNMatrix, KNVector> R)
{
  BLAS_daxpy(n, R.vecB.alpha, R.vecB.vec.v, 1, v, 1);
  BLAS_dgemv(R.op.tr == Trans ? 'T' : 'N', R.op.vec.r, R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
  std::cout << "+=__op_mul_vec_plus_vec<KNMatrix,KNVector>\n";
  return *this;
}

inline KNVector& KNVector::operator-=(const __op_mul_vec_plus_vec<KNMatrix, KNVector> R)
{
  BLAS_daxpy(n, -R.vecB.alpha, R.vecB.vec.v, 1, v, 1);
  BLAS_dgemv(R.op.tr == Trans ? 'T' : 'N', R.op.vec.r, R.op.vec.c, -R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.v, 1, 1.0, this->v, 1);
  std::cout << "-=__op_mul_vec_plus_vec<KNMatrix,KNVector>\n";
  return *this;
}

/// WITH RANGES

inline KNVector& KNVector::operator=(const __scal_vec_trans_rng<KNVector>)
{
  P_ASSERT_X(false, "__scal_vec_rng\n");
  return *this;
}

inline KNVector& KNVector::operator=(const __op_mul_vec_rng<KNMatrix, KNVector>)
{
  // cblas_dgemv( ... )
  P_ASSERT_X(false, "__op_mul_vec_rngMatrix,KNVector>\n");
  return *this;
}

inline KNVector& KNVector::operator=(const __op_mul_vec_plus_vec_rng<KNMatrix, KNVector>)
{
  P_ASSERT_X(false, "__op_mul_vec_plus_vec_rng<KNMatrix,KNVector>\n");
  return *this;
}

// End of implementation of KNVector

// Implementation of KNMatrix

inline KNMatrix& KNMatrix::operator=(const __scal_vec_trans<KNMatrix>)
{
  P_ASSERT_X(false, "__scal_vec\n");
  return *this;
}

inline KNMatrix& KNMatrix::operator=(const __op_mul_vec<KNMatrix, KNMatrix> R)
{
  BLAS_dgemm(R.op.tr == Trans ? 'T' : 'N', 'N',
              this->r, this->c, R.op.tr == Trans ? R.op.vec.r : R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.m, R.vecA.vec.r, 0.0, this->m, this->r);
  P_ASSERT_X(false, "= __op_mul_vec<KNMatrix,Matrix>\n");
  return *this;
}

inline KNMatrix& KNMatrix::operator+=(const __op_mul_vec<KNMatrix, KNMatrix> R)
{
  BLAS_dgemm(R.op.tr == Trans ? 'T' : 'N', 'N',
              this->r, this->c, R.op.tr == Trans ? R.op.vec.r : R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.m, R.vecA.vec.r, 1.0, this->m, this->r);
  P_ASSERT_X(false, "+= __op_mul_vec<KNMatrix,Matrix>\n");
  return *this;
}

inline KNMatrix& KNMatrix::operator-=(const __op_mul_vec<KNMatrix, KNMatrix> R)
{
  BLAS_dgemm(R.op.tr == Trans ? 'T' : 'N', 'N',
              this->r, this->c, R.op.tr == Trans ? R.op.vec.r : R.op.vec.c, R.op.alpha, R.op.vec.m, R.op.vec.r,
              R.vecA.vec.m, R.vecA.vec.r, -1.0, this->m, this->r);
  P_ASSERT_X(false, "-= __op_mul_vec<KNMatrix,Matrix>\n");
  return *this;
}

// obsolote : op_mul_vec_plus_vec

inline KNMatrix& KNMatrix::operator=(const __op_mul_vec_plus_vec<KNMatrix, KNMatrix>)
{
  P_ASSERT_X(false, "__op_mul_vec_plus_vec<KNMatrix,Matrix>\n");
  return *this;
}

/// WITH RANGES

inline KNMatrix& KNMatrix::operator=(const __scal_vec_trans_rng<KNMatrix>)
{
  P_ASSERT_X(false, "__scal_vec_rng\n");
  return *this;
}

inline KNMatrix& KNMatrix::operator=(const __op_mul_vec_rng<KNMatrix, KNMatrix>)
{
  // cblas_dgemm( ... );
  P_ASSERT_X(false, "__op_mul_vec_rngMatrix,KNMatrix>\n");
  return *this;
}

// obsolote
inline KNMatrix& KNMatrix::operator=(const __op_mul_vec_plus_vec_rng<KNMatrix, KNMatrix>)
{
  P_ASSERT_X(false, "__op_mul_vec_plus_vec_rng<KNMatrix,Matrix>\n");
  return *this;
}

// End of implementation of KNMatrix

#endif // KNUTSYS_H

#endif
