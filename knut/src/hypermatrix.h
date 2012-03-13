// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef HYPERMATRIX_H
#define HYPERMATRIX_H

#include "matrix.h"
#include "spmatrix.h"
#include "plot.h"
#include <cmath>

class KNBlockVector
{
  public:

    inline KNBlockVector(size_t i, size_t, size_t k) : V1(i), V3(k)
    { }
    inline KNBlockVector(const KNBlockVector& hv) : V1(hv.V1), V3(hv.V3)
    { }
    inline ~KNBlockVector()
    { }

    inline KNVector& getV1()
    {
      return V1;
    }
    inline const KNVector& getV1() const
    {
      return V1;
    }
    inline KNVector& getV3()
    {
      return V3;
    }
    inline const KNVector& getV3() const
    {
      return V3;
    }

  private:

    KNVector V1;
    KNVector V3;
};

// #define FACT KNLuSparseMatrix
template< class FACT > class KNBlockMatrix
{
  public:

    // constructing Sparse-Dense type hypermatrix
    KNBlockMatrix(size_t i, size_t, size_t k, size_t nz);
    ~KNBlockMatrix();

    inline FACT&        getA11()
    {
      return *A11;
    }
    inline JagVector2D& getA13()
    {
      return *A13;
    }
    inline KNVector&      getA13(size_t i)
    {
      return (*A13)(i);
    }
    inline JagVector2D& getA31()
    {
      return *A31;
    }
    inline KNVector&      getA31(size_t i)
    {
      return (*A31)(i);
    }
    inline KNMatrix&      getA33()
    {
      return *A33;
    }
    inline double&      getA33(size_t i, size_t j)
    {
      return (*A33)(i, j);
    }

    template< class T, bool trans > inline void __BEM
    (
      T& _A, KNVector& _b, KNVector& _bStar, double& _d,
      KNVector& x, double& y, const KNVector& f, const double& g,
      KNVector& bem_v, KNVector& bem_vStar, KNVector& bem_w, KNVector& bem_f1
    );
    template< class T, bool DTR > inline void __BEMWS
    (
      T&           A,  JagVector2D& B,      JagVector2D&  BStar, KNMatrix&       D,
      JagVector2D& X,  KNVector&      Y,      JagVector2D&  F,     KNVector&       G,
      size_t j, const JagVector2D& V, const JagVector2D& VStar, const KNVector& delta, const KNVector& deltaStar, bool trans
    );
    template< class T, bool DTR > inline void __BEMWF
    (
      size_t bord,
      T&           A,  JagVector2D& B,      JagVector2D&  BStar, KNMatrix& D,
      JagVector2D& X,  KNVector&      Y,      JagVector2D&  F,     KNVector& G,
      JagVector2D& V,  JagVector2D& VStar,  KNVector&       delta, KNVector& deltaStar
    );
    template< class T, bool trans > inline void __BEMW
    (
      size_t bord,
      T&           A,  JagVector2D& B,      JagVector2D&  BStar, KNMatrix&       D,
      JagVector2D& X,  KNVector&      Y,      JagVector2D&  F,     KNVector&       G,
      JagVector2D& V,  JagVector2D& VStar,  KNVector&       delta, KNVector&       deltaStar,
      KNVector&      x,  KNVector&      y,      const KNVector& f,     const KNVector& g
    );

    void solve(KNBlockVector& X, const KNBlockVector& F);

    void solve(KNBlockVector& X, const KNBlockVector& F, size_t bord);

    template<bool trans> void multiply(KNBlockVector& X, const KNBlockVector& F, size_t bord);

    void check(const KNBlockVector& X, const KNBlockVector& F, size_t bord);

    void solveDirect(KNBlockVector& X, const KNBlockVector& F);

    void solve(KNVector& x, const KNVector& f)
    {
      A11->solve(x, f);
    }

    template<bool trans> void multiply(KNVector& R1, double& R3, const KNVector& X1, const double& X3);

    template<bool trans> void check(const KNVector& x, const double& z, const KNVector& f, const double& h);

    void solve(KNVector& x, double& z, const KNVector& f, const double& h);   // BEM

    void solveTr(KNVector& x, double& z, const KNVector& f, const double& h);   // BEM

    template<bool trans> void multiply(size_t bord, KNVector& R1, KNVector& R3, const KNVector& X1, const KNVector& X3);

    template<bool trans> void check(size_t bord, const KNVector& x, const KNVector& z, const KNVector& f, const KNVector& h);

    void solve(size_t bord, KNVector& x, KNVector& z, const KNVector& f, const KNVector& h);   // BEMW

    void solveTr(size_t bord, KNVector& x, KNVector& z, const KNVector& f, const KNVector& h);   // BEMW

  protected:

    // matrix blocks

    FACT*    A11;
    JagVector2D* A13;
    JagVector2D* A31;
    KNMatrix*   A33;

  private:
    // temporary storage for factorization

    // BEMW 1.
    JagVector2D* VV1;
    JagVector2D* VV1Star;
    JagVector2D* FF1;
    KNVector*   GG1;
    JagVector2D* XX1;
    KNVector*   YY1;
    KNVector*   delta1;
    KNVector*   delta1Star;

    // BEM temporaries
    KNVector         bem_v_1;
    KNVector         bem_vStar_1;
    KNVector         bem_w_1;
    KNVector         bem_f1_1;
    KNVector         bem_v_2;
    KNVector         bem_vStar_2;
    KNVector         bem_w_2;
    KNVector         bem_f1_2;
};

typedef KNBlockMatrix<KNLuSparseMatrix> KNSparseBlockMatrix;

///--------------------------------
///
/// Now the implementation
///
///--------------------------------

// constructing Sparse-Dense type hypermatrix
template<class FACT> KNBlockMatrix<FACT> :: KNBlockMatrix(size_t i, size_t j, size_t k, size_t nz) :
    bem_v_1(i),
    bem_vStar_1(i),
    bem_w_1(i),
    bem_f1_1(i),
    bem_v_2(j),
    bem_vStar_2(j),
    bem_w_2(j),
    bem_f1_2(j)
{
  // constructing the data members
  if (i != 0) A11 = new FACT('R', i, nz);
  else A11 = 0;
  if (k != 0) A33 = new KNMatrix(k, k);
  else A33 = 0;

  if ((i != 0) && (k != 0))
  {
    A13 = new JagVector2D(k, i);
    A31 = new JagVector2D(k, i);

    VV1 = new JagVector2D(k + 1);
    VV1Star = new JagVector2D(k + 1);
    FF1 = new JagVector2D(k + 1);
    GG1 = new KNVector(k + 1);
    XX1 = new JagVector2D(k + 1);
    YY1 = new KNVector(k + 1);
    delta1 = new KNVector(k + 1);
    delta1Star = new KNVector(k + 1);

    for (size_t r = 0; r < k + 1; r++)
    {
      (*VV1)(r).init(i + r);
      (*VV1Star)(r).init(i + r);
      (*XX1)(r).init(i + r);
      (*FF1)(r).init(i + r);
    }
  }
  else
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
}

template<class FACT> KNBlockMatrix<FACT> :: ~KNBlockMatrix()
{
  // deleting the data members
  delete A11;
  delete A13;
  delete A31;
  delete A33;

  delete VV1;
  delete VV1Star;
  delete FF1;
  delete GG1;
  delete XX1;
  delete YY1;
  delete delta1;
  delete delta1Star;
}

//  BEM temporaries
//  KNVector v( _b.size() ); == A.size()                     bem_v
//  KNVector vStar( _b.size() );                             bem_vStar
//  KNVector w( _b.size() );                                 bem_w
//  KNVector f1( _b.size() );                                bem_f1;

//  GMBE temporaries
//  KNVector xi( _A32.size() ); == dim2                      gmbe_xi;
//  KNVector cc( _A31 ); == dim1                             gmbe_cc
//  KNVector fr( F2 ); == dim2                               cm_fr

// GMBEW temporaries
//  JagVector2D xi( bord, _A32(0).size() );  -> different  gmbew_xi
//  JagVector2D cc( _A31 );                  -> different  gmbew_cc
//  KNMatrix dd( bord, bord );                               gmbew_dd
//  KNVector gg( bord );                                     gmbew_gg
//  KNVector fr( F2.size() ); -> same as GMBE                cm_fr
//  KNVector gr( bord );                                     gmbew_gr
//  KNVector yy( bord );                                     gmbew_yy
//  KNMatrix m_beta( bord, bord );                           gmbew_m_beta

// BEM
// swap b and bStar if trans==true
template<class FACT> template<class T, bool trans>
inline void KNBlockMatrix<FACT> :: __BEM(T& _A, KNVector& _b, KNVector& _bStar, double& _d,
                                    KNVector& x, double& y, const KNVector& f, const double& g,
                                    KNVector& bem_v, KNVector& bem_vStar, KNVector& bem_w, KNVector& bem_f1)
{

  double delta, deltaStar;
  double g1;
  double y1, y2;

  // Doolittle
  _A.solve(bem_vStar, _bStar, !trans);                              // Step 1.
  deltaStar = _d - bem_vStar * _b;                                  // Step 2. Scalar + ddot
  // Crout
  _A.solve(bem_v, _b, trans);                                       // Step 3.
  delta = _d - _bStar * bem_v;                                      // Step 4. Scalar + ddot
  // approx Y
  y1 = (g - bem_vStar * f) / deltaStar;                             // Step 5. Scalar + ddot
  // residuals
  bem_f1 = f;
  bem_f1 -= y1 * _b;                                                // Step 6. daxpy
  g1 = g - _d * y1;                                                 // Step 7. Scalar

  // residual corrections
  _A.solve(bem_w, bem_f1, trans);                                   // Step 8.
  y2 = (g1 - _bStar * bem_w) / delta;                               // Step 9. Scalar + ddot
  bem_w -= y2 * bem_v;                                              // Step 10. daxpy
  x = bem_w;
  y = y1 + y2;                                                      // Step 11. Scalar

}

// interchange B and BStar when DTR is true
template<class FACT> template<class T, bool DTR>
inline void KNBlockMatrix<FACT> :: __BEMWS
(
  T&           A,  JagVector2D& B,      JagVector2D&  BStar, KNMatrix&       D,
  JagVector2D& X,  KNVector&      Y,      JagVector2D&  F,     KNVector&       G,
  size_t j, const JagVector2D& V, const JagVector2D& VStar, const KNVector& delta, const KNVector& deltaStar, bool trans
)
{

  for (size_t k = j; k > 0; k--)
  {
    const size_t len   = A.col() + k - 1;
    const size_t aalen = A.col();
    const size_t ddlen = k - 1;
    const size_t idx   = k - 1;

    // Y_k = ( -v*_k 1 )^T . f_k
    Y(k) = F(k)(len);
    for (size_t r = 0; r < len; r++) Y(k) -= VStar(k)(r) * F(k)(r);  // ddot with range
    Y(k) /= deltaStar(k);

    // residuals
    for (size_t r = 0; r < aalen; r++)          F(k - 1)(r)         = F(k)(r)         - B(idx)(r) * Y(k);  // ???
    if (!DTR) for (size_t r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = F(k)(aalen + r) - D(r, idx) * Y(k);
    else     for (size_t r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = F(k)(aalen + r) - D(idx, r) * Y(k);
    G(k) = F(k)(len) - D(idx, idx) * Y(k);
  }

  A.solve(X(0), F(0), trans);
  double ypp;

  for (size_t k = 1; k <= j; k++)
  {
    const size_t len   = A.col() + k - 1;
    const size_t aalen = A.col();
    const size_t ddlen = k - 1;
    const size_t idx   = k - 1;

    ypp = G(k);
    for (size_t r = 0; r < aalen; r++)          ypp -= BStar(idx)(r) * X(k - 1)(r);     // ddot with range
    if (!DTR) for (size_t r = 0; r < ddlen; r++) ypp -= D(idx, r) * X(k - 1)(aalen + r);  // ddot with range this won't be vectorized.
    else     for (size_t r = 0; r < ddlen; r++) ypp -= D(r, idx) * X(k - 1)(aalen + r);
    ypp /= delta(k);

    for (size_t r = 0; r < len; r++) X(k)(r) = X(k - 1)(r) - V(k)(r) * ypp;    // daxpy with range
    X(k)(len) = Y(k) + ypp;
  }
}

// interchange B and BStar, when DTR is true
template<class FACT> template<class T, bool DTR>
inline void KNBlockMatrix<FACT> :: __BEMWF
(
  size_t bord,
  T&           A,  JagVector2D& B,      JagVector2D&  BStar, KNMatrix& D,
  JagVector2D& X,  KNVector&      Y,      JagVector2D&  F,     KNVector& G,
  JagVector2D& V,  JagVector2D& VStar,  KNVector&       delta, KNVector& deltaStar
)
{
  for (size_t k = 1; k < bord + 1; k++)
  {
    const size_t len   = A.col() + k - 1;
    const size_t aalen = A.col();
    const size_t ddlen = k - 1;
    const size_t idx   = k - 1;

    // with the ordinary solver
    for (size_t r = 0; r < aalen; r++)          F(k - 1)(r)         = B(idx)(r);   // dcopy with range
    if (!DTR) for (size_t r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(r, idx);  // dcopy with range
    else     for (size_t r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(idx, r);

    __BEMWS<T, DTR>(A, B, BStar, D, X, Y, F, G, k - 1, V, VStar, delta, deltaStar, DTR);

    // VV
    for (size_t r = 0; r < len; r++) V(k)(r) = X(k - 1)(r);  // dcopy with range
    // delta
    delta(k) = D(idx, idx);
    for (size_t r = 0; r < aalen; r++)          delta(k) -= BStar(idx)(r) * V(k)(r);   // ddot with range
    if (!DTR) for (size_t r = 0; r < ddlen; r++) delta(k) -= D(idx, r) * V(k)(aalen + r);
    else     for (size_t r = 0; r < ddlen; r++) delta(k) -= D(r, idx) * V(k)(aalen + r);

    // with the adjoint solver
    for (size_t r = 0; r < aalen; r++)          F(k - 1)(r)         = BStar(idx)(r);   // dcopy with range
    if (!DTR) for (size_t r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(idx, r);
    else     for (size_t r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(r, idx);

    __BEMWS<T, DTR>(A, B, BStar, D, X, Y, F, G, k - 1, VStar, V, deltaStar, delta, !DTR);

    // VV
    for (size_t r = 0; r < len; r++) VStar(k)(r) = X(k - 1)(r);             // dcopy with range
    // delta
    deltaStar(k) = D(idx, idx);
    for (size_t r = 0; r < aalen; r++)          deltaStar(k) -= B(idx)(r) * VStar(k)(r);          // ddot with range
    if (!DTR) for (size_t r = 0; r < ddlen; r++) deltaStar(k) -= D(r, idx) * VStar(k)(aalen + r);
    else     for (size_t r = 0; r < ddlen; r++) deltaStar(k) -= D(idx, r) * VStar(k)(aalen + r);
  }
}

// BEMW
// for trans == true interchange the borders
template<class FACT> template<class T, bool trans>
inline void KNBlockMatrix<FACT> :: __BEMW
(
  size_t bord,
  T&           A,  JagVector2D& B,      JagVector2D&  BStar, KNMatrix&       D,
  JagVector2D& X,  KNVector&      Y,      JagVector2D&  F,     KNVector&       G,
  JagVector2D& V,  JagVector2D& VStar,  KNVector&       delta, KNVector&       deltaStar,
  KNVector&      x,  KNVector&      y,      const KNVector& f,     const KNVector& g
)
{
  __BEMWF<T, trans>(bord, A, B, BStar, D, X, Y, F, G, V, VStar, delta, deltaStar);

  for (size_t i = 0; i < A.row(); i++) F(bord)(i)           = f(i);      // dcopy with range : F(bord)[rng(0,A.Row)] = f;
  for (size_t i = 0; i < bord; i++)    F(bord)(A.row() + i) = g(i);      // dcopy with range : F(bord)[rng(A.Row,A.Row+bord)] = g[rng(0,bord)];

  __BEMWS<T, trans>(A, B, BStar, D, X, Y, F, G, bord, V, VStar, delta, deltaStar, trans);

  for (size_t i = 0; i < A.row(); i++) x(i) = X(bord)(i);                // dcopy with range : x = X(bord)[rng(0,A.Row)];
  for (size_t i = 0; i < bord; i++)    y(i) = X(bord)(A.row() + i);      // dcopy with range : y = X(bord)[rng(A.Row,A.Row+bord)];
}

template<class FACT> template<bool trans>
void KNBlockMatrix<FACT>::multiply(KNVector& R1, double& R3, const KNVector& X1, const double& X3)
{
  const FACT&         _A11 = *A11;
  const JagVector2D&  _A13 = *A13;
  const JagVector2D&  _A31 = *A31;
  const KNMatrix&       _A33 = *A33;

  if (A11)
  {
    if (!trans) R1 = _A11 * X1;
    else         R1 = !_A11 * X1;
    if (A33)
    {
      if (!trans) R1 += X3 * _A13(0);
      else         R1 += X3 * _A31(0);
    }
  }
  if (A33)
  {
    if (!trans) R3 = _A31(0) * X1 + _A33(0, 0) * X3;
    else         R3 = _A13(0) * X1 + _A33(0, 0) * X3;
  }
}

template<class FACT> template<bool trans>
void KNBlockMatrix<FACT>::check(const KNVector& x, const double& z, const KNVector& f, const double& h)
{
  const FACT&         _A11 = *A11;
  const JagVector2D&  _A13 = *A13;
  const JagVector2D&  _A31 = *A31;
  const KNMatrix&       _A33 = *A33;
  KNVector              R1(x.size());
  double              R3;
  multiply<trans>(R1, R3, x, z);
  if (A11)
  {
    R1 -= f;
    std::cout << "F1: " << sqrt(f*f) << "R1: " << sqrt(R1*R1) << "\n";
  }
  if (A33)
  {
    R3 -= h;
    std::cout << "F3: " << sqrt(h*h) << "R3: " << sqrt(R3*R3) << "\n";
  }
}

template<class FACT> template<bool trans>
void KNBlockMatrix<FACT>::multiply(size_t bord, KNVector& R1, KNVector& R3, const KNVector& X1, const KNVector& X3)
{
  const FACT&         _A11 = *A11;
  const JagVector2D&  _A13 = *A13;
  const JagVector2D&  _A31 = *A31;
  const KNMatrix&       _A33 = *A33;

  if (A11)
  {
    if (!trans) R1 = _A11 * X1;
    else         R1 = !_A11 * X1;
    if (A33)
    {
      for (size_t i = 0; i < X1.size(); i++)
      {
        for (size_t j = 0; j < bord; j++)
        {
          if (!trans) R1(i) += _A13(j)(i) * X3(j);
          else         R1(i) += _A31(j)(i) * X3(j);
        }
      }
    }
  }
  else
  {
    R1.clear();
  }
  for (size_t i = 0; i < bord; i++)
  {
    if (A11)
    {
      if (!trans) R3(i) = _A31(i) * X1;
      else         R3(i) = _A13(i) * X1;
    }
    else
    {
      R3(i) = 0.0;
    }
  }
  if (A33)
  {
    for (size_t i = 0; i < bord; i++)
    {
      for (size_t j = 0; j < bord; j++)
      {
        if (!trans) R3(i) += _A33(i, j) * X3(j);
        else         R3(i) += _A33(j, i) * X3(j);
      }
    }
  }
}

template<class FACT> template<bool trans>
void KNBlockMatrix<FACT>::multiply(KNBlockVector& X, const KNBlockVector& F, size_t bord)
{
  multiply<trans> (bord, X.getV1(), X.getV3(), F.getV1(), F.getV3());
}

template<class FACT> template<bool trans>
void KNBlockMatrix<FACT>::check(size_t bord, const KNVector& X1, const KNVector& X3, const KNVector& F1, const KNVector& F3)
{
  // this may be buggy with some compilers
  const FACT&         _A11 = *A11;
  const KNMatrix&       _A33 = *A33;
  //multiply back...
  KNVector R1(X1.size());
  KNVector R3(X3.size());
  multiply<trans>(bord, R1, R3, X1, X3);
  if (A11)
  {
    R1 -= F1;
    std::cout << "F1: " << sqrt(F1*F1) << " R1:\t" << sqrt(R1*R1) << "\n";
  }
  if (A33)
  {
    R3 -= F3;
    std::cout << "F3: " << sqrt(F3*F3) << " R3(" << bord << ")\t" << sqrt(R3*R3) << "\n";
  }
}

template<class FACT>
void KNBlockMatrix<FACT>::check(const KNBlockVector& X, const KNBlockVector& F, size_t bord)
{
  check<false>(bord, X.getV1(), X.getV3(), F.getV1(), F.getV3());
}

// Wrapper functions
template<class FACT>
void KNBlockMatrix<FACT>::solve(KNVector& x, double& z, const KNVector& f, const double& h)   // BEM
{
  __BEM<FACT, false>(*A11, (*A13)(0), (*A31)(0), (*A33)(0, 0), x, z, f, h, bem_v_1, bem_vStar_1, bem_w_1, bem_f1_1);
}

template<class FACT>
void KNBlockMatrix<FACT>::solveTr(KNVector& x, double& z, const KNVector& f, const double& h)   // BEM
{
  __BEM<FACT, true>(*A11, (*A31)(0), (*A13)(0), (*A33)(0, 0), x, z, f, h, bem_v_1, bem_vStar_1, bem_w_1, bem_f1_1);
}

template<class FACT>
void KNBlockMatrix<FACT>::solve(size_t bord, KNVector& X1, KNVector& X3, const KNVector& F1, const KNVector& F3)   // BEMW
{
  __BEMW<FACT, false>(bord, *A11, *A13, *A31, *A33, *XX1, *YY1, *FF1, *GG1, *VV1, *VV1Star, *delta1, *delta1Star, X1, X3, F1, F3);
  // Check<false>( bord, X1, X3, F1, F3 );
}

template<class FACT>
void KNBlockMatrix<FACT>::solveTr(size_t bord, KNVector& X1, KNVector& X3, const KNVector& F1, const KNVector& F3)   // BEMW
{
  __BEMW<FACT, true>(bord, *A11, *A31, *A13, *A33, *XX1, *YY1, *FF1, *GG1, *VV1, *VV1Star, *delta1, *delta1Star, X1, X3, F1, F3);
  // Check<true>( bord, X1, X3, F1, F3 );
}

template<class FACT>
void KNBlockMatrix<FACT>::solve(KNBlockVector& X, const KNBlockVector& F)
{
  P_ASSERT_X(A11 != 0, "KNBlockMatrix::Solve Error");
  if (A33 != 0)
  {
    // MBEW v BEM
    if (A33->col() == 1)
    {
      solve(X.getV1(), X.getV3()(0), F.getV1(), F.getV3()(0));
    }
    else
    {
      solve(A33->col(), X.getV1(), X.getV3(), F.getV1(), F.getV3());
    }
  }
  else
  {
    A11->solve(X.getV1(), F.getV1());
  }
  // check( X, F, F.getV3().size() );
}

template<class FACT>
void KNBlockMatrix<FACT>::solve(KNBlockVector& X, const KNBlockVector& F, size_t bord)
{
  P_ASSERT_X(A11 != 0, "KNBlockMatrix::Solve Error");
  if (A33 != 0)
  {
    // MBEW v BEM
    if (bord == 1)
    {
      solve(X.getV1(), X.getV3()(0), F.getV1(), F.getV3()(0));
    }
    else
    {
      solve(bord, X.getV1(), X.getV3(), F.getV1(), F.getV3());
    }
  }
  else
  {
    A11->solve(X.getV1(), F.getV1());
  }
  // check( X, F, F.getV3().size() );
}

template<class FACT>
void KNBlockMatrix<FACT>::solveDirect(KNBlockVector& X, const KNBlockVector& F)
{
  size_t dim1, dim3;
  if (A11) dim1 = A11->col();
  else dim1 = 0;
  if (A33) dim3 = A33->col();
  else dim3 = 0;

  KNLuSparseMatrix* A11S = A11;

  KNLuSparseMatrix  AA('R', dim1 + dim3, A11S->nonzeros() + (dim1 + dim3)*(dim3) + dim1*dim3 + 10);
  KNVector  XX(dim1 + dim3), FF(dim1 + dim3);

  for (size_t i = 0; i < dim1; i++)
  {
    AA.newLine(A11S->lineLength(i) + dim3);
    for (size_t j = 0; j < A11S->lineLength(i); j++)
    {
      AA.writeIndex(i, j) = A11S->writeIndex(i, j);
      AA.writeData(i, j) = A11S->writeData(i, j);
    }
    for (size_t j = 0; j < dim3; j++)
    {
      AA.writeIndex(i, A11S->lineLength(i) + j) = static_cast<int>(dim1 + j);
      AA.writeData(i, A11S->lineLength(i) + j) = getA13(j)(i);
    }
    FF(i) = F.getV1()(i);
  }
  for (size_t i = 0; i < dim3; i++)
  {
    AA.newLine(dim1 + dim3);
    for (size_t j = 0;j < dim1;j++)
    {
      AA.writeIndex(dim1 + i, j) = static_cast<int>(j);
      AA.writeData(dim1 + i, j) = getA31(i)(j);
    }
    for (size_t j = 0;j < dim3;j++)
    {
      AA.writeIndex(dim1 + i, dim1 + j) = static_cast<int>(dim1 + j);
      AA.writeData(dim1 + i, dim1 + j) = getA33(i, j);
    }
    FF(dim1 + i) = F.getV3()(i);
  }

  AA.solve(XX, FF);

  for (size_t i = 0; i < dim1; i++)  X.getV1()(i) = XX(i);
  for (size_t i = 0; i < dim3; i++)  X.getV3()(i) = XX(dim1 + i);

//  check( X, F );
  std::cout << "DR ";

}

#endif
