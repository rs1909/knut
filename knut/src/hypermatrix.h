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

class HyperVector
{
  public:

    inline HyperVector(int i, int, int k) : V1(i), V3(k)
    { }
    inline HyperVector(const HyperVector& hv) : V1(hv.V1), V3(hv.V3)
    { }
    inline ~HyperVector()
    { }

    inline Vector& getV1()
    {
      return V1;
    }
    inline const Vector& getV1() const
    {
      return V1;
    }
    inline Vector& getV3()
    {
      return V3;
    }
    inline const Vector& getV3() const
    {
      return V3;
    }

  private:

    Vector V1;
    Vector V3;
};

// #define FACT SpFact
template< class FACT > class HyMatrix
{
  public:

    // constructing Sparse-Dense type hypermatrix
    HyMatrix(int i, int, int k, int nz);
    ~HyMatrix();

    inline FACT&        getA11()
    {
      return *A11;
    }
    inline JagVector2D& getA13()
    {
      return *A13;
    }
    inline Vector&      getA13(int i)
    {
      return (*A13)(i);
    }
    inline JagVector2D& getA31()
    {
      return *A31;
    }
    inline Vector&      getA31(int i)
    {
      return (*A31)(i);
    }
    inline Matrix&      getA33()
    {
      return *A33;
    }
    inline double&      getA33(int i, int j)
    {
      return (*A33)(i, j);
    }

    template< class T, bool trans > inline void __BEM
    (
      T& _A, Vector& _b, Vector& _bStar, double& _d,
      Vector& x, double& y, const Vector& f, const double& g,
      Vector& bem_v, Vector& bem_vStar, Vector& bem_w, Vector& bem_f1
    );
    template< class T, bool DTR > inline void __BEMWS
    (
      T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
      JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
      int j, const JagVector2D& V, const JagVector2D& VStar, const Vector& delta, const Vector& deltaStar, bool trans
    );
    template< class T, bool DTR > inline void __BEMWF
    (
      int bord,
      T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix& D,
      JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector& G,
      JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector& deltaStar
    );
    template< class T, bool trans > inline void __BEMW
    (
      int bord,
      T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
      JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
      JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector&       deltaStar,
      Vector&      x,  Vector&      y,      const Vector& f,     const Vector& g
    );

    void Solve(HyperVector& X, const HyperVector& F);

    void Solve(HyperVector& X, const HyperVector& F, int bord);

    template<bool trans> void Multiply(HyperVector& X, const HyperVector& F, int bord);

    void Check(const HyperVector& X, const HyperVector& F, int bord);

    void SolveDIRECT(HyperVector& X, const HyperVector& F);

    void Solve(Vector& x, const Vector& f)
    {
      A11->Solve(x, f);
    }

    template<bool trans> void Multiply(Vector& R1, double& R3, const Vector& X1, const double& X3);

    template<bool trans> void Check(const Vector& x, const double& z, const Vector& f, const double& h);

    void Solve(Vector& x, double& z, const Vector& f, const double& h);   // BEM

    void SolveTR(Vector& x, double& z, const Vector& f, const double& h);   // BEM

    template<bool trans> void Multiply(int bord, Vector& R1, Vector& R3, const Vector& X1, const Vector& X3);

    template<bool trans> void Check(int bord, const Vector& x, const Vector& z, const Vector& f, const Vector& h);

    void Solve(int bord, Vector& x, Vector& z, const Vector& f, const Vector& h);   // BEMW

    void SolveTR(int bord, Vector& x, Vector& z, const Vector& f, const Vector& h);   // BEMW

  protected:

    // matrix blocks

    FACT*    A11;
    JagVector2D* A13;
    JagVector2D* A31;
    Matrix*   A33;

  private:
    // temporary storage for factorization

    // BEMW 1.
    JagVector2D* VV1;
    JagVector2D* VV1Star;
    JagVector2D* FF1;
    Vector*   GG1;
    JagVector2D* XX1;
    Vector*   YY1;
    Vector*   delta1;
    Vector*   delta1Star;

    // BEM temporaries
    Vector         bem_v_1;
    Vector         bem_vStar_1;
    Vector         bem_w_1;
    Vector         bem_f1_1;
    Vector         bem_v_2;
    Vector         bem_vStar_2;
    Vector         bem_w_2;
    Vector         bem_f1_2;
};

typedef HyMatrix<SpFact> HyperMatrix;

///--------------------------------
///
/// Now the implementation
///
///--------------------------------

// constructing Sparse-Dense type hypermatrix
template<class FACT> HyMatrix<FACT> :: HyMatrix(int i, int j, int k, int nz) :
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
  if (k != 0) A33 = new Matrix(k, k);
  else A33 = 0;

  if ((i != 0) && (k != 0))
  {
    A13 = new JagVector2D(k, i);
    A31 = new JagVector2D(k, i);

    VV1 = new JagVector2D(k + 1);
    VV1Star = new JagVector2D(k + 1);
    FF1 = new JagVector2D(k + 1);
    GG1 = new Vector(k + 1);
    XX1 = new JagVector2D(k + 1);
    YY1 = new Vector(k + 1);
    delta1 = new Vector(k + 1);
    delta1Star = new Vector(k + 1);

    for (int r = 0; r < k + 1; r++)
    {
      (*VV1)(r).Init(i + r);
      (*VV1Star)(r).Init(i + r);
      (*XX1)(r).Init(i + r);
      (*FF1)(r).Init(i + r);
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

template<class FACT> HyMatrix<FACT> :: ~HyMatrix()
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
//  Vector v( _b.size() ); == A.size()                     bem_v
//  Vector vStar( _b.size() );                             bem_vStar
//  Vector w( _b.size() );                                 bem_w
//  Vector f1( _b.size() );                                bem_f1;

//  GMBE temporaries
//  Vector xi( _A32.size() ); == dim2                      gmbe_xi;
//  Vector cc( _A31 ); == dim1                             gmbe_cc
//  Vector fr( F2 ); == dim2                               cm_fr

// GMBEW temporaries
//  JagVector2D xi( bord, _A32(0).size() );  -> different  gmbew_xi
//  JagVector2D cc( _A31 );                  -> different  gmbew_cc
//  Matrix dd( bord, bord );                               gmbew_dd
//  Vector gg( bord );                                     gmbew_gg
//  Vector fr( F2.size() ); -> same as GMBE                cm_fr
//  Vector gr( bord );                                     gmbew_gr
//  Vector yy( bord );                                     gmbew_yy
//  Matrix m_beta( bord, bord );                           gmbew_m_beta

// BEM
// swap b and bStar if trans==true
template<class FACT> template<class T, bool trans>
inline void HyMatrix<FACT> :: __BEM(T& _A, Vector& _b, Vector& _bStar, double& _d,
                                    Vector& x, double& y, const Vector& f, const double& g,
                                    Vector& bem_v, Vector& bem_vStar, Vector& bem_w, Vector& bem_f1)
{

  double delta, deltaStar;
  double g1;
  double y1, y2;

  // Doolittle
  _A.Solve(bem_vStar, _bStar, !trans);                              // Step 1.
  deltaStar = _d - bem_vStar * _b;                                  // Step 2. Scalar + ddot
  // Crout
  _A.Solve(bem_v, _b, trans);                                       // Step 3.
  delta = _d - _bStar * bem_v;                                      // Step 4. Scalar + ddot
  // approx Y
  y1 = (g - bem_vStar * f) / deltaStar;                             // Step 5. Scalar + ddot
  // residuals
  bem_f1 = f;
  bem_f1 -= y1 * _b;                                                // Step 6. daxpy
  g1 = g - _d * y1;                                                 // Step 7. Scalar

  // residual corrections
  _A.Solve(bem_w, bem_f1, trans);                                   // Step 8.
  y2 = (g1 - _bStar * bem_w) / delta;                               // Step 9. Scalar + ddot
  bem_w -= y2 * bem_v;                                              // Step 10. daxpy
  x = bem_w;
  y = y1 + y2;                                                      // Step 11. Scalar

}

// interchange B and BStar when DTR is true
template<class FACT> template<class T, bool DTR>
inline void HyMatrix<FACT> :: __BEMWS
(
  T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
  JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
  int j, const JagVector2D& V, const JagVector2D& VStar, const Vector& delta, const Vector& deltaStar, bool trans
)
{

  for (int k = j; k > 0; k--)
  {
    const int len   = A.Col() + k - 1;
    const int aalen = A.Col();
    const int ddlen = k - 1;
    const int idx   = k - 1;

    // Y_k = ( -v*_k 1 )^T . f_k
    Y(k) = F(k)(len);
    for (int r = 0; r < len; r++) Y(k) -= VStar(k)(r) * F(k)(r);  // ddot with range
    Y(k) /= deltaStar(k);

    // residuals
    for (int r = 0; r < aalen; r++)          F(k - 1)(r)         = F(k)(r)         - B(idx)(r) * Y(k);  // ???
    if (!DTR) for (int r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = F(k)(aalen + r) - D(r, idx) * Y(k);
    else     for (int r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = F(k)(aalen + r) - D(idx, r) * Y(k);
    G(k) = F(k)(len) - D(idx, idx) * Y(k);
  }

  A.Solve(X(0), F(0), trans);
  double ypp;

  for (int k = 1; k <= j; k++)
  {
    const int len   = A.Col() + k - 1;
    const int aalen = A.Col();
    const int ddlen = k - 1;
    const int idx   = k - 1;

    ypp = G(k);
    for (int r = 0; r < aalen; r++)          ypp -= BStar(idx)(r) * X(k - 1)(r);     // ddot with range
    if (!DTR) for (int r = 0; r < ddlen; r++) ypp -= D(idx, r) * X(k - 1)(aalen + r);  // ddot with range this won't be vectorized.
    else     for (int r = 0; r < ddlen; r++) ypp -= D(r, idx) * X(k - 1)(aalen + r);
    ypp /= delta(k);

    for (int r = 0; r < len; r++) X(k)(r) = X(k - 1)(r) - V(k)(r) * ypp;    // daxpy with range
    X(k)(len) = Y(k) + ypp;
  }
}

// interchange B and BStar, when DTR is true
template<class FACT> template<class T, bool DTR>
inline void HyMatrix<FACT> :: __BEMWF
(
  int bord,
  T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix& D,
  JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector& G,
  JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector& deltaStar
)
{
  for (int k = 1; k < bord + 1; k++)
  {
    const int len   = A.Col() + k - 1;
    const int aalen = A.Col();
    const int ddlen = k - 1;
    const int idx   = k - 1;

    // with the ordinary solver
    for (int r = 0; r < aalen; r++)          F(k - 1)(r)         = B(idx)(r);   // dcopy with range
    if (!DTR) for (int r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(r, idx);  // dcopy with range
    else     for (int r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(idx, r);

    __BEMWS<T, DTR>(A, B, BStar, D, X, Y, F, G, k - 1, V, VStar, delta, deltaStar, DTR);

    // VV
    for (int r = 0; r < len; r++) V(k)(r) = X(k - 1)(r);  // dcopy with range
    // delta
    delta(k) = D(idx, idx);
    for (int r = 0; r < aalen; r++)          delta(k) -= BStar(idx)(r) * V(k)(r);   // ddot with range
    if (!DTR) for (int r = 0; r < ddlen; r++) delta(k) -= D(idx, r) * V(k)(aalen + r);
    else     for (int r = 0; r < ddlen; r++) delta(k) -= D(r, idx) * V(k)(aalen + r);

    // with the adjoint solver
    for (int r = 0; r < aalen; r++)          F(k - 1)(r)         = BStar(idx)(r);   // dcopy with range
    if (!DTR) for (int r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(idx, r);
    else     for (int r = 0; r < ddlen; r++) F(k - 1)(aalen + r) = D(r, idx);

    __BEMWS<T, DTR>(A, B, BStar, D, X, Y, F, G, k - 1, VStar, V, deltaStar, delta, !DTR);

    // VV
    for (int r = 0; r < len; r++) VStar(k)(r) = X(k - 1)(r);             // dcopy with range
    // delta
    deltaStar(k) = D(idx, idx);
    for (int r = 0; r < aalen; r++)          deltaStar(k) -= B(idx)(r) * VStar(k)(r);          // ddot with range
    if (!DTR) for (int r = 0; r < ddlen; r++) deltaStar(k) -= D(r, idx) * VStar(k)(aalen + r);
    else     for (int r = 0; r < ddlen; r++) deltaStar(k) -= D(idx, r) * VStar(k)(aalen + r);
  }
}

// BEMW
// for trans == true interchange the borders
template<class FACT> template<class T, bool trans>
inline void HyMatrix<FACT> :: __BEMW
(
  int bord,
  T&           A,  JagVector2D& B,      JagVector2D&  BStar, Matrix&       D,
  JagVector2D& X,  Vector&      Y,      JagVector2D&  F,     Vector&       G,
  JagVector2D& V,  JagVector2D& VStar,  Vector&       delta, Vector&       deltaStar,
  Vector&      x,  Vector&      y,      const Vector& f,     const Vector& g
)
{
  __BEMWF<T, trans>(bord, A, B, BStar, D, X, Y, F, G, V, VStar, delta, deltaStar);

  for (int i = 0; i < A.Row(); i++) F(bord)(i)           = f(i);      // dcopy with range : F(bord)[rng(0,A.Row)] = f;
  for (int i = 0; i < bord; i++)    F(bord)(A.Row() + i) = g(i);      // dcopy with range : F(bord)[rng(A.Row,A.Row+bord)] = g[rng(0,bord)];

  __BEMWS<T, trans>(A, B, BStar, D, X, Y, F, G, bord, V, VStar, delta, deltaStar, trans);

  for (int i = 0; i < A.Row(); i++) x(i) = X(bord)(i);                // dcopy with range : x = X(bord)[rng(0,A.Row)];
  for (int i = 0; i < bord; i++)    y(i) = X(bord)(A.Row() + i);      // dcopy with range : y = X(bord)[rng(A.Row,A.Row+bord)];
}

template<class FACT> template<bool trans>
void HyMatrix<FACT>::Multiply(Vector& R1, double& R3, const Vector& X1, const double& X3)
{
  const FACT&         _A11 = *A11;
  const JagVector2D&  _A13 = *A13;
  const JagVector2D&  _A31 = *A31;
  const Matrix&       _A33 = *A33;

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
void HyMatrix<FACT>::Check(const Vector& x, const double& z, const Vector& f, const double& h)
{
  const FACT&         _A11 = *A11;
  const JagVector2D&  _A13 = *A13;
  const JagVector2D&  _A31 = *A31;
  const Matrix&       _A33 = *A33;
  Vector              R1(x.size());
  double              R3;
  Multiply<trans>(R1, R3, x, z);
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
void HyMatrix<FACT>::Multiply(int bord, Vector& R1, Vector& R3, const Vector& X1, const Vector& X3)
{
  const FACT&         _A11 = *A11;
  const JagVector2D&  _A13 = *A13;
  const JagVector2D&  _A31 = *A31;
  const Matrix&       _A33 = *A33;

  if (A11)
  {
    if (!trans) R1 = _A11 * X1;
    else         R1 = !_A11 * X1;
    if (A33)
    {
      for (int i = 0; i < X1.size(); i++)
      {
        for (int j = 0; j < bord; j++)
        {
          if (!trans) R1(i) += _A13(j)(i) * X3(j);
          else         R1(i) += _A31(j)(i) * X3(j);
        }
      }
    }
  }
  else
  {
    R1.Clear();
  }
  for (int i = 0; i < bord; i++)
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
    for (int i = 0; i < bord; i++)
    {
      for (int j = 0; j < bord; j++)
      {
        if (!trans) R3(i) += _A33(i, j) * X3(j);
        else         R3(i) += _A33(j, i) * X3(j);
      }
    }
  }
}

template<class FACT> template<bool trans>
void HyMatrix<FACT>::Multiply(HyperVector& X, const HyperVector& F, int bord)
{
  Multiply<trans> (bord, X.getV1(), X.getV3(), F.getV1(), F.getV3());
}

template<class FACT> template<bool trans>
void HyMatrix<FACT>::Check(int bord, const Vector& X1, const Vector& X3, const Vector& F1, const Vector& F3)
{
  // this may be buggy with some compilers
  const FACT&         _A11 = *A11;
  const Matrix&       _A33 = *A33;
  //multiply back...
  Vector R1(X1.size());
  Vector R3(X3.size());
  Multiply<trans>(bord, R1, R3, X1, X3);
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
void HyMatrix<FACT>::Check(const HyperVector& X, const HyperVector& F, int bord)
{
  Check<false>(bord, X.getV1(), X.getV3(), F.getV1(), F.getV3());
}

// Wrapper functions
template<class FACT>
void HyMatrix<FACT>::Solve(Vector& x, double& z, const Vector& f, const double& h)   // BEM
{
  __BEM<FACT, false>(*A11, (*A13)(0), (*A31)(0), (*A33)(0, 0), x, z, f, h, bem_v_1, bem_vStar_1, bem_w_1, bem_f1_1);
}

template<class FACT>
void HyMatrix<FACT>::SolveTR(Vector& x, double& z, const Vector& f, const double& h)   // BEM
{
  __BEM<FACT, true>(*A11, (*A31)(0), (*A13)(0), (*A33)(0, 0), x, z, f, h, bem_v_1, bem_vStar_1, bem_w_1, bem_f1_1);
}

template<class FACT>
void HyMatrix<FACT>::Solve(int bord, Vector& X1, Vector& X3, const Vector& F1, const Vector& F3)   // BEMW
{
  __BEMW<FACT, false>(bord, *A11, *A13, *A31, *A33, *XX1, *YY1, *FF1, *GG1, *VV1, *VV1Star, *delta1, *delta1Star, X1, X3, F1, F3);
  // Check<false>( bord, X1, X3, F1, F3 );
}

template<class FACT>
void HyMatrix<FACT>::SolveTR(int bord, Vector& X1, Vector& X3, const Vector& F1, const Vector& F3)   // BEMW
{
  __BEMW<FACT, true>(bord, *A11, *A31, *A13, *A33, *XX1, *YY1, *FF1, *GG1, *VV1, *VV1Star, *delta1, *delta1Star, X1, X3, F1, F3);
  // Check<true>( bord, X1, X3, F1, F3 );
}

template<class FACT>
void HyMatrix<FACT>::Solve(HyperVector& X, const HyperVector& F)
{
  P_ASSERT_X(A11 != 0, "HyMatrix::Solve Error");
  if (A33 != 0)
  {
    // MBEW v BEM
    if (A33->Col() == 1)
    {
      Solve(X.getV1(), X.getV3()(0), F.getV1(), F.getV3()(0));
    }
    else
    {
      Solve(A33->Col(), X.getV1(), X.getV3(), F.getV1(), F.getV3());
    }
  }
  else
  {
    A11->Solve(X.getV1(), F.getV1());
  }
  // Check( X, F, F.getV3().size() );
}

template<class FACT>
void HyMatrix<FACT>::Solve(HyperVector& X, const HyperVector& F, int bord)
{
  P_ASSERT_X(A11 != 0, "HyMatrix::Solve Error");
  if (A33 != 0)
  {
    // MBEW v BEM
    if (bord == 1)
    {
      Solve(X.getV1(), X.getV3()(0), F.getV1(), F.getV3()(0));
    }
    else
    {
      Solve(bord, X.getV1(), X.getV3(), F.getV1(), F.getV3());
    }
  }
  else
  {
    A11->Solve(X.getV1(), F.getV1());
  }
  // Check( X, F, F.getV3().size() );
}

template<class FACT>
void HyMatrix<FACT>::SolveDIRECT(HyperVector& X, const HyperVector& F)
{
  int dim1, dim3;
  if (A11) dim1 = A11->Col();
  else dim1 = 0;
  if (A33) dim3 = A33->Col();
  else dim3 = 0;

  SpFact* A11S = A11;

  SpFact  AA('R', dim1 + dim3, A11S->GetNZ() + (dim1 + dim3)*(dim3) + dim1*dim3 + 10);
  Vector  XX(dim1 + dim3), FF(dim1 + dim3);

  for (int i = 0; i < dim1; i++)
  {
    AA.NewL(A11S->GetL(i) + dim3);
    for (int j = 0; j < A11S->GetL(i); j++)
    {
      AA.WrLi(i, j) = A11S->WrLi(i, j);
      AA.WrLx(i, j) = A11S->WrLx(i, j);
    }
    for (int j = 0; j < dim3; j++)
    {
      AA.WrLi(i, A11S->GetL(i) + j) = dim1 + j;
      AA.WrLx(i, A11S->GetL(i) + j) = getA13(j)(i);
    }
    FF(i) = F.getV1()(i);
  }
  for (int i = 0; i < dim3; i++)
  {
    AA.NewL(dim1 + dim3);
    for (int j = 0;j < dim1;j++)
    {
      AA.WrLi(dim1 + i, j) = j;
      AA.WrLx(dim1 + i, j) = getA31(i)(j);
    }
    for (int j = 0;j < dim3;j++)
    {
      AA.WrLi(dim1 + i, dim1 + j) = dim1 + j;
      AA.WrLx(dim1 + i, dim1 + j) = getA33(i, j);
    }
    FF(dim1 + i) = F.getV3()(i);
  }

  AA.Solve(XX, FF);

  for (int i = 0; i < dim1; i++)  X.getV1()(i) = XX(i);
  for (int i = 0; i < dim3; i++)  X.getV3()(i) = XX(dim1 + i);

//  Check( X, F );
  std::cout << "DR ";

}

#endif
