// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef SPMATRIX_H
#define SPMATRIX_H

#include "matrix.h"
#include "plot.h"
#include "knerror.h"
#include <vector>

extern "C"
{

#include "cspblas.h"
#include "umfpack.h"

}

// **************************************************************************//
//                                                                           //
// **********************        SpMatrix Class        **********************//
//                                                                           //
// **************************************************************************//

class SpMatrix
{

  protected:

    char format;
    int n;
    int m;
    int size;
    int *Ap;
    int *Ai;
    double *Ax;

  public:

    inline SpMatrix() : format('R'), n(0), m(0), size(0), Ap(0), Ai(0), Ax(0)
    { }
    inline SpMatrix(char F, int n_, int m_, int nz)
    {
      init(F, n_, m_, nz);
    }
    inline SpMatrix(char F, int n_, int nz)
    {
      init(F, n_, n_, nz);
    }
    inline SpMatrix(const SpMatrix& M)
    {
      init(M);
    }
    inline virtual ~SpMatrix();

    inline void init(char F, int n_, int m_, int nz);
    inline void init(const SpMatrix& M);

    inline void mmx(enum cspblas_Trans trans, double* out, const double* in, double alpha) const
    {
      cspblas_mmx(format, trans, n, m, Ap, Ai, Ax, out, in, alpha);
    }
    inline void mmxpy(enum cspblas_Trans trans, double* out, const double* in, double alpha, const double* Y, double beta) const
    {
      cspblas_mmxpy(format, trans, n, m, Ap, Ai, Ax, out, in, alpha, Y, beta);
    }
    inline void mmxm(enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha, int nrhs) const
    {
      cspblas_mmxm(format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, nrhs);
    }
    inline void mmxmpym(enum cspblas_Trans trans, double* out, int ldout, const double* in, int ldin, double alpha,
                        const double* Y, int ldY, double beta, int nrhs) const
    {
      cspblas_mmxmpym(format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, Y, ldY, beta, nrhs);
    }

    inline __op_mul_vec<SpMatrix, Vector>     operator*(const Vector& v) const
    {
      return __op_mul_vec<SpMatrix, Vector>(__scal_vec_trans<SpMatrix>(*this), __scal_vec_trans<Vector>(v));
    }
    inline __op_mul_vec<SpMatrix, Matrix>     operator*(const Matrix& v) const
    {
      return __op_mul_vec<SpMatrix, Matrix>(__scal_vec_trans<SpMatrix>(*this), __scal_vec_trans<Matrix>(v));
    }
    inline __op_mul_vec_rng<SpMatrix, Vector> operator*(const __scal_vec_trans_rng<Vector> v) const
    {
      return __op_mul_vec_rng<SpMatrix, Vector>(__scal_vec_trans_rng<SpMatrix>(*this), v);
    }
    inline __op_mul_vec_rng<SpMatrix, Matrix> operator*(const __scal_vec_trans_rng<Matrix> v) const
    {
      return __op_mul_vec_rng<SpMatrix, Matrix>(__scal_vec_trans_rng<SpMatrix>(*this), v);
    }
    inline __scal_vec_trans<SpMatrix>        operator!() const
    {
      return __scal_vec_trans<SpMatrix>(*this, 1.0, Trans);
    }
    inline __scal_vec_trans_rng<SpMatrix>    operator[ ](const rng r) const
    {
      return __scal_vec_trans_rng<SpMatrix>(*this, r);
    }

    inline void AX(double* out, const double* in, double alpha, bool trans) const;
    inline void AXpY(double* out, const double* in, const double* Y, double alpha, double beta, bool trans) const;

    inline int row() const
    {
      if (format == 'R') return n;
      else return m;
    }
    inline int col() const
    {
      if (format == 'R') return m;
      else return n;
    }

    /// clear the matrix
    /// these have to be virtual, because these might be called as SpFact
    virtual void clear();
    virtual void clear(char F);
    /// checks the structure of the matrix
    void Check();

    // Fill in routines
    /// Creates a new line in the matrix
    inline int NewL(int size_);
    /// Writes or returns the row or column index
    /// into line `l' and the `e'-th element
    inline int& WrLi(int l, int e);
    /// Writes into line `l' and the `e'-th element
    inline double& WrLx(int l, int e);
    /// returns the length of the n_ -th line in the matrix
    inline int GetL(int n_);
    /// returns the nonzero elements in the matrix
    inline int GetNZ()
    {
      return Ap[n];
    }
    /// returns the number of lines, e.g. columns or rows in the matrix depending on format
    inline int GetN()
    {
      return n;
    }
    // Computation routines
    /// transposes the matrix into the other format
    void Swap();
    // these are used for debugging only
    /// plots the structure of the matrix
    void sparsityPlot(GnuPlot& pl);
    /// prints out Ap
    void PrintAp()
    {
      for (int i = 0; i < n + 1; i++) std::cout << Ap[i] << '\t';
      std::cout << '\n';
    }
    /// prints the whole matrix onto the screen
    void print();
};

// **************************************************************************//
//                                                                           //
// **********************         SpFact Class         **********************//
//                                                                           //
// **************************************************************************//

class SpFact : public SpMatrix
{
  private:

    bool   fact;
    void   *Numeric;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    // temporary space
    int    *Wi;
    double *W;

  public:
    void init(int nn_);
    SpFact(char F, int n_, int m_, int nz);
    SpFact(char F, int n_, int nz);
    SpFact(SpMatrix& M);
    inline SpFact(SpFact&) : SpMatrix('F', 1, 1)
    {
      std::cout << "SpFact::SpFact(SpFact&): not implemented\n";
    }
    ~SpFact();

    inline void modified()
    {
      fact = false;
    }

    inline void SetIter(int N)
    {
      Control[UMFPACK_IRSTEP] = N;
    }
    void clear();
    void clear(char F);

    void solve(double* x, double* b, bool trans = false);
    void solve(Vector& x, const Vector& b, bool trans = false);
    void solve(Matrix& x, const Matrix& b, bool trans = false);
    /// get the diagonal of the Upper factor
    int GetDX(Vector& V)
    {
      if (!fact) luFactorize();
      return umfpack_di_get_numeric(0, 0, 0, 0, 0, 0, 0, 0, V.v, 0, 0, Numeric);
    }
  private:
    void luFactorize();
};

// void NColloc::StabJac( StabMat& AB, const Vector& par, const JagMatrix3D& solData );
// this computes its own matrix structures, because it is nowhere else needed: kkSI, eeSI, rrSI, ddSI, etc.
// However the variables will be contained within the NColloc class

class StabMatrix
{
  public:

    StabMatrix(int nmat_, int n_, int nz) : A0('R', n_,  nz), AI(nmat_)
    {
      for (int i = 0; i < AI.size(); i++) AI(i).init('R', n_, n_, nz);
      RESID = new double[ AI.size() * n_ + 1 ];
      isINIT = false;
    }
    ~StabMatrix()
    {
      delete[] RESID;
    }

    int                nmat()
    {
      return AI.size();
    }
    SpFact&            getA0()
    {
      return A0;
    }
    Array1D<SpMatrix>& getAI()
    {
      return AI;
    }
    SpMatrix&          getAI(int i)
    {
      return AI(i);
    }

    void eigenvalues(Vector& wr, Vector& wi);

  private:

    SpFact             A0;
    Array1D<SpMatrix>  AI;

    double* RESID;
    bool    isINIT;
};

// Implementation of SpMatrix

inline void SpMatrix::init(char F, int n_, int m_, int nz)
{
  if ((F != 'R') && (F != 'C')) std::cout << "SpMatrix::CONSTRUCTOR: invalid format specification.\n";
  format = F;
  n = 0;
  m = m_;
  size = nz;
  Ap = new int[n_+1];
  Ai = new int[nz];
  Ax = new double[nz];
  for (int i = 0; i < nz; i++)
  {
    Ai[i] = 0;
    Ax[i] = 0.0;
  }
  for (int i = 0; i < n_ + 1; i++)
  {
    Ap[i] = 0;
  }
}

inline void SpMatrix::init(const SpMatrix& M)
{
  format = M.format;
  n = M.n;
  m = M.m;
  size = M.size;

  Ap = new int[M.n+1];
  Ai = new int[size];
  Ax = new double[size];
  for (int i = 0; i < size; i++)
  {
    Ai[i] = M.Ai[i];
    Ax[i] = M.Ax[i];
  }
  for (int i = 0; i < n + 1; i++)
  {
    Ap[i] = M.Ap[i];
  }
}

inline SpMatrix::~SpMatrix()
{
  delete []Ap;
  delete []Ai;
  delete []Ax;
}

inline void SpMatrix::clear()
{

  for (int i = 0; i < n + 1; i++)
  {
    Ap[i] = 0;
  }
  n = 0;
  for (int i = 0; i < size; i++)
  {
    Ai[i] = 0;
    Ax[i] = 0.0;
  }
}

inline void SpMatrix::clear(char F)
{

  for (int i = 0; i < n + 1; i++)
  {
    Ap[i] = 0;
  }
  n = 0;
  format = F;
  for (int i = 0; i < size; i++)
  {
    Ai[i] = 0;
    Ax[i] = 0.0;
  }
}

inline int SpMatrix::NewL(int size_)
{
  n++;
  Ap[n] = Ap[n-1] + size_;
  return n -1;
}

inline int& SpMatrix::WrLi(int l, int e)
{
#ifdef DEBUG
//   std::cout << "n:" << n << " l:" << l << " lsz:" << Ap[l+1] - Ap[l] << " e:" << e << " size:" << size << " Ap[l] + e:"<<Ap[l] + e<<"\n";
  P_ASSERT_X5(l < n, "WrLi bound: n=", n, ", l=", l, ".");
  P_ASSERT_X5(e < Ap[l+1]-Ap[l], "WrLi bound: e=", e, " Ap[l+1]-Ap[l]=", Ap[l+1]-Ap[l], ".");
  P_ASSERT_X3(e >= 0, "WrLi bound: e=", e, ".");
  P_ASSERT_X3(l >= 0, "WrLi bound: l=", l, ".");
  P_ASSERT_X7(Ap[l] + e < size, "WrLi bound: Ap[l]=", Ap[l], ", e=", e, ", size=", size, ".");
#endif
  return Ai[Ap[l] + e];
//  std::cout<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n";
}

inline double& SpMatrix::WrLx(int l, int e)
{
#ifdef DEBUG
  P_ASSERT_X5(l < n, "WrLi bound: n=", n, ", l=", l, ".");
  P_ASSERT_X5(e < Ap[l+1]-Ap[l], "WrLi bound: e=", e, " Ap[l+1]-Ap[l]=", Ap[l+1]-Ap[l], ".");
  P_ASSERT_X3(e >= 0, "WrLi bound: e=", e, ".");
  P_ASSERT_X3(l >= 0, "WrLi bound: l=", l, ".");
  P_ASSERT_X7(Ap[l] + e < size, "WrLi bound: Ap[l]=", Ap[l], ", e=", e, ", size=", size, ".");
#endif
  return Ax[Ap[l] + e];
//  std::cout<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n";
}

inline int SpMatrix::GetL(int n_)
{
#ifdef DEBUG
  P_ASSERT_X(n_ < n, "SpMatrix::GetL: Error\n");
#endif
  return Ap[n_+1] - Ap[n_];
}

inline void SpMatrix::AX(double* out, const double* in, double alpha, bool trans) const
{
  mmx(trans ? Trans : NoTrans, out, in, alpha);
}

inline void SpMatrix::AXpY(double* out, const double* in, const double* Y, double alpha, double beta, bool trans) const
{
  mmxpy(trans ? Trans : NoTrans, out, in, alpha, Y, beta);
}

// End of implementation of SpMatrix


// Implementation of Vector

inline Vector& Vector::operator=(const __op_mul_vec<SpMatrix, Vector> R)
{
  R.op.vec.mmx(R.op.tr, this->v, R.vecA.vec.v, R.op.alpha);
//  std::cout<<"__op_mul_vec<Matrix,Vector>\n";
  return *this;
}

inline Vector& Vector::operator=(const __op_mul_vec_plus_vec<SpMatrix, Vector> R)
{
  R.op.vec.mmxpy(R.op.tr, this->v, R.vecA.vec.v, R.op.alpha, R.vecB.vec.v, R.vecB.alpha);
//  std::cout<<"__op_mul_vec_plus_vec<Matrix,Vector>\n";
  return *this;
}

inline Matrix& Matrix::operator=(const __op_mul_vec<SpMatrix, Matrix> R)
{
  R.op.vec.mmxm(R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecA.vec.c);
//  std::cout<<"__op_mul_vec<Matrix,Matrix>\n";
  return *this;
}

inline Matrix& Matrix::operator=(const __op_mul_vec_plus_vec<SpMatrix, Matrix> R)
{
  R.op.vec.mmxmpym(R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecB.vec.m, R.vecB.vec.r, R.vecB.alpha, R.vecB.vec.c);
//  std::cout<<"__op_mul_vec_plus_vec<Matrix,Matrix>\n";
  return *this;
}

#endif
