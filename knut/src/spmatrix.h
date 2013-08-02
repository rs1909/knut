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
// **********************        KNSparseMatrix Class        **********************//
//                                                                           //
// **************************************************************************//

class KNSparseMatrix
{

  protected:

    char format;
    size_t n;
    size_t m;
    size_t size;
    int *Ap;
    int *Ai;
    double *Ax;

  public:

    inline KNSparseMatrix() : format('R'), n(0), m(0), size(0), Ap(0), Ai(0), Ax(0)
    { }
    inline KNSparseMatrix(char F, size_t n_, size_t m_, size_t nz)
    {
      init(F, n_, m_, nz);
    }
    inline KNSparseMatrix(char F, size_t n_, size_t nz)
    {
      init(F, n_, n_, nz);
    }
    inline KNSparseMatrix(const KNSparseMatrix& M)
    {
      init(M);
    }
    inline virtual ~KNSparseMatrix();

    inline void init(char F, size_t n_, size_t m_, size_t nz);
    inline void init(const KNSparseMatrix& M);

    inline void mmx(enum cspblas_Trans trans, double* out, const double* in, double alpha) const
    {
      cspblas_mmx(format, trans, n, m, Ap, Ai, Ax, out, in, alpha);
    }
    inline void mmxpy(enum cspblas_Trans trans, double* out, const double* in, double alpha, const double* Y, double beta) const
    {
      cspblas_mmxpy(format, trans, n, m, Ap, Ai, Ax, out, in, alpha, Y, beta);
    }
    inline void mmxm(enum cspblas_Trans trans, double* out, size_t ldout, const double* in, size_t ldin, double alpha, size_t nrhs) const
    {
      cspblas_mmxm(format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, nrhs);
    }
    inline void mmxmpym(enum cspblas_Trans trans, double* out, size_t ldout, const double* in, size_t ldin, double alpha,
                        const double* Y, size_t ldY, double beta, size_t nrhs) const
    {
      cspblas_mmxmpym(format, trans, n, m, Ap, Ai, Ax, out, ldout, in, ldin, alpha, Y, ldY, beta, nrhs);
    }

    inline __op_mul_vec<KNSparseMatrix, KNVector>     operator*(const KNVector& v) const
    {
      return __op_mul_vec<KNSparseMatrix, KNVector>(__scal_vec_trans<KNSparseMatrix>(*this), __scal_vec_trans<KNVector>(v));
    }
    inline __op_mul_vec<KNSparseMatrix, KNMatrix>     operator*(const KNMatrix& v) const
    {
      return __op_mul_vec<KNSparseMatrix, KNMatrix>(__scal_vec_trans<KNSparseMatrix>(*this), __scal_vec_trans<KNMatrix>(v));
    }
    inline __op_mul_vec_rng<KNSparseMatrix, KNVector> operator*(const __scal_vec_trans_rng<KNVector> v) const
    {
      return __op_mul_vec_rng<KNSparseMatrix, KNVector>(__scal_vec_trans_rng<KNSparseMatrix>(*this), v);
    }
    inline __op_mul_vec_rng<KNSparseMatrix, KNMatrix> operator*(const __scal_vec_trans_rng<KNMatrix> v) const
    {
      return __op_mul_vec_rng<KNSparseMatrix, KNMatrix>(__scal_vec_trans_rng<KNSparseMatrix>(*this), v);
    }
    inline __scal_vec_trans<KNSparseMatrix>        operator!() const
    {
      return __scal_vec_trans<KNSparseMatrix>(*this, 1.0, Trans);
    }
    inline __scal_vec_trans_rng<KNSparseMatrix>    operator[ ](const rng r) const
    {
      return __scal_vec_trans_rng<KNSparseMatrix>(*this, r);
    }

    inline void timesX(double* out, const double* in, double alpha, bool trans) const;
    inline void timesXPlusY(double* out, const double* in, const double* Y, double alpha, double beta, bool trans) const;

    inline size_t row() const
    {
      if (format == 'R') return n;
      else return m;
    }
    inline size_t col() const
    {
      if (format == 'R') return m;
      else return n;
    }

    /// clear the matrix
    /// these have to be virtual, because these might be called as KNLuSparseMatrix
    virtual void clear();
    virtual void clear(char F);
    /// checks the structure of the matrix
    void check();

    // Fill in routines
    /// Creates a new line in the matrix
    inline size_t newLine(size_t size_);
    /// Writes or returns the row or column index
    /// into line `l' and the `e'-th element
    inline int& writeIndex(size_t l, size_t e);
    /// Writes into line `l' and the `e'-th element
    inline double& writeData(size_t l, size_t e);
    /// returns the length of the n_ -th line in the matrix
    inline size_t lineLength(size_t n_);
    /// returns the nonzero elements in the matrix
    inline size_t nonzeros()
    {
      return static_cast<size_t>(Ap[n]);
    }
    /// returns the number of lines, e.g. columns or rows in the matrix depending on format
    inline size_t lines()
    {
      return n;
    }
    // Computation routines
    /// transposes the matrix into the other format
    void swap();
    // these are used for debugging only
    /// plots the structure of the matrix
    void sparsityPlot(GnuPlot& pl);
    /// prints out Ap
    void printAp()
    {
      for (size_t i = 0; i < n + 1; i++) std::cout << Ap[i] << '\t';
      std::cout << '\n';
    }
    /// prints the whole matrix onto the screen
    void print();
};

// **************************************************************************//
//                                                                           //
// **********************         KNLuSparseMatrix Class         **********************//
//                                                                           //
// **************************************************************************//

class KNLuSparseMatrix : public KNSparseMatrix
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
    void init(size_t nn_);
    KNLuSparseMatrix(char F, size_t n_, size_t m_, size_t nz);
    KNLuSparseMatrix(char F, size_t n_, size_t nz);
    KNLuSparseMatrix(KNSparseMatrix& M);
    inline KNLuSparseMatrix(KNLuSparseMatrix&) : KNSparseMatrix('F', 1, 1)
    {
      std::cout << "KNLuSparseMatrix::KNLuSparseMatrix(SpFact&): not implemented\n";
    }
    ~KNLuSparseMatrix();

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
    void solve(KNVector& x, const KNVector& b, bool trans = false);
    void solve(KNMatrix& x, const KNMatrix& b, bool trans = false);
    /// get the diagonal of the Upper factor
    int GetDX(KNVector& V)
    {
      if (!fact) luFactorize();
      return umfpack_di_get_numeric(0, 0, 0, 0, 0, 0, 0, 0, V.v, 0, 0, Numeric);
    }
  private:
    void luFactorize();
};

// void KNDdeBvpCollocation::jacobianOfStability( StabMat& AB, const KNVector& par, const JagMatrix3D& solData );
// this computes its own matrix structures, because it is nowhere else needed: kkSI, eeSI, rrSI, ddSI, etc.
// However the variables will be contained within the KNDdeBvpCollocation class

class KNSparseMatrixPolynomial
{
  public:

    KNSparseMatrixPolynomial(size_t nmat_, size_t n_, size_t nz) : A0('R', n_,  nz), AI(nmat_)
    {
      for (size_t i = 0; i < AI.size(); i++) AI(i).init('R', n_, n_, nz);
      RESID = new double[ AI.size() * n_ + 1 ];
      isINIT = false;
    }
    ~KNSparseMatrixPolynomial()
    {
      delete[] RESID;
    }

    size_t                nmat()
    {
      return AI.size();
    }
    KNLuSparseMatrix&            getA0()
    {
      return A0;
    }
    KNArray1D<KNSparseMatrix>& getAI()
    {
      return AI;
    }
    KNSparseMatrix&          getAI(size_t i)
    {
      return AI(i);
    }

    void eigenvalues(KNVector& wr, KNVector& wi);

  private:

    KNLuSparseMatrix             A0;
    KNArray1D<KNSparseMatrix>  AI;

    double* RESID;
    bool    isINIT;
};

// Implementation of KNSparseMatrix

inline void KNSparseMatrix::init(char F, size_t n_, size_t m_, size_t nz)
{
  if ((F != 'R') && (F != 'C')) std::cout << "KNSparseMatrix::CONSTRUCTOR: invalid format specification.\n";
  format = F;
  n = 0;
  m = m_;
  size = nz;
  Ap = new int[n_+1];
  Ai = new int[nz];
  Ax = new double[nz];
  for (size_t i = 0; i < nz; i++)
  {
    Ai[i] = 0;
    Ax[i] = 0.0;
  }
  for (size_t i = 0; i < n_ + 1; i++)
  {
    Ap[i] = 0;
  }
}

inline void KNSparseMatrix::init(const KNSparseMatrix& M)
{
  format = M.format;
  n = M.n;
  m = M.m;
  size = M.size;

  Ap = new int[M.n+1];
  Ai = new int[size];
  Ax = new double[size];
  for (size_t i = 0; i < size; i++)
  {
    Ai[i] = M.Ai[i];
    Ax[i] = M.Ax[i];
  }
  for (size_t i = 0; i < n + 1; i++)
  {
    Ap[i] = M.Ap[i];
  }
}

inline KNSparseMatrix::~KNSparseMatrix()
{
  delete []Ax;
  delete []Ai;
  delete []Ap;
}

inline void KNSparseMatrix::clear()
{

  for (size_t i = 0; i < n + 1; i++)
  {
    Ap[i] = 0;
  }
  n = 0;
  for (size_t i = 0; i < size; i++)
  {
    Ai[i] = 0;
    Ax[i] = 0.0;
  }
}

inline void KNSparseMatrix::clear(char F)
{

  for (size_t i = 0; i < n + 1; i++)
  {
    Ap[i] = 0;
  }
  n = 0;
  format = F;
  for (size_t i = 0; i < size; i++)
  {
    Ai[i] = 0;
    Ax[i] = 0.0;
  }
}

inline size_t KNSparseMatrix::newLine(size_t size_)
{
  n++;
  Ap[n] = Ap[n-1] + static_cast<int>(size_);
#ifdef DEBUG
  P_ASSERT_X5(static_cast<size_t>(Ap[n]) <= size, "WrLi bound: Ap[n]=", Ap[n], ", size=", size, ".");
#endif
  return n -1;
}

inline int& KNSparseMatrix::writeIndex(size_t l, size_t e)
{
#ifdef DEBUG
//   std::cout << "n:" << n << " l:" << l << " lsz:" << Ap[l+1] - Ap[l] << " e:" << e << " size:" << size << " Ap[l] + e:"<<Ap[l] + e<<"\n";
  P_ASSERT_X5(l < n, "WrLi bound: n=", n, ", l=", l, ".");
  P_ASSERT_X5(e < static_cast<size_t>(Ap[l+1]-Ap[l]), "WrLi bound: e=", e, " Ap[l+1]-Ap[l]=", Ap[l+1]-Ap[l], ".");
//  P_ASSERT_X3(e >= 0, "WrLi bound: e=", e, ".");
//  P_ASSERT_X3(l >= 0, "WrLi bound: l=", l, ".");
  P_ASSERT_X7(Ap[l] + e < size, "WrLi bound: Ap[l]=", Ap[l], ", e=", e, ", size=", size, ".");
#endif
  return Ai[static_cast<size_t>(Ap[l]) + e];
//  std::cout<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n";
}

inline double& KNSparseMatrix::writeData(size_t l, size_t e)
{
#ifdef DEBUG
  P_ASSERT_X5(l < n, "WrLi bound: n=", n, ", l=", l, ".");
  P_ASSERT_X5(e < static_cast<size_t>(Ap[l+1]-Ap[l]), "WrLi bound: e=", e, " Ap[l+1]-Ap[l]=", Ap[l+1]-Ap[l], ".");
//  P_ASSERT_X3(e >= 0, "WrLi bound: e=", e, ".");
//  P_ASSERT_X3(l >= 0, "WrLi bound: l=", l, ".");
  P_ASSERT_X7(static_cast<size_t>(Ap[l] + e) < size, "WrLi bound: Ap[l]=", Ap[l], ", e=", e, ", size=", size, ".");
#endif
  return Ax[static_cast<size_t>(Ap[l]) + e];
//  std::cout<<l<<","<<n<<"-"<<Ap[l] + e<<","<<Ap[l+1]<<"\n";
}

inline size_t KNSparseMatrix::lineLength(size_t n_)
{
#ifdef DEBUG
  P_ASSERT_X(n_ < n, "KNSparseMatrix::GetL: Error\n");
#endif
  return static_cast<size_t>(Ap[n_+1] - Ap[n_]);
}

inline void KNSparseMatrix::timesX(double* out, const double* in, double alpha, bool trans) const
{
  mmx(trans ? Trans : NoTrans, out, in, alpha);
}

inline void KNSparseMatrix::timesXPlusY(double* out, const double* in, const double* Y, double alpha, double beta, bool trans) const
{
  mmxpy(trans ? Trans : NoTrans, out, in, alpha, Y, beta);
}

// End of implementation of KNSparseMatrix


// Implementation of KNVector

inline KNVector& KNVector::operator=(const __op_mul_vec<KNSparseMatrix, KNVector> R)
{
  R.op.vec.mmx(R.op.tr, this->v, R.vecA.vec.v, R.op.alpha);
//  std::cout<<"__op_mul_vec<KNMatrix,KNVector>\n";
  return *this;
}

inline KNVector& KNVector::operator=(const __op_mul_vec_plus_vec<KNSparseMatrix, KNVector> R)
{
  R.op.vec.mmxpy(R.op.tr, this->v, R.vecA.vec.v, R.op.alpha, R.vecB.vec.v, R.vecB.alpha);
//  std::cout<<"__op_mul_vec_plus_vec<KNMatrix,KNVector>\n";
  return *this;
}

inline KNMatrix& KNMatrix::operator=(const __op_mul_vec<KNSparseMatrix, KNMatrix> R)
{
  R.op.vec.mmxm(R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecA.vec.c);
//  std::cout<<"__op_mul_vec<KNMatrix,Matrix>\n";
  return *this;
}

inline KNMatrix& KNMatrix::operator=(const __op_mul_vec_plus_vec<KNSparseMatrix, KNMatrix> R)
{
  R.op.vec.mmxmpym(R.op.tr, this->m, this->r, R.vecA.vec.m, R.vecA.vec.r, R.op.alpha, R.vecB.vec.m, R.vecB.vec.r, R.vecB.alpha, R.vecB.vec.c);
//  std::cout<<"__op_mul_vec_plus_vec<KNMatrix,Matrix>\n";
  return *this;
}

#endif
