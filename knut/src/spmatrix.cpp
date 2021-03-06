// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages packages root directory
//
// ------------------------------------------------------------------------- //

#include <iostream>
#include "spmatrix.h" // itt mar meghivtuk <vector>-t es "matrix.h"-t is
#include <cmath>
#include <functional>
#include <algorithm>

using namespace std;

// **************************************************************************//
//                                                                           //
// **********************        KNSparseMatrix Class        **********************//
//                                                                           //
// **************************************************************************//

void KNSparseMatrix::check()
{
  bool ok = true;
  for (size_t j = 0; j < n && ok; j++)
  {
    for (int k = Ap[j] + 1; k < Ap[j+1] && ok; k++)
    {
      if (Ai[k-1] >= Ai[k] || Ai[k] < 0 || Ai[k-1] < 0)
      {
        std::cerr << "KNSparseMatrix::check() in row=" << j << ", col=(" << Ai[k-1] << ">=" << Ai[k] << ")\n";
        ok = false;
      }
    }
  }
  for (size_t j = 0; j < n && ok; j++)
  {
    for (int k = Ap[j]; k < Ap[j+1] && ok; k++)
    {
      if (!isfinite(Ax[k]))
      {
        std::cerr << "KNSparseMatrix::check(): NaN at (" << j << "," << Ai[k] << ")\n";
        ok = false;
      }
    }
  }
  if (ok) std::cerr << "KNMatrix tested VALID.\n";
  else std::cerr << "KNMatrix tested INVALID. (Only the first error is shown.)\n";
  std::cerr.flush();
}

void KNSparseMatrix::swap()
{
  auto *Rp = new int[n+1];
  auto *Ri = new int[size];
  auto *Rx = new double[size];

  umfpack_di_transpose((int)n, (int)m, Ap, Ai, Ax, nullptr, nullptr, Rp, Ri, Rx);

  delete []Ap;
  delete []Ai;
  delete []Ax;

  Ap = Rp;
  Ai = Ri;
  Ax = Rx;

  if (format == 'C') format = 'R';
  else format = 'C';
  size_t nt = n;
  n = m;
  m = nt;
}

void KNSparseMatrix::sparsityPlot(GnuPlot& pl)
{
  pl.SetPointSize(0.8);
  pl.Plot(0, "with points");

  for (size_t i = 0; i < n; i++)
  {
    for (int j = Ap[i]; j < Ap[i+1]; j++)
    {
      if (Ax[j] != 0.0) pl.AddData(0, (double)i, (double)Ai[j]);
    }
  }
  pl.Show();
}

void KNSparseMatrix::print(std::ostream& os)
{
  for (size_t i = 0; i < n; i++)
  {
    for (int j = Ap[i]; j < Ap[i+1]; j++)
    {
      os << "(" << i << "," << Ai[j] << ")=" << Ax[j] << '\t';
      os.flush();
    }
    os << "\n";
  }
}

// **************************************************************************//
//                                                                           //
// **********************         KNLuSparseMatrix Class         **********************//
//                                                                           //
// **************************************************************************//

void KNLuSparseMatrix::init(size_t nn_)
{
  fact = false;
  Numeric = nullptr;
  Wi = new int[5*nn_+1];
  W  = new double[10*nn_+1];

  umfpack_di_defaults(Control);
  Control[UMFPACK_IRSTEP] = 0;
  Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
//  Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
}

KNLuSparseMatrix::KNLuSparseMatrix(char F, size_t nn_, size_t mm_, size_t nz) : KNSparseMatrix(F, nn_, mm_, nz)
{
  init(nn_);
}

KNLuSparseMatrix::KNLuSparseMatrix(char F, size_t nn_, size_t nz) : KNSparseMatrix(F, nn_, nn_, nz)
{
  init(nn_);
}

KNLuSparseMatrix::KNLuSparseMatrix(KNSparseMatrix& M) : KNSparseMatrix(M)
{
  init(n);
}

KNLuSparseMatrix::~KNLuSparseMatrix()
{
  if (Numeric != nullptr) umfpack_di_free_numeric(&Numeric);
  Numeric = nullptr;

  delete[] W;
  delete[] Wi;
}

void KNLuSparseMatrix::clear()
{
  KNSparseMatrix::clear();
  if (Numeric != nullptr)
  {
    P_ASSERT_X(fact, "KNMatrix was not factorized.");
    umfpack_di_free_numeric(&Numeric);
    Numeric = nullptr;
  }
  fact = false;
}

void KNLuSparseMatrix::clear(char F)
{
  KNSparseMatrix::clear(F);
  if (Numeric != nullptr)
  {
    P_ASSERT_X(fact, "KNMatrix was not factorized.");
    umfpack_di_free_numeric(&Numeric);
    Numeric = nullptr;
  }
  fact = false;
}

static const char* sp_umf_error(int status)
{
  switch (status)
  {
    case UMFPACK_OK:
      return "No problem.";
      break; 
    case UMFPACK_WARNING_singular_matrix:
      return "Singular matrix.";
      break;
    case UMFPACK_ERROR_out_of_memory:
      return "Out of memory.";
      break;
    case UMFPACK_ERROR_invalid_Numeric_object:
      return "Invalid Numeric object.";
      break;
    case UMFPACK_ERROR_invalid_Symbolic_object:
      return "Invalid Symbolic object.";
      break;
    case UMFPACK_ERROR_argument_missing:
      return "Missing argument.";
      break;
    case UMFPACK_ERROR_n_nonpositive:
      return "n is less than or equal to zero.";
      break;
    case UMFPACK_ERROR_invalid_matrix:
      return "Invalid matrix.";
      break;
    case UMFPACK_ERROR_different_pattern:
      return "Sparsity pattern has changed between calls.";
      break; 
    case UMFPACK_ERROR_invalid_system:
      return "Invalid system.";
      break;
    case UMFPACK_ERROR_invalid_permutation:
      return "Invalid permutation.";
      break;
    case UMFPACK_ERROR_internal_error:
      return "Internal Error.";
    case UMFPACK_ERROR_file_IO:
      return "File input/output error.";
      break;
    default:
      return "Unknown error.";
      break;
  }
}

void KNLuSparseMatrix::luFactorize()
{
  if (!fact)
  {
    void  *Symbolic = nullptr;
    P_ASSERT_X(Numeric == nullptr, "This is a bug. The sparse matrix was already factorized.");

    int   status = umfpack_di_symbolic((int)n, (int)n, this->Ap, this->Ai, this->Ax, &Symbolic, Control, nullptr);
    P_ERROR_X4(status == UMFPACK_OK, "Error report from 'umfpack_di_symbolic()': (", status, ") ", sp_umf_error(status));
    status = umfpack_di_numeric(this->Ap, this->Ai, this->Ax, Symbolic, &Numeric, Control, nullptr);
    fact = true;
    if (status != UMFPACK_OK)
    {
      if (status == UMFPACK_WARNING_singular_matrix)
      {
        std::cerr << "The diagonal of the factorized matrix is being dumped:\n";
        KNVector DX(n);
        GetDX(DX);
        DX.print(std::cout);
      }
      P_ERROR_X4(status == UMFPACK_OK, "Error report from 'umfpack_di_numeric()': (", status, ") ", sp_umf_error(status));
    }
    umfpack_di_free_symbolic(&Symbolic);
    Symbolic = nullptr;
  }
}

void KNLuSparseMatrix::solve(double* x, double* b, bool trans)
{
  if (!fact) luFactorize();
  int status, sys;

  if (format == 'C')
  {
    if (trans) sys = UMFPACK_At;
    else sys = UMFPACK_A;
  }
  else
  {
    if (trans) sys = UMFPACK_A;
    else sys = UMFPACK_At;
    P_ERROR_X1(format == 'R', "Unknown sparse matrix format.");
  }
  status = umfpack_di_wsolve(sys, Ap, Ai, Ax, x, b, Numeric, Control, nullptr, Wi, W);
  P_ERROR_X2(status == UMFPACK_OK, "Error report from 'umfpack_di_wsolve()': ", sp_umf_error(status));
}

void KNLuSparseMatrix::solve(KNVector& x, const KNVector& b, bool trans)
{
  if (!fact) luFactorize();
  int status, sys;

  if (format == 'C')
  {
    if (trans) sys = UMFPACK_At;
    else sys = UMFPACK_A;
  }
  else
  {
    if (trans) sys = UMFPACK_A;
    else sys = UMFPACK_At;
    P_ERROR_X1(format == 'R', "Unknown sparse matrix format.");
  }

  status = umfpack_di_wsolve(sys, Ap, Ai, Ax, x.pointer(), b.v, Numeric, Control, nullptr, Wi, W);
  P_ERROR_X2(status == UMFPACK_OK, "Error report from 'umfpack_di_wsolve()': ", sp_umf_error(status));
}

void KNLuSparseMatrix::solve(KNMatrix& x, const KNMatrix &b, bool trans)
{
  if (!fact) luFactorize();
  int status, sys;
  if (format == 'C')
  {
    if (trans) sys = UMFPACK_At;
    else sys = UMFPACK_A;
  }
  else
  {
    if (trans) sys = UMFPACK_A;
    else sys = UMFPACK_At;
    P_ERROR_X1(format == 'R', "Unknown sparse matrix format.");
  }

  for (size_t i = 0; i < b.c; i++)
  {
    status = umfpack_di_wsolve(sys, Ap, Ai, Ax, x.pointer(0, i), b.pointer(0, i), Numeric, Control, nullptr, Wi, W);
    P_ERROR_X2(status == UMFPACK_OK, "Error report from 'umfpack_di_numeric()': ", sp_umf_error(status));
  }
}

// for sorting the eigenvalues
struct eigvalcomp : public std::binary_function<size_t, size_t, bool>
{
  const double* re;
  const double* im;
  eigvalcomp(double *re_, double *im_) : re(re_), im(im_) {}
  bool operator()(size_t i, size_t j)
  {
    return re[i]*re[i] + im[i]*im[i] < re[j]*re[j] + im[j]*im[j];
  }
};

void KNSparseMatrixPolynomial::eigenvalues(KNVector& wr, KNVector& wi)
{
  P_ERROR_X1((wr.size() == wi.size()) && (wr.size() > 1), "Real and imaginary vectors are not of the same size.");

  blasint IDO      = 0;
  char    BMAT     = 'I';
  size_t  N        = AI.size() * A0.col();
  char    WHICH[]  = {'L', 'M'};
  size_t  NEV      = wr.size() - 1;
  double  TOL      = knut_dlamch("EPS", 3);
//  RESID is defined in the class
  size_t  NCV      = 2 * NEV + 1;
  auto* V        = new double[(N+1)*(NCV+1)];
  size_t  LDV      = N;
  blasint IPARAM[] = {1,     // ISHIFT
                  0,     // n.a.ff
                  300,   // MXITER
                  1,     // NB
                  0,     // NCONV (out)
                  0,     // IUPD n.a.
                  1,     // MODE
                  0,     // NP if there is shift ()
                  0,     // NUMOP
                  0,     // NUMOPB
                  0,     // NUMREO
                  0 };   // padding
  blasint IPNTR[14];
  auto* WORKD     = new double[4*N+1]; //3*N
  size_t  LWORKL    = 3 * NCV * NCV + 6 * NCV;
  auto* WORKL     = new double[LWORKL+1];
  blasint INFO;
  // whether we need to create a new RESID or we can use it from a previous run
  if (isINIT) INFO = 1;
  else INFO = 0;

  auto* tvec = new double[A0.col()+1];
  auto* tvec2 = new double[A0.col()+1];

  // int it=0;
  do
  {
    knut_dnaupd(&IDO, BMAT, N, WHICH[0], WHICH[1], NEV, &TOL, RESID, NCV, V, LDV,
                IPARAM, IPNTR, WORKD, WORKL, LWORKL, &INFO);
    if ((IDO == 1) || (IDO == -1))
    {
      double* in = &(WORKD[IPNTR[0] - 1]);
      double* out = &(WORKD[IPNTR[1] - 1]);
      // these are the unit operators above the diagonal
      for (size_t i = 0; i < AI.size() - 1; i++)
      {
        for (size_t j = 0; j < A0.col(); j++) out[j + i*A0.col()] = in[j + (i+1)*A0.col()];
      }
      // the last row: multiplication and solution
      for (size_t i = 0; i < AI.size(); i++)
      {
        if (i == 0)
        {
          AI(AI.size() - 1).timesX(tvec, in, 1.0, false);
        }
        else
        {
          if (AI(AI.size() - 1 - i).nonzeros() != 0)
          {
            AI(AI.size() - 1 - i).timesX(tvec2, in + i*A0.col(), 1.0, false);
            for (size_t j = 0; j < A0.col(); j++) tvec[j] += tvec2[j];
          }
        }
      }
      A0.solve(out + (AI.size() - 1)*A0.col(), tvec);
    }
  }
  while ((IDO == 1) || (IDO == -1));
  P_ERROR_X3(IDO == 99, "IDO=", (int)IDO, " is not expected.");
  delete[] tvec2;
  delete[] tvec;

  const size_t NCONV = NEV; // IPARAM[4];
//   blasint  RVEC     = 0;
  char     HOWMNY   = 'A';
  auto*    SELECT   = new bool[NCV+1];
  auto*  DR       = new double[NCONV+2];
  auto*  DI       = new double[NCONV+2];
  auto*  Z        = new double[N*(NCONV+2)];
  size_t   LDZ      = N;
  double   SIGMAR   = 0.0;
  double   SIGMAI   = 0.0;
  auto*  WORKEV   = new double[4*NCV]; // 3*NCV

//   std::cout<<"INFO1:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
//   std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();

  knut_dneupd(false, HOWMNY, SELECT, DR, DI, Z, LDZ,
    &SIGMAR, &SIGMAI, WORKEV,
    BMAT, N, WHICH[0], WHICH[1], NEV, &TOL, RESID, NCV, V, LDV,
    IPARAM, IPNTR, WORKD, WORKL, LWORKL, &INFO);

//   std::cout<<"INFO2:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
//   std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();

  wr.clear();
  wi.clear();
  // sorting eigenvalues;
  const auto evals = static_cast<size_t>(IPARAM[4]);
  auto* sortindex = new size_t[evals];
  for (size_t i = 0; i < evals; i++) sortindex[i] = i;
  eigvalcomp aa(DR, DI);
  std::sort(sortindex, sortindex + evals, aa);
  const size_t eigstart = wr.size() - evals;
  for (size_t i = 0; i < evals; i++)
  {
    wr(eigstart+i) = DR[sortindex[i]];
    wi(eigstart+i) = DI[sortindex[i]];
  }
  delete[] sortindex;

  // from the second
  delete[] WORKEV;
  delete[] Z;
  delete[] DI;
  delete[] DR;
  delete[] SELECT;
  // from the first
  delete[] WORKL;
  delete[] WORKD;
  delete[] V;
  // now we have a RESID
  isINIT = true;
}
