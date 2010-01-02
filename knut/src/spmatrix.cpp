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
// **********************        SpMatrix Class        **********************//
//                                                                           //
// **************************************************************************//

void SpMatrix::Check()
{
  bool ok = true;
  for (int j = 0; j < n && ok; j++)
  {
    for (int k = Ap[j] + 1; k < Ap[j+1] && ok; k++)
    {
      if (Ai[k-1] >= Ai[k] || Ai[k] < 0 || Ai[k-1] < 0)
      {
        std::cerr << "SpMatrix::Check() in row=" << j << ", col=(" << Ai[k-1] << ">=" << Ai[k] << ")\n";
        ok = false;
      }
    }
  }
  for (int j = 0; j < n && ok; j++)
  {
    for (int k = Ap[j]; k < Ap[j+1] && ok; k++)
    {
      if (!isfinite(Ax[k]))
      {
        std::cerr << "SpMatrix::Check(): NaN at (" << j << "," << Ai[k] << ")\n";
        ok = false;
      }
    }
  }
  if (ok) std::cerr << "Matrix tested VALID.\n";
  else std::cerr << "Matrix tested INVALID. (Only the first error is shown.)\n";
  std::cerr.flush();
}

void SpMatrix::Swap()
{
  int *Rp = new int[n+1];
  int *Ri = new int[size];
  double *Rx = new double[size];

  umfpack_di_transpose(n, m, Ap, Ai, Ax, 0, 0, Rp, Ri, Rx);

  delete []Ap;
  delete []Ai;
  delete []Ax;

  Ap = Rp;
  Ai = Ri;
  Ax = Rx;

  if (format == 'C') format = 'R';
  else format = 'C';
  int nt = n;
  n = m;
  m = nt;
}

void SpMatrix::StrPlot(GnuPlot& pl)
{
  pl.SetPointSize(0.8);
  pl.Plot(0, "with points");

  for (int i = 0; i < n; i++)
  {
    for (int j = Ap[i]; j < Ap[i+1]; j++)
    {
      if (Ax[j] != 0.0) pl.AddData(0, (double)i, (double)Ai[j]);
    }
  }
  pl.Show();
}

void SpMatrix::Print()
{
  for (int i = 0; i < n; i++)
  {
    for (int j = Ap[i]; j < Ap[i+1]; j++)
    {
      cout << "(" << i << "," << Ai[j] << ")=" << Ax[j] << '\t';
      cout.flush();
    }
    cout << "\n";
  }
}

// **************************************************************************//
//                                                                           //
// **********************         SpFact Class         **********************//
//                                                                           //
// **************************************************************************//

void SpFact::init(int nn_)
{
  fact = false;
  Numeric = 0;
  Wi = new int[5*nn_+1];
  W  = new double[10*nn_+1];

  umfpack_di_defaults(Control);
  Control[UMFPACK_IRSTEP] = 0;
  Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
//  Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
}

SpFact::SpFact(char F, int nn_, int mm_, int nz) : SpMatrix(F, nn_, mm_, nz)
{
  init(nn_);
}

SpFact::SpFact(char F, int nn_, int nz) : SpMatrix(F, nn_, nn_, nz)
{
  init(nn_);
}

SpFact::SpFact(SpMatrix& M) : SpMatrix(M)
{
  init(n);
}

SpFact::~SpFact()
{
  if (Numeric != 0) umfpack_di_free_numeric(&Numeric);
  Numeric = 0;

  delete[] W;
  delete[] Wi;
}

void SpFact::clear()
{
  SpMatrix::clear();
  if (Numeric != 0)
  {
    P_ASSERT_X(fact, "Matrix was not factorized.");
    umfpack_di_free_numeric(&Numeric);
    Numeric = 0;
  }
  fact = false;
}

void SpFact::clear(char F)
{
  SpMatrix::clear(F);
  if (Numeric != 0)
  {
    P_ASSERT_X(fact, "Matrix was not factorized.");
    umfpack_di_free_numeric(&Numeric);
    Numeric = 0;
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

void SpFact::Fact()
{
  if (!fact)
  {
    void  *Symbolic = 0;
    P_ASSERT_X(Numeric == 0, "This is a bug. The sparse matrix was already factorized.");

    int   status = umfpack_di_symbolic(n, n, this->Ap, this->Ai, this->Ax, &Symbolic, Control, 0);
    P_ERROR_X4(status == UMFPACK_OK, "Error report from 'umfpack_di_symbolic()': (", status, ") ", sp_umf_error(status));
    status = umfpack_di_numeric(this->Ap, this->Ai, this->Ax, Symbolic, &Numeric, Control, 0);
    fact = true;
    if (status != UMFPACK_OK)
    {
      if (status == UMFPACK_WARNING_singular_matrix)
      {
        std::cerr << "The diagonal of the factorized matrix is being dumped:\n";
        Vector DX(n);
        GetDX(DX);
        DX.Print();
      }
      P_ERROR_X4(status == UMFPACK_OK, "Error report from 'umfpack_di_numeric()': (", status, ") ", sp_umf_error(status));
    }
    umfpack_di_free_symbolic(&Symbolic);
    Symbolic = 0;
  }
}

void SpFact::Solve(double* x, double* b, bool trans)
{
  if (!fact) Fact();
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
  status = umfpack_di_wsolve(sys, Ap, Ai, Ax, x, b, Numeric, Control, 0, Wi, W);
}

void SpFact::Solve(Vector& x, const Vector& b, bool trans)
{
  if (!fact) Fact();
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

  status = umfpack_di_wsolve(sys, Ap, Ai, Ax, x.Pointer(), b.v, Numeric, Control, 0, Wi, W);
  P_ERROR_X2(status == UMFPACK_OK, "Error report from 'umfpack_di_numeric()': ", sp_umf_error(status));
}

void SpFact::Solve(Matrix& x, const Matrix &b, bool trans)
{
  if (!fact) Fact();
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

  for (int i = 0; i < b.c; i++)
  {
    status = umfpack_di_wsolve(sys, Ap, Ai, Ax, x.Pointer(0, i), b.Pointer(0, i), Numeric, Control, 0, Wi, W);
    P_ERROR_X2(status == UMFPACK_OK, "Error report from 'umfpack_di_numeric()': ", sp_umf_error(status));
  }
}

// for sorting the eigenvalues
struct eigvalcomp : public std::binary_function<int, int, bool>
{
  const double* re;
  const double* im;
  eigvalcomp(double *re_, double *im_) : re(re_), im(im_) {}
  bool operator()(int i, int j)
  {
    return re[i]*re[i] + im[i]*im[i] < re[j]*re[j] + im[j]*im[j];
  }
};

void StabMatrix::Eigval(Vector& wr, Vector& wi)
{
  P_ERROR_X1((wr.size() == wi.size()) && (wr.size() > 1), "Real and imaginary vectors are not of the same size.");

  int     IDO      = 0;
  char    BMAT     = 'I';
  int     N        = AI.size() * A0.Col();
  char    WHICH[]  = {'L', 'M'};
  int     NEV      = wr.size() - 1;
  double  TOL      = knut_dlamch("EPS", 3);
//  RESID is defined in the class
  int     NCV      = 2 * NEV + 1;
  double* V        = new double[(N+1)*(NCV+1)];
  int     LDV      = N;
  int IPARAM[] = {1,     // ISHIFT
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
  int IPNTR[14];
  double* WORKD     = new double[4*N+1]; //3*N
  int     LWORKL    = 3 * NCV * NCV + 6 * NCV;
  double* WORKL     = new double[LWORKL+1];
  int INFO;
  // whether we need to create a new RESID or we can use it from a previous run
  if (isINIT) INFO = 1;
  else INFO = 0;

  double* tvec = new double[A0.Col()+1];
  double* tvec2 = new double[A0.Col()+1];

  // int it=0;
  do
  {
    knut_dnaupd(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV,
                IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO, 1, 2);
    if ((IDO == 1) || (IDO == -1))
    {
      double* in = &(WORKD[IPNTR[0] - 1]);
      double* out = &(WORKD[IPNTR[1] - 1]);
      // these are the unit operators above the diagonal
      for (int i = 0; i < AI.size() - 1; i++)
      {
        for (int j = 0; j < A0.Col(); j++) out[j + i*A0.Col()] = in[j + (i+1)*A0.Col()];
      }
      // the last row: multiplication and solution
      for (int i = 0; i < AI.size(); i++)
      {
        if (i == 0)
        {
          AI(AI.size() - 1).AX(tvec, in, 1.0, false);
        }
        else
        {
          if (AI(AI.size() - 1 - i).GetNZ() != 0)
          {
            AI(AI.size() - 1 - i).AX(tvec2, in + i*A0.Col(), 1.0, false);
            for (int j = 0; j < A0.Col(); j++) tvec[j] += tvec2[j];
          }
        }
      }
      A0.Solve(out + (AI.size() - 1)*A0.Col(), tvec);
    }
  }
  while ((IDO == 1) || (IDO == -1));
  P_ERROR_X3(IDO == 99, "IDO=", (int)IDO, " is not expected.");
  delete[] tvec2;
  delete[] tvec;

  const int NCONV = NEV; // IPARAM[4];
  int      RVEC     = 0;
  char     HOWMNY   = 'A';
  bool*    SELECT   = new bool[NCV+1];
  double*  DR       = new double[NCONV+2];
  double*  DI       = new double[NCONV+2];
  double*  Z        = new double[N*(NCONV+2)];
  int      LDZ      = N;
  double   SIGMAR   = 0.0;
  double   SIGMAI   = 0.0;
  double*  WORKEV   = new double[4*NCV]; // 3*NCV

//   std::cout<<"INFO1:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
//   std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();

  knut_dneupd(&RVEC, &HOWMNY, SELECT, DR, DI, Z, &LDZ,
    &SIGMAR, &SIGMAI, WORKEV,
    &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV,
    IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO, 1, 1, 1);

//   std::cout<<"INFO2:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
//   std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();

  wr.clear();
  wi.clear();
  // sorting eigenvalues;
  int* sortindex = new int[IPARAM[4]+1];
  for (int i = 0; i < IPARAM[4]; i++) sortindex[i] = i;
  eigvalcomp aa(DR, DI);
  std::sort(sortindex, sortindex + IPARAM[4], aa);
  const int eigstart = wr.size() - IPARAM[4];
  for (int i = 0; i < IPARAM[4]; i++)
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
