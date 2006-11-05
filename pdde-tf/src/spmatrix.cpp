// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the packages packages root directory
//
// ------------------------------------------------------------------------- //

#include <iostream>
#include "spmatrix.h" // itt mar meghivtuk <vector>-t
#include <cmath>

using namespace std;

extern "C"
{
  // LAPACK
  double pdde_dlamch(char *cmach, int cmach_len);
  /* ARPACK */
  int pdde_dnaupd(int *ido, char *bmat, int *n, char* which,
                  int *nev, double *tol, double *resid, int *ncv,
                  double *v, int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl, int *lworkl, int *info,
                  int bmat_len, int which_len);

  int pdde_dneupd(bool *rvec, char *howmny, bool *select,
                  double *dr, double *di, double *z__, int *ldz,
                  double *sigmar, double *sigmai, double *workev,
                  char *bmat, int *n, char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *v, int *ldv,
                  int *iparam, int *ipntr, double *workd, double *workl,
                  int *lworkl, int *info, int howmny_len, int bmat_len,
                  int which_len);
}

// **************************************************************************//
//                                                                           //
// **********************        SpMatrix Class        **********************//
//                                                                           //
// **************************************************************************//

void SpMatrix::Check()
{
  for (int j = 0; j < n; j++)
  {
    for (int k = Ap[j] + 1; k < Ap[j+1]; k++)
    {
      if (Ai[k-1] >= Ai[k]) cout << "SpMatrix::Check(): " << j << "," << Ai[k-1] << "\n";
    }
  }
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

void SpFact::Init(int nn_)
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
  Init(nn_);
}

SpFact::SpFact(char F, int nn_, int nz) : SpMatrix(F, nn_, nn_, nz)
{
  Init(nn_);
}

SpFact::SpFact(SpMatrix& M) : SpMatrix(M)
{
  Init(n);
}

SpFact::~SpFact()
{
  if (Numeric != 0) umfpack_di_free_numeric(&Numeric);
  Numeric = 0;

  delete[] W;
  delete[] Wi;
}

void SpFact::Clear()
{
  SpMatrix::Clear();
  if (Numeric != 0)
  {
    P_ASSERT_X(fact, "Matrix was not factorized.");
    umfpack_di_free_numeric(&Numeric);
    Numeric = 0;
  }
  fact = false;
}

void SpFact::Clear(char F)
{
  SpMatrix::Clear(F);
  if (Numeric != 0)
  {
    P_ASSERT_X(fact, "Matrix was not factorized.");
    umfpack_di_free_numeric(&Numeric);
    Numeric = 0;
  }
  fact = false;
}


void SpFact::Fact()
{
  if (!fact)
  {
    void  *Symbolic = 0;
    P_ASSERT_X(Numeric == 0, "This is a bug. The sparse matrix was already factorized.");

    int   status = umfpack_di_symbolic(n, n, this->Ap, this->Ai, this->Ax, &Symbolic, Control, 0);
    P_ERROR_X2(status == 0, "Error report from \"umfpack_di_symbolic()\":", status);
    status = umfpack_di_numeric(this->Ap, this->Ai, this->Ax, Symbolic, &Numeric, Control, 0);
    fact = true;
    if (status != 0)
    {
      if (status == 1)
      {
        std::cout << "The matrix is singular. "
        << "The diagonal of the factorized matrix is being dumped: "
        << status << "\n";
        Vector DX(n);
        GetDX(DX);
        DX.Print();
      }
      P_ERROR_X2(status == 0, "Error report from \"umfpack_di_numeric()\":", status);
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
    status = umfpack_di_wsolve(sys, Ap, Ai, Ax, x.m + i * x.r, b.m + i * b.r, Numeric, Control, 0, Wi, W);
  }
}

void StabMatrix::Eigval(Vector& wr, Vector& wi)
{
  P_ERROR_X1((wr.Size() == wi.Size()) && (wr.Size() > 1), "Real and imaginary vectors are not of the same size.");

  int     IDO      = 0;
  char    BMAT     = 'I';
  int     N        = AI.Size() * A0.Col();
  char    WHICH[]  = {'L', 'M'};
  int     NEV      = wr.Size() - 1;
  double  TOL      = pdde_dlamch("EPS", 3);
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
    pdde_dnaupd(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV,
                IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO, 1, 2);
    if ((IDO == 1) || (IDO == -1))
    {
      double* in = &(WORKD[IPNTR[0] - 1]);
      double* out = &(WORKD[IPNTR[1] - 1]);
      // these are the unit operators above the diagonal
      for (int i = 0; i < AI.Size() - 1; i++)
      {
        for (int j = 0; j < A0.Col(); j++) out[j + i*A0.Col()] = in[j + (i+1)*A0.Col()];
      }
      // the last row: multiplication and solution
      for (int i = 0; i < AI.Size(); i++)
      {
        if (i == 0)
        {
          AI(AI.Size() - 1).AX(tvec, in, 1.0, false);
        }
        else
        {
          if (AI(AI.Size() - 1 - i).GetNZ() != 0)
          {
            AI(AI.Size() - 1 - i).AX(tvec2, in + i*A0.Col(), 1.0, false);
            for (int j = 0; j < A0.Col(); j++) tvec[j] += tvec2[j];
          }
        }
      }
      A0.Solve(out + (AI.Size() - 1)*A0.Col(), tvec);
    }
  }
  while ((IDO == 1) || (IDO == -1));
  P_ERROR_X3(IDO == 99, "IDO = ", (int)IDO, " is not expected\n");
  delete[] tvec2;
  delete[] tvec;

  bool     RVEC     = false;
  char     HOWMNY   = 'A';
  bool*    SELECT   = new bool[NCV+1];
  double*  DR       = new double[NEV+2];
  double*  DI       = new double[NEV+2];
  double*  Z        = new double[N*(NEV+2)];
  int      LDZ      = N;
  double   SIGMAR   = 0.0;
  double   SIGMAI   = 0.0;
  double*  WORKEV   = new double[4*NCV]; // 3*NCV

//  std::cout<<"INFO1:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
//  std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();

  if (true /*IPARAM[4] >= NEV*/)
  {
    pdde_dneupd(&RVEC, &HOWMNY, SELECT, DR, DI, Z, &LDZ,
                &SIGMAR, &SIGMAI, WORKEV,
                &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV,
                IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO, 1, 1, 2);
//   std::cout<<"INFO2:"<<INFO<<" N: "<<N<<" NEV: "<<NEV<<'\n';
//   std::cout<<"converged: "<<IPARAM[4]<<"\n"; std::cout.flush();
    wr.Clear();
    wi.Clear();
    if (IPARAM[4] > wr.Size())
    {
      std::cout << "More eigenvalues were calculated than expected: N=" << IPARAM[4] << "\n";
      std::cout.flush();
      IPARAM[4] = wr.Size();
    }
    for (int i = 0; i < IPARAM[4]; i++)
    {
      wr(i) = DR[i];
      wi(i) = DI[i];
    }
  }

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
