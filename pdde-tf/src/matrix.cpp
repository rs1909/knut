// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "plot.h"

using namespace std;

// **************************************************************************//
//                                                                           //
// **********************          Matrix Class        **********************//
//                                                                           //
// **************************************************************************//


void Matrix::StrPlot(GnuPlot& pl)
{
  pl.Plot(0, "with points ps 0.1");

  for (int i = 0; i < this->r; i++)
  {
    for (int j = 0; j < this->c; j++)
    {
      if ((*this)(i, j) != 0.0)
      {
        pl.AddData(0, i, j);
      }
    }
  }
  pl.Show();
}

void Matrix::Eigval(Vector& wr, Vector& wi)
{
  if ((this->r != this->c) || (this->r != wr.Size()) || (this->r != wi.Size()))
  {
    cout << "eigval: wrong dimensions" << '\n';
    exit(-1);
  }
  const Matrix& AA = *this;

  if (this->r == 1)
  {
    wr(0) = AA(0, 0);
    wi(0) = 0.0;
    return;
  }
  if (this->r == 2)
  {
    double tr = AA(0, 0) + AA(1, 1), det = AA(0, 0) * AA(1, 1) - AA(0, 1) * AA(1, 0);
    if (tr*tr - 4.0*det >= 0)
    {
      wr(0) = 0.5 * (tr - sqrt(tr * tr - 4.0 * det));
      wi(0) = 0.0;
      wr(1) = 0.5 * (tr + sqrt(tr * tr - 4.0 * det));
      wi(1) = 0.0;
    }
    else
    {
      wr(0) = 0.5 * tr;
      wi(0) = 0.5 * sqrt(-tr * tr + 4.0 * det);
      wr(1) = 0.5 * tr;
      wi(1) = -0.5 * sqrt(-tr * tr + 4.0 * det);
    }
    return;
  }

  char    jobvl = 'N', jobvr = 'N';
  Matrix  vl(this->r, 1), vr(this->r, 1);
  int n = this->r, lda = this->r;
  int ldvl = 1, ldvr = 1;
  Matrix  work(this->r, 4);
  int lwork = 4 * (this->r);
  int info;

  knut_dgeev(&jobvl, &jobvr, &n, this->m, &lda,
             wr.Pointer(), wi.Pointer(), vl.m, &ldvl, vr.m, &ldvr,
             work.m, &lwork, &info, 1, 1);
  if (info != 0) cout << "Error code" << info << '\n';
}

void Matrix::Eigval(Vector& wr, Vector& wi, Matrix& vl, Matrix& vr)
{
  if ((this->r != this->c) || (this->r != wr.Size()) || (this->r != wi.Size()))
  {
    cout << "eigval: wrong dimensions" << '\n';
    exit(-1);
  }

  if (this->r == 1)
  {
    wr(0) = *(this->m);
    wi(0) = 0.0;
    vl(0, 0) = 1.0;
    vr(0, 0) = 1.0;
    return;
  }

  char jobvl = 'V', jobvr = 'V';
  int  n = this->r, lda = this->r;
  int  ldvl = vl.Row(), ldvr = vr.Row();
  Matrix work(this->r, 10);
  int  lwork = 10 * (this->r);
  int  info;

  knut_dgeev(&jobvl, &jobvr, &n, this->m, &lda,
             wr.Pointer(), wi.Pointer(), vl.m, &ldvl, vr.m, &ldvr,
             work.m, &lwork, &info, 1, 1);
  if (info != 0) cout << "Error code" << info << '\n';
}

// MatFact routines

void MatFact::Fact()
{
  P_ASSERT_X(this->r == this->c, "MatFact::Fact wrong dimensions");

  char FACT = 'N', trans = 'N', equed = 'N';
  int n = this->r;
  int nrhs = 0;

  knut_dgesvx(&FACT, &trans, &n, &nrhs, this->m, &n, this->mf, &n,
              ipiv, &equed, NULL, NULL, NULL, &n, NULL, &n,
              &rcond, NULL, NULL, work, iwork, &info, 1, 1, 1);

  // cout<<"work(1): "<<work[0]<<" rcond< "<<rcond<<"\n";
  if (info != 0) cout << "MatFact::Fact: Error code" << info << '\n';
  fact = true;
}

void MatFact::Solve(Vector& x, const Vector& b, bool trans)
{
  P_ASSERT_X(this->r == this->c, "MatFact::Solve not a square matrix\n");
  P_ASSERT_X(b.Size() == this->c, "MatFact::Solve Vector sizes differ: M-V\n");
  P_ASSERT_X(b.Size() == x.Size(), "MatFact::Solve Vector sizes differ: V-V\n");

  const Matrix& AA = *this;
  double det;
  switch (this->r)
  {
    case 1:
      x(0) = b(0) / AA(0, 0);
      break;
    case 2:
      det = 1.0 / (AA(0, 0) * AA(1, 1) - AA(0, 1) * AA(1, 0));
      if (!trans)
      {
        x(0) = det * (AA(1, 1) * b(0) - AA(0, 1) * b(1));
        x(1) = det * (-AA(1, 0) * b(0) + AA(0, 0) * b(1));
      }
      else
      {
        x(0) = det * (AA(1, 1) * b(0) - AA(1, 0) * b(1));
        x(1) = det * (-AA(0, 1) * b(0) + AA(0, 0) * b(1));
      }
      break;
    default:
      if (!fact) Fact();

      char FACT = 'F', equed = 'N', TRANS;
      if (trans) TRANS = 'T';
      else TRANS = 'N';
      int n = this->r;
      int nrhs = 1;

      knut_dgesvx(&FACT, &TRANS, &n, &nrhs, this->m, &n, this->mf, &n,
                  ipiv, &equed, NULL, NULL, b.v, &n, x.Pointer(), &n,
                  &rcond, &ferr, &berr, work, iwork, &info, 1, 1, 1);
      if (info != 0) cout << "MatFact::Solve: Error code" << info << '\n';
      break;
  }
//  Vector err(b);
//  AXpY( err, x, b, 1.0, -1.0, trans );
//  cout<<"dim: "<<this->r<<" b: "<<sqrt(b*b)<<" e: "<<sqrt(err*err)<<"\n";
}

void MatFact::Solve(Matrix& x, const Matrix& b, bool trans)
{
  P_ASSERT_X(this->r == this->c, "MatFact::Solve not a square matrix\n");
  P_ASSERT_X(b.Col() != x.Col(), "MatFact::Solve Matrix columns differ\n");
  P_ASSERT_X((b.Row() == this->c) && (b.Row() == x.Row()), "MatFact::Solve Matrix rows differ\n");

  const Matrix& AA = *this;
  double det;
  switch (this->r)
  {
    case 1:
      for (int i = 0; i < b.Col(); i++)
      {
        x(0, i) = b(0, i) / AA(0, 0);
      }
      break;
    case 2:
      det = 1.0 / (AA(0, 0) * AA(1, 1) - AA(0, 1) * AA(1, 0));
      if (!trans)
      {
        for (int i = 0; i < b.Col(); i++)
        {
          x(0, i) = det * (AA(1, 1) * b(0, i) - AA(0, 1) * b(1, i));
          x(1, i) = det * (-AA(1, 0) * b(0, i) + AA(0, 0) * b(1, i));
        }
      }
      else
      {
        for (int i = 0; i < b.Col(); i++)
        {
          x(0, i) = det * (AA(1, 1) * b(0, i) - AA(1, 0) * b(1, i));
          x(1, i) = det * (-AA(0, 1) * b(0, i) + AA(0, 0) * b(1, i));
        }
      }
      break;
    default:
      if (!fact) Fact();

      char FACT = 'F', equed = 'N', TRANS;
      if (trans) TRANS = 'T';
      else TRANS = 'N';
      double *fberr = new double[2*b.Col()+1];
      int n = this->r;
      int nrhs = b.Col();

      knut_dgesvx(&FACT, &TRANS, &n, &nrhs, this->m, &n, this->mf, &n,
                  ipiv, &equed, NULL, NULL, b.m, &n, x.m, &n,
                  &rcond, fberr, fberr + b.Col(), work, iwork, &info, 1, 1, 1);

      if (info != 0) cout << "MatFact::Solve: Error code" << info << '\n';
      delete[] fberr;
      break;
  }
}
