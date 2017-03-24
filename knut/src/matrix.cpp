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
// **********************          KNMatrix Class        **********************//
//                                                                           //
// **************************************************************************//


void KNMatrix::sparsityPlot(GnuPlot& pl)
{
  pl.Plot(0, "with points ps 0.1");

  for (size_t i = 0; i < this->r; i++)
  {
    for (size_t j = 0; j < this->c; j++)
    {
      if ((*this)(i, j) != 0.0)
      {
        pl.AddData(0, i, j);
      }
    }
  }
  pl.Show();
}

void KNMatrix::eigenvalues(KNVector& wr, KNVector& wi)
{
  if ((this->r != this->c) || (this->r != wr.size()) || (this->r != wi.size()))
  {
    cout << "eigval: wrong dimensions" << '\n';
    exit(-1);
  }
  const KNMatrix& AA = *this;

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
  KNMatrix  vl(this->r, 1), vr(this->r, 1);
  ptrdiff_t n = this->r, lda = this->r;
  ptrdiff_t ldvl = 1, ldvr = 1;
  KNMatrix  work(this->r, 4);
  ptrdiff_t lwork = 4 * (this->r);
  blasint info;

  knut_dgeev(jobvl, jobvr, n, this->m, lda,
             wr.pointer(), wi.pointer(), vl.m, ldvl, vr.m, ldvr,
             work.m, lwork, &info);
  if (info != 0) cout << "Error code" << info << '\n';
}

void KNMatrix::eigenvalues(KNVector& wr, KNVector& wi, KNMatrix& vl, KNMatrix& vr)
{
  if ((this->r != this->c) || (this->r != wr.size()) || (this->r != wi.size()))
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
  ptrdiff_t  n = this->r, lda = this->r;
  ptrdiff_t  ldvl = vl.row(), ldvr = vr.row();
  KNMatrix work(this->r, 10);
  ptrdiff_t  lwork = 10 * (this->r);
  blasint  info;

  knut_dgeev(jobvl, jobvr, n, this->m, lda,
             wr.pointer(), wi.pointer(), vl.m, ldvl, vr.m, ldvr,
             work.m, lwork, &info);
  if (info != 0) cout << "Error code" << info << '\n';
}

// KNLuMatrix routines

void KNLuMatrix::luFactorize()
{
  P_ASSERT_X(this->r == this->c, "KNLuMatrix::Fact wrong dimensions");

  char FACT = 'N', trans = 'N', equed = 'N';
  ptrdiff_t n = this->r;
  ptrdiff_t nrhs = 0;

  knut_dgesvx(FACT, trans, n, nrhs, this->m, n, this->mf, n,
              ipiv, equed, nullptr, nullptr, nullptr, n, nullptr, n,
              &rcond, nullptr, nullptr, work, iwork, &info);

  // cout<<"work(1): "<<work[0]<<" rcond< "<<rcond<<"\n";
  if (info != 0) cout << "KNLuMatrix::Fact: Error code" << info << '\n';
  fact = true;
}

void KNLuMatrix::solve(KNVector& x, const KNVector& b, bool trans)
{
  P_ASSERT_X(this->r == this->c, "KNLuMatrix::Solve not a square matrix\n");
  P_ASSERT_X(b.size() == this->c, "KNLuMatrix::Solve KNVector sizes differ: M-V\n");
  P_ASSERT_X(b.size() == x.size(), "KNLuMatrix::Solve KNVector sizes differ: V-V\n");

  const KNMatrix& AA = *this;
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
      if (!fact) luFactorize();

      char FACT = 'F', equed = 'N', TRANS;
      if (trans) TRANS = 'T';
      else TRANS = 'N';
      size_t n = this->r;
      size_t nrhs = 1;

      knut_dgesvx(FACT, TRANS, n, nrhs, this->m, n, this->mf, n,
                  ipiv, equed, nullptr, nullptr, b.v, n, x.pointer(), n,
                  &rcond, &ferr, &berr, work, iwork, &info);
      if (info != 0) cout << "KNLuMatrix::Solve: Error code" << info << '\n';
      break;
  }
//  KNVector err(b);
//  timesXPlusY( err, x, b, 1.0, -1.0, trans );
//  cout<<"dim: "<<this->r<<" b: "<<sqrt(b*b)<<" e: "<<sqrt(err*err)<<"\n";
}

void KNLuMatrix::solve(KNMatrix& x, const KNMatrix& b, bool trans)
{
  P_ASSERT_X(this->r == this->c, "KNLuMatrix::Solve not a square matrix\n");
  P_ASSERT_X(b.col() != x.col(), "KNLuMatrix::Solve KNMatrix columns differ\n");
  P_ASSERT_X((b.row() == this->c) && (b.row() == x.row()), "KNLuMatrix::Solve KNMatrix rows differ\n");

  const KNMatrix& AA = *this;
  double det;
  switch (this->r)
  {
    case 1:
      for (size_t i = 0; i < b.col(); i++)
      {
        x(0, i) = b(0, i) / AA(0, 0);
      }
      break;
    case 2:
      det = 1.0 / (AA(0, 0) * AA(1, 1) - AA(0, 1) * AA(1, 0));
      if (!trans)
      {
        for (size_t i = 0; i < b.col(); i++)
        {
          x(0, i) = det * (AA(1, 1) * b(0, i) - AA(0, 1) * b(1, i));
          x(1, i) = det * (-AA(1, 0) * b(0, i) + AA(0, 0) * b(1, i));
        }
      }
      else
      {
        for (size_t i = 0; i < b.col(); i++)
        {
          x(0, i) = det * (AA(1, 1) * b(0, i) - AA(1, 0) * b(1, i));
          x(1, i) = det * (-AA(0, 1) * b(0, i) + AA(0, 0) * b(1, i));
        }
      }
      break;
    default:
      if (!fact) luFactorize();

      char FACT = 'F', equed = 'N', TRANS;
      if (trans) TRANS = 'T';
      else TRANS = 'N';
      double *fberr = new double[2*b.col()+1];
      size_t n = this->r;
      size_t nrhs = b.col();

      knut_dgesvx(FACT, TRANS, n, nrhs, this->m, n, this->mf, n,
                  ipiv, equed, nullptr, nullptr, b.m, n, x.m, n,
                  &rcond, fberr, fberr + b.col(), work, iwork, &info);

      if (info != 0) cout << "KNLuMatrix::Solve: Error code" << info << '\n';
      delete[] fberr;
      break;
  }
}
