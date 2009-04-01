// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2007 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include <cmath>
#include <cfloat>
#include "multipliers.h"
#include "matrix.h"
#include "pointtype.h"

static inline void findTrivialIndices(const Vector& mulRe, const Vector& mulIm,
                                      const int lp, const int pd, const int ns,
                                      int *imin, double* dmin)
{
  for (int i = 0; i < lp+pd+ns; ++i)
  {
    imin[i] = -1;
    dmin[i] = DBL_MAX;
  }
  for (int j = 0; j < lp+pd+ns; ++j)
  {
    for (int i = 0; i < mulRe.Size(); ++i)
    {
      bool skip = false;
      for (int k = 0; k < j; ++k)
      {
        if (imin[k] == i) skip = true;
      }
      if (skip) continue;
      else
      {
        double mabs;
        if (j < lp) mabs = fabs(sqrt((mulRe(i) - 1.0) * (mulRe(i) - 1.0) + mulIm(i) * mulIm(i)));
        else if (j < lp+pd) mabs = fabs(sqrt((mulRe(i) + 1.0) * (mulRe(i) + 1.0) + mulIm(i) * mulIm(i)));
        else mabs = fabs(sqrt(mulRe(i) * mulRe(i) + mulIm(i) * mulIm(i)) - 1.0);
        if ( mabs < dmin[j] )
        {
          dmin[j] = mabs;
          imin[j] = i;
        }
      }
    }
  }
}

int unstableMultipliers(const Vector& mulRe, const Vector& mulIm, const int lp, const int pd, const int ns)
{
#define NCRIT 6
  P_ERROR_X1(lp + pd + ns < NCRIT, "Too many critical multipliers. Change the value of NCRIT on the previous line.");
  int imin[NCRIT];
  double dmin[NCRIT];
  findTrivialIndices(mulRe, mulIm, lp, pd, ns, imin, dmin);
  int ustab = 0;
  for (int i = 0; i < mulRe.Size(); ++i)
  {
    const double mabs = (mulRe(i) * mulRe(i) + mulIm(i) * mulIm(i));
    bool ok = true;
    for (int j = 0; j < lp+pd+ns; ++j) if (i == imin[j]) ok = false;
    if (ok && mabs >= 1.0)  ++ustab;
  }
  return ustab;
#undef NCRIT
}

PtType bifurcationType(const Vector& mulRe, const Vector& mulIm, const int lp, const int pd, const int ns)
{
#define NCRIT 6
  const int aut = lp+pd+ns;
  P_ERROR_X1(aut < NCRIT, "Too many critical multipliers. Change the value of NCRIT on the previous line.");
  int imin[NCRIT];
  double dmin[NCRIT];
  findTrivialIndices(mulRe, mulIm, lp, pd, ns, imin, dmin);
  double dminLP = DBL_MAX, dminPD = DBL_MAX, dminNS = DBL_MAX;
  int iminLP = -1, iminPD = -1, iminNS = -1;
  for (int i = 0; i < mulRe.Size(); ++i)
  {
    const double mre = mulRe(i);
    const double mim = mulIm(i);
    bool ok = true;
    for (int j = 0; j < aut; ++j) if (i == imin[j]) ok = false;
    if (ok)
    {
      const double LPabs = fabs(sqrt((mre - 1.0) * (mre - 1.0) + mim * mim));
      const double PDabs = fabs(sqrt((mre + 1.0) * (mre + 1.0) + mim * mim));
      const double NSabs = fabs(sqrt(mre * mre + mim * mim) - 1.0);
      if ((dminLP > LPabs) && (mim == 0.0))
      {
        dminLP = LPabs;
        iminLP = i;
      }
      if ((dminPD > PDabs) && (mim == 0.0))
      {
        dminPD = PDabs;
        iminPD = i;
      }
      if ((dminNS > NSabs) && (mim != 0.0))
      {
        dminNS = NSabs;
        iminNS = i;
      }
    }
  }
  if ((dminLP < dminPD) && (dminLP < dminNS)) return BifTFLP;
  else if ((dminPD < dminLP) && (dminPD < dminNS)) return BifTFPD;
  else if ((dminNS < dminPD) && (dminNS < dminLP)) return BifTFNS;
  else return SolTF;
#undef NCRIT
}
