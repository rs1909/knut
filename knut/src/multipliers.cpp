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
#include <cstring>
#include <algorithm>
#include <vector>
#include <list>
#include "multipliers.h"
#include "matrix.h"
#include "pointtype.h"

// findTrivialIndices :
// 	Input :
//		mulRe : real part of the multipliers
//		mulIm : imaginary part of the multipliers
//		lp    : number of +1s to find
//		pd    : number of -1s to find
//		ns    : number of unit norms other than -1 or +1 to find
//	Output :
//		imin : array of indices of the found multipliers
//		dmin : the distance from the target (should be small, but not checked)
static inline size_t findTrivialIndices(const KNVector& mulRe, const KNVector& mulIm,
                                     const size_t lp, const size_t pd, const size_t ns,
                                     size_t *imin, double* dmin, size_t )
{
  for (size_t i = 0; i < lp+pd+ns; ++i)
  {
    imin[i] = mulRe.size();
    dmin[i] = DBL_MAX;
  }
  for (size_t j = 0; j < lp+pd+ns; ++j)
  {
    for (size_t i = 0; i < mulRe.size(); ++i)
    {
      bool skip = false;
      for (size_t k = 0; k < j; ++k)
      {
        if (imin[k] == i) skip = true;
      }
      if (skip) continue;
      else
      {
        double mabs = DBL_MAX;
        if (j < lp) { mabs = fabs(mulRe(i) - 1.0) + fabs(mulIm(i)); }
        else if (j < lp+pd) { mabs = fabs(mulRe(i) + 1.0) + fabs(mulIm(i)); }
        else { mabs = fabs(sqrt(mulRe(i) * mulRe(i) + mulIm(i) * mulIm(i)) - 1.0); }
        if ( mabs < dmin[j] && mabs < 0.1)
        {
          dmin[j] = mabs;
          imin[j] = i;
        }
      }
    }
  }
  size_t ustab = 0;
  for (size_t i = 0; i < mulRe.size(); ++i)
  {
    const double mabs = (mulRe(i) * mulRe(i) + mulIm(i) * mulIm(i));
    bool ok = true;
    for (size_t j = 0; j < lp+pd+ns; ++j) if (i == imin[j]) ok = false;
    if (ok && (mabs > 1.0))  ++ustab;
  }
//  std::cerr << pt << ". ";
//  for (int j = 0; j < lp+pd+ns; ++j) 
//  {
//    if (imin[j] < 0) std::cerr << "didn't find " << j << "\n";
//    std::cerr << "(" << mulRe(imin[j]) << "," << mulIm(imin[j]) << ") ";
//  }
//  std::cerr << " -U- ";
//  for (int i = 0; i < mulRe.size(); ++i)
//  {
//    const double mabs = (mulRe(i) * mulRe(i) + mulIm(i) * mulIm(i));
//    bool ok = true;
//    for (int j = 0; j < lp+pd+ns; ++j) if (i == imin[j]) ok = false;
//    if (ok && (mabs > 1.0))  { std::cerr << "(" << mulRe(i) << "," << mulIm(i) << ") ";}
//  }
//  std::cerr << "\n";
//  for (int j = lp+pd; j < lp+pd+ns; ++j)
//  {
//    if (mulIm(imin[j]) == 0) std::cerr << "findTrivialIndices : NS It is a real multiplier\n";
//  }
  return ustab;
}

size_t unstableMultipliers(const KNVector& mulRe, const KNVector& mulIm, const size_t lp, const size_t pd, const size_t ns, size_t pt)
{
#define NCRIT 6
  P_ERROR_X1(lp + pd + ns < NCRIT, "Too many critical multipliers. Change the value of NCRIT on the previous line.");
  size_t imin[NCRIT];
  double dmin[NCRIT];
  return findTrivialIndices(mulRe, mulIm, lp, pd, ns, imin, dmin, pt);
#undef NCRIT
}

// Re0,Im0 must have more unstable eigenvalues!
BifType bifurcationType(const KNVector& mr0, const KNVector& mi0, 
                        const KNVector& mr1, const KNVector& mi1,
                        const size_t lp, const size_t pd, const size_t ns, size_t pt0, size_t pt1)
{
#define NCRIT 6
  const size_t aut = lp+pd+ns;
  P_ERROR_X1(aut < NCRIT, "Too many critical multipliers. Change the value of NCRIT on the previous line.");
  size_t imin0[NCRIT];
  size_t imin1[NCRIT];
  size_t itmp[NCRIT];
  double dmin0[NCRIT];
  double dmin1[NCRIT];
  double dtmp[NCRIT];
  const size_t us0 = findTrivialIndices(mr0, mi0, lp, pd, ns, imin0, dmin0, pt0);
  const size_t us1 = findTrivialIndices(mr1, mi1, lp, pd, ns, imin1, dmin1, pt1);
  P_ERROR_X1(us0 != us1, "No bifurcation!");
  const KNVector* mulRe0;
  const KNVector* mulRe1;
  const KNVector* mulIm0;
  const KNVector* mulIm1;
  if (us0 > us1)
  {
    mulRe0 = &mr0; mulRe1 = &mr1;
    mulIm0 = &mi0; mulIm1 = &mi1;
  } else
  {
    mulRe1 = &mr0; mulRe0 = &mr1;
    mulIm1 = &mi0; mulIm0 = &mi1;
    memcpy(itmp, imin0, sizeof(itmp)); memcpy(imin0, imin1, sizeof(itmp)); memcpy(imin1, itmp, sizeof(itmp));
    memcpy(dtmp, dmin0, sizeof(dtmp)); memcpy(dmin0, dmin1, sizeof(dtmp)); memcpy(dmin1, dtmp, sizeof(dtmp));
  }
  // the maximum eigenvalue
  size_t imaxa = mulRe0->size();
  double dmaxa = 0;
  // the closest eigenvalue to "0" in "1"
  std::vector<size_t> closest(mulRe0->size());
  for (size_t i = 0; i < mulRe0->size(); ++i)
  {
    const double mre0 = (*mulRe0)(i);
    const double mim0 = (*mulIm0)(i);
    bool ok0 = true;
    for (size_t j = 0; j < aut; ++j) if (i == imin0[j]) ok0 = false;
    if (ok0 && mre0*mre0+mim0*mim0 > 1.0)
    {
      if ( dmaxa < mre0*mre0+mim0*mim0) { dmaxa = mre0*mre0+mim0*mim0; imaxa = i; }
      // second loop : find the closest multiplier to mulRe0(i), mulIm0(i)
      size_t imina = mulRe1->size();
      double dmina = DBL_MAX;
      for (size_t k = 0; k < mulRe1->size(); ++k)
      {
        const double mre1 = (*mulRe1)(k);
        const double mim1 = (*mulIm1)(k);
        bool ok1 = true;
        for (size_t j = 0; j < aut; ++j) if (k == imin1[j]) ok1 = false;
        if (ok1 && mre1*mre1+mim1*mim1 > 1.0)
        {
          double mabsa = (mre0-mre1)*(mre0-mre1)+(mim0-mim1)*(mim0-mim1);
          if (dmina > mabsa) { dmina = mabsa; imina = k; }
        }
      }
      closest[i] = imina;
    } else closest[i] = mulRe1->size();
  }
  // find the duplicate in closest
  std::vector<std::list<size_t> > duplicate(mulRe0->size());
  std::vector<bool> valid(mulRe0->size(), true);
  for (size_t i=0; i<closest.size(); ++i)
  {
    duplicate[i] = std::list<size_t>(0);
    for (size_t j=i+1; j<closest.size(); ++j)
    {  
      if (valid[i]&&(closest[i] < mulRe1->size())&&(closest[i] == closest[j]))
      {
        duplicate[i].push_back(j);
        valid[j] = false;
      }
    }
    if (!duplicate[i].empty()) duplicate[i].push_back(i);
  }
  size_t dnum = 0;
  size_t imax[3] = {mulRe0->size(), mulRe0->size(), mulRe0->size()};
  double dmax[3] = {0, 0, 0};
  for (size_t i=0; i<duplicate.size(); ++i)
  {
    if (!duplicate[i].empty()) 
    {
      const double mre0 = (*mulRe1)(closest[i]);
      const double mim0 = (*mulIm1)(closest[i]);
//      std :: cerr << "DNUM " << dnum << " (" << mre0 << "," << mim0 << ") ==>> ";
      for (std::list<size_t>::iterator it = duplicate[i].begin(); it != duplicate[i].end(); it++)
      {
        const double mre1 = (*mulRe0)(*it);
        const double mim1 = (*mulIm0)(*it);
//        std :: cerr << " (" << mre1 << "," << mim1 << ")";
        double mabs = (mre0-mre1)*(mre0-mre1)+(mim0-mim1)*(mim0-mim1);
        if (dmax[dnum] < mabs) { dmax[dnum] = mabs; imax[dnum] = *it; }
      }
//      std :: cerr << "\n";
      dnum++;
    }
    if (dnum > 3) break;
  }
//  std :: cerr << "DNUM " << dnum << "\n";
  if (dnum == 0)
  {
    // nothing to do
//    std :: cerr << "Bif Mul 0 = (" << (*mulRe0)(imaxa) << "," << (*mulIm0)(imaxa) << ")\n";
  } else if (dnum == 1)
  {
//    std :: cerr << "Bif Mul 1 = (" << (*mulRe0)(imax[0]) << "," << (*mulIm0)(imax[0]) << ")\n";
    imaxa = imax[0];
  } else if (dnum == 2)
  {
//    std :: cerr << "Bif Mul 2 = (" << (*mulRe0)(imax[0]) << "," << (*mulIm0)(imax[0]) << ") - (" << (*mulRe0)(imax[1]) << "," << (*mulIm0)(imax[1]) << ")\n";
//    if (imax[0] == imax[1]) std :: cerr << "The same index\n";
    // this must be a complex pair drawn to another complex pair
    if (((*mulRe0)(imax[0]) == (*mulRe0)(imax[1]))&&((*mulIm0)(imax[0]) == -(*mulIm0)(imax[1])))
    {
      imaxa = imax[0];
    } else return BifUN;
  } else return BifUN;
  
  const double mrea = (*mulRe0)(imaxa);
  const double mima = (*mulIm0)(imaxa);
  // find out where did the eigenvalue come from ?
  size_t iminb = mulRe1->size();
  double dminb = DBL_MAX;
  for (size_t k = 0; k < mulRe1->size(); ++k)
  {
    const double mre1 = (*mulRe1)(k);
    const double mim1 = (*mulIm1)(k);
    bool ok1 = true;
    for (size_t j = 0; j < aut; ++j) if (k == imin1[j]) ok1 = false;
    if (ok1 && mre1*mre1+mim1*mim1 <= 1.0)
    {
      double mabsb = (mrea-mre1)*(mrea-mre1)+(mima-mim1)*(mima-mim1);
      if (dminb > mabsb) { dminb = mabsb; iminb = k; }
    }
  }
  const double LPabs = fabs(sqrt((mrea - 1.0) * (mrea - 1.0) + mima * mima));
  const double PDabs = fabs(sqrt((mrea + 1.0) * (mrea + 1.0) + mima * mima));
  const double NSabs = fabs(sqrt(mrea * mrea + mima * mima) - 1.0);
  if (iminb < mulRe1->size())
  {
    const double mreb = (*mulRe1)(iminb);
    const double mimb = (*mulIm1)(iminb);

  //  std :: cerr << "From = (" << mreb << "," << mimb << ")\n";
    if ( (mimb != 0.0)&&(mima != 0.0) )
    {
      return BifNS;
    } else if ((mimb == 0.0)&&(mima == 0.0))
    {
      if (LPabs < PDabs) return BifLP;
      else return BifPD;
    } else return BifUN;
  } else
  {
    if (mima != 0.0) return BifNS;
    else
    {
      if (LPabs < PDabs) return BifLP;
      else return BifPD;
    }
  }
#undef NCRIT
}
