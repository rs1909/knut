// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "config.h"

#include "testfunct.h"
#include "matrix.h"
#include "spmatrix.h"
#include "hypermatrix.h"
#include "ncolloc.h"
#include "pointtype.h"

#ifdef DEBUG
#include <iomanip>
#include <fstream>
#endif

#define NDIM (col.nDim())
#define NDEG (col.nDeg())
#define NINT (col.nInt())
#define NTAU (col.nTau())
#define NPAR (col.nPar())

template<bool trans> inline void rotbord(KNVector& V, KNDdeBvpCollocation& col, const KNVector& IN, KNArray1D<size_t>& Re, KNArray1D<size_t>& Im)
{
  V.clear();
  for (size_t idx = 0; idx < NINT*NDEG + 1; idx++)
  {
    for (size_t k = 0; k < Re.size(); k++)
    {
      if (trans)
      {
        V(Re(k) + NDIM*idx) = - IN(Im(k) + NDIM * idx);
        V(Im(k) + NDIM*idx) = IN(Re(k) + NDIM * idx);
      }
      else
      {
        V(Re(k) + NDIM*idx) = IN(Im(k) + NDIM * idx);
        V(Im(k) + NDIM*idx) = - IN(Re(k) + NDIM * idx);
      }
    }
  }
}

inline void conjugate(KNVector& out, const KNVector& inp)
{
  for (size_t i = 0; i < out.size() / 2; ++i)
  {
    out(2*i)   = -inp(2 * i + 1);
    out(2*i + 1) =  inp(2 * i);
  }
}

KNTestFunctional::KNTestFunctional(KNDdeBvpCollocation& col, double Z) :
    ZZ(Z),
    AHAT(NDIM*(NDEG*NINT + 1), 0, 1, NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    A_p(NDIM*(NDEG*NINT + 1)),
    A_x('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    rhs(NDIM*(NDEG*NINT + 1)),
    uu(NDIM*(NDEG*NINT + 1)), vv(NDIM*(NDEG*NINT + 1)),
    uudiff(NDIM*(NDEG*NINT + 1)), vvdiff(NDIM*(NDEG*NINT + 1)),
    vvData(NDIM, 2*NTAU + 1, NDEG*NINT)
{
  first = true;
}

KNTestFunctional::~KNTestFunctional()
{}

double KNTestFunctional::initStep()
{
  double one = 1.0;
  double gg = 0.0;
  double hh = 0.0;
  double ggdiff, hhdiff;

  AHAT.multiply<false>(rhs, one, vv, gg);
  one -= 1.0;
  AHAT.solve(vvdiff, ggdiff, rhs, one);
  AHAT.multiply<true>(rhs, one, uu, hh);
  one -= 1.0;
  AHAT.solveTr(uudiff, hhdiff, rhs, one);
  vv -= vvdiff;
  gg -= ggdiff;
  uu -= uudiff;
  hh -= hhdiff;
  vv /= sqrt(vv * vv);
  uu /= sqrt(uu * uu);
  AHAT.getA31(0) = vv;
  AHAT.getA13(0) = uu;
  const double diffnorm = std::max<double>(sqrt(vvdiff * vvdiff + ggdiff * ggdiff), sqrt(uudiff * uudiff + hhdiff * hhdiff));
#ifdef DEBUG
  std::cout << "dnor " << diffnorm << "\n";
#endif
  return diffnorm;
}

void KNTestFunctional::init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  // creating the matrix
  col.jotf_x(AHAT.getA11(), par, ZZ);
  AHAT.getA13(0).random();
  AHAT.getA31(0).random();
  AHAT.getA33()(0, 0) = 0.0;
  // norming the borders
  AHAT.getA31(0) /= sqrt(AHAT.getA31(0) * AHAT.getA31(0));
  AHAT.getA13(0) /= sqrt(AHAT.getA13(0) * AHAT.getA13(0));
  vv = AHAT.getA31(0);
  uu = AHAT.getA13(0);

  double diffnorm = 1.0;
  size_t it = 0;
  do
  {
    diffnorm = initStep();
  }
  while ((++it < kernIter) && (diffnorm > kernEps));
  if (diffnorm > kernEps) std::cout << "KNTestFunctional::Init: warning: No convergence in finding the singular vector. Residual = " << diffnorm << "\n";
//  std::cout<<"TF: "<<gg<<", "<<hh<<"\n";
#ifdef DEBUG
  std::ofstream uufile("eigv");
  std::ofstream vvfile("eigvstar");
  uufile << std::scientific;
  vvfile << std::scientific;
  uufile.precision(12);
  vvfile.precision(12);

  for (size_t i = 0; i < NINT; i++)
  {
    for (size_t j = 0; j < NDEG + 1; j++)
    {
      const double t = col.Profile(i, j);
      for (size_t p = 0; p < NDIM; p++)
      {
        vvfile << vv(p + (j + i*NDEG)*NDIM) << "\t";
        uufile << uu(p + (j + i*NDEG)*NDIM) << "\t";
      }
      vvfile << par(0)*t << "\n";
      uufile << par(0)*t << "\n";
    }
  }
#endif
}

double KNTestFunctional::funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol)
{
  double one = 1.0;
  double gg = 0.0;
  double hh = 0.0;
  double ggdiff, hhdiff;

  if (first)
  {
    init(col, par, sol);
    first = false;
  }

  rhs.clear();
  col.jotf_x(AHAT.getA11(), par, ZZ);
  AHAT.solve(vv, gg, rhs, one);
  AHAT.solveTr(uu, hh, rhs, one);
  // vvdiff is temporary ...
//  std::cout<<"TF: "<<gg<<", "<<hh;
//  if( gg > 0.0 ) std::cout<<"\t+++\n";
//  else           std::cout<<"\t---\n";
  return gg;
}

double KNTestFunctional::funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/, size_t alpha)
{
  col.interpolate(vvData, vv);
  col.jotf_x_p(A_p, par, vvData, ZZ, alpha);
  const double gg_p = (uu * A_p);
//  std::cout<<"GP "<<alpha<<":"<<gg_p;
  return gg_p; // positive, because the test function is not negated in the Jacobian
}

void   KNTestFunctional::funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  col.interpolate(vvData, vv);
  col.jotf_x_x(A_x, par, vvData, ZZ);
  func = !A_x * uu; // positive, because the test function is not negated in the Jacobian
//  func *= -1;
}

void   KNTestFunctional::kernel(KNVector& phi)
{
  phi = vv;
}

/// ---------------------------------------------------------
/// test function for TORUS BIFURCATIONS
/// ---------------------------------------------------------

KNComplexTestFunctional::KNComplexTestFunctional(KNDdeBvpCollocation& col) :
    first(true),
    AHAT(2*NDIM*(NDEG*NINT + 1), 0, 2, 4*NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    A_p(2*NDIM*(NDEG*NINT + 1)),
    A_x('R', 2*NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1), 2*NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    rhs(2*NDIM*(NDEG*NINT + 1)),
    one((size_t)2),
    uu(2*NDIM*(NDEG*NINT + 1)),
    vv(2*NDIM*(NDEG*NINT + 1)),
    uudiff(2*NDIM*(NDEG*NINT + 1)),
    vvdiff(2*NDIM*(NDEG*NINT + 1)),
    gg((size_t)2),
    hh((size_t)2),
    ggdiff((size_t)2),
    hhdiff((size_t)2),
    vvDataRe(NDIM, NTAU + 1, NDEG*NINT),
    vvDataIm(NDIM, NTAU + 1, NDEG*NINT)
{}

KNComplexTestFunctional::~KNComplexTestFunctional()
{}

double KNComplexTestFunctional::initStep()
{
  AHAT.multiply<false>(2, rhs, one, vv, gg);
  one(0) -= 1.0;
  AHAT.solve(2, vvdiff, ggdiff, rhs, one);
  AHAT.multiply<true>(2, rhs, one, uu, hh);
  one(0) -= 1.0;
  AHAT.solveTr(2, uudiff, hhdiff, rhs, one);
  vv -= vvdiff;
  gg -= ggdiff;
  uu -= uudiff;
  hh -= hhdiff;
  vv /= sqrt(vv * vv);
  uu /= sqrt(uu * uu);
  AHAT.getA31(0) = vv;
  AHAT.getA13(0) = uu;
  conjugate(AHAT.getA13(1), AHAT.getA13(0));
  conjugate(AHAT.getA31(1), AHAT.getA31(0));
  const double diffnorm = std::max<double>(sqrt(uudiff * uudiff), sqrt(vvdiff * vvdiff));
#ifdef DEBUG
  std::cout << "dnorCX " << diffnorm << "\n";
#endif
  return diffnorm;
}

void KNComplexTestFunctional::init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/,
                         double Re, double Im)
{
  ZRe = Re;
  ZIm = Im;
  col.jotf_x(AHAT.getA11(), par, Re, Im);
  AHAT.getA13(0).random();
  AHAT.getA31(0).random();
  AHAT.getA33().clear();
  // norming the borders
  AHAT.getA31(0) /= sqrt(AHAT.getA31(0) * AHAT.getA31(0));
  AHAT.getA13(0) /= sqrt(AHAT.getA13(0) * AHAT.getA13(0));
  // conjugate
  conjugate(AHAT.getA31(1), AHAT.getA31(0));
  conjugate(AHAT.getA13(1), AHAT.getA13(0));

  double diffnorm = 1.0;
  size_t it = 0;
  do
  {
    diffnorm = initStep();
  }
  while ((++it < kernIter) && (diffnorm > kernEps));
  if (diffnorm > kernEps) std::cout << "KNComplexTestFunctional::Init: warning: No convergence in finding the singular vector. Residual = " << diffnorm << "\n";
//  std::cout<<"TF: "<<gg<<", "<<hh<<"\n";
}

void KNComplexTestFunctional::funct(double& f1, double& f2,
                          KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, double Re, double Im)
{
  if (first)
  {
    init(col, par, sol, Re, Im);
    first = false;
  }

  ZRe = Re;
  ZIm = Im;
  one(0) = 1.0; one(1) = 0.0;
  rhs.clear();
  col.jotf_x(AHAT.getA11(), par, Re, Im);
  AHAT.solve(2, vv, gg, rhs, one);
  AHAT.solveTr(2, uu, hh, rhs, one);
  // for later use
  col.interpolateComplex(vvDataRe, vvDataIm, vv);

  f1 = gg(0);
  f2 = gg(1);
//  std::cout<<"gg: "<<gg(0)<<", "<<gg(1)<<", "<<hh(0)<<", "<<hh(1)<<"\n";
//  if( gg(0) > 0.0 ) std::cout<<"\t+++\n";
//  else               std::cout<<"\t---\n";
//  if( gg(1) > 0.0 ) std::cout<<"\t+++\n";
//  else               std::cout<<"\t---\n";
}

void KNComplexTestFunctional::funct_p(double& f1, double& f2,
                            KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/,
                            size_t alpha)
{
  col.jotf_x_p(A_p, par, vvDataRe, vvDataIm, ZRe, ZIm, alpha);
  conjugate(uudiff, uu);
  f1 = (uu * A_p);
  f2 = (uudiff/*uu^conj*/ * A_p);
}

void KNComplexTestFunctional::funct_z(double& f1, double& f2,
                            KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  col.jotf_x_z(A_p, par, AHAT.getA31(0), vvDataRe, vvDataIm, ZRe, ZIm);
  
  conjugate(uudiff, uu);
  const double dzre = (uu * A_p);
  const double dzim = (uudiff/*uu^conj*/ * A_p);
  f1 = (-dzim * ZRe - dzre * ZIm);
  f2 = (dzre * ZRe - dzim * ZIm);
}

void KNComplexTestFunctional::funct_x(KNVector& func1, KNVector& func2,
                            KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  col.jotf_x_x(A_x, par, vvDataRe, vvDataIm, ZRe, ZIm);
  conjugate(uudiff, uu);
  func1 = !A_x * uu; /*uu*/
  func2 = !A_x * uudiff; /*uu^conj*/
}

void KNComplexTestFunctional::kernel(KNVector& Re, KNVector& Im, double& alpha)
{
  P_ASSERT_X((2*Re.size() == vv.size()) && (2*Im.size() == vv.size()), "KNComplexTestFunctional::Switch: Bad sizes\n");
  std::cout << "zRe=" << ZRe << ", zIm=" << ZIm << "\n";
  for (size_t i = 0; i < Re.size(); i++)
  {
    Re(i) = vv(2 * i);
    Im(i) = vv(2 * i + 1);
  }
  if (ZRe > 0.0)
  {
    alpha = atan(fabs(ZIm / ZRe));
  }
  else
  {
    alpha = atan(fabs(ZRe / ZIm)) + M_PI / 2.0;
  }
}

double KNComplexTestFunctional::kernelComplex(double& newperiod, KNVector& Re, KNVector& Im, KNDdeBvpCollocation& col, const KNVector& par)
{
  double norm1 = 0.0;
  for (size_t i = 0; i < NDIM; i++)
  {
    Re(i) = vv(2 * i);
    Im(i) = vv(2 * i + 1);
    if (norm1 < fabs(Re(i))) norm1 = fabs(Re(i));
  }

#ifdef DEBUG
  std::cout << "Printing: eigenvec\n";
  std::ofstream file("eigenvec");
  file << std::scientific;
  file.precision(12);
#endif
  // computing the period of the new branch
  double period = par(0);
  double angle = par(VarToIndex(VarAngle,NPAR));
  if (par(VarToIndex(VarAngle,NPAR)) < 0) angle += 2*M_PI;
  
  double dist;
  size_t it = 0;
  do {
  	dist = 0.0;
    period = par(0)*(2*M_PI/angle);
  
    for (size_t i = 0; i < NINT; i++)
    {
      for (size_t j = 0; j < NDEG + 1; j++)
      {
        const double t = col.Profile(i, j);
        for (size_t p = 0; p < NDIM; p++)
        {
          const double t2 = t*par(0)/period;
          const double d1 =(cos(2.0 * M_PI * t2) * vv(2 * p) + sin(2.0 * M_PI * t2) * vv(2 * p + 1));
          const double d2 = vv(2*(p + (j + i*NDEG)*NDIM));
          const double d = d2 - d1;
          if (dist < fabs(d)) dist = fabs(d);
        #ifdef DEBUG
          file << d1 << "\t" << d2 << "\t";
        #endif
        }
      #ifdef DEBUG
        file << par(0)*t << "\n";
      #endif
      }
    }
    angle += 2*M_PI;
    ++it;
  } while ( (dist/norm1 > 0.001) && (it < 12) );
  
  newperiod = period;
  return dist/norm1;
}

/// ---------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems
/// ---------------------------------------------------------

KNLpAutTestFunctional::KNLpAutTestFunctional(KNDdeBvpCollocation& col, double Z) :
    ZZ(Z),
    AHAT(NDIM*(NDEG*NINT + 1), 0, 2, NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    A_p(NDIM*(NDEG*NINT + 1)),
    A_x('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    mB('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    mB_p(NDIM*(NDEG*NINT + 1)),
    mB_x('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    rhs(NDIM*(NDEG*NINT + 1)),
    phi(NDIM*(NDEG*NINT + 1)),
    temp(NDIM*(NDEG*NINT + 1)),
    DpPhi(NDIM*(NDEG*NINT + 1)),
    DxPhi(NDIM*(NDEG*NINT + 1)),
    uu2(NDIM*(NDEG*NINT + 1)),
    vv2(NDIM*(NDEG*NINT + 1)),
    gg2((size_t)2),
    hh2((size_t)2),
    one2((size_t)2),
    phiData(NDIM, 2*NTAU + 1, NDEG*NINT),
    vv2Data(NDIM, 2*NTAU + 1, NDEG*NINT),
    solMSHData(NDIM, 2*NTAU + 1, NDEG*NINT + 1)
{
  first = true;
}

KNLpAutTestFunctional::~KNLpAutTestFunctional()
{}

double KNLpAutTestFunctional::initStep()
{
  return 0.0;
}

void KNLpAutTestFunctional::init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol)
{
  // creating the matrix
  col.jotf_x(AHAT.getA11(), par, ZZ);
  col.jotf_mB(mB, par, ZZ);

  col.interpolateOnMesh(solMSHData, sol);
  col.jotf_trivialKernelOnMesh(phi, par, solMSHData);

  col.star(AHAT.getA31(0), phi);
  AHAT.getA13(0) = mB * phi;

  AHAT.getA13(1).random();
  AHAT.getA31(1).random();
  AHAT.getA33().clear(); // 2x2 matrix
  // norming the borders
  AHAT.getA31(1) /= sqrt(AHAT.getA31(1) * AHAT.getA31(1));
  AHAT.getA13(1) /= sqrt(AHAT.getA13(1) * AHAT.getA13(1));

  // generating the non-trivial kernel
  one2(0) = 0.0;
  one2(1) = 1.0;

  KNVector v2diff(vv2), u2diff(uu2);
  double g2diff, h2diff;
  double diffnorm = 1.0;
  size_t it = 0;
  do
  {
    AHAT.solve(2, vv2, gg2, rhs, one2);
    AHAT.solveTr(2, uu2, hh2, rhs, one2);
    u2diff = AHAT.getA13(1);
    u2diff -= uu2;
    h2diff = AHAT.getA33(0, 1) - hh2(0);
    v2diff = AHAT.getA31(1);
    v2diff -= vv2;
    g2diff = AHAT.getA33(1, 0) - gg2(0);
    diffnorm = std::max<double>(sqrt(u2diff * u2diff + h2diff * h2diff), sqrt(v2diff * v2diff + g2diff * g2diff));
//   std::cout<<"std::max(u2diff, v2diff) "<<diffnorm<<"\n";
    const double nrm_v = (1.0 / sqrt(vv2 * vv2 + gg2(0) * gg2(0)));
    const double nrm_u = (1.0 / sqrt(uu2 * uu2 + hh2(0) * hh2(0)));
    AHAT.getA31(1) = nrm_v * vv2;
    AHAT.getA13(1) = nrm_u * uu2;
    AHAT.getA33(1, 0) = nrm_v * gg2(0);
    AHAT.getA33(0, 1) = nrm_v * hh2(0);
  }
  while ((++it < kernIter) && (diffnorm > kernEps));
  if (diffnorm > kernEps) std::cout << "KNLpAutTestFunctional::Init: warning: No convergence in finding the singular vector. Residual = " << diffnorm << "\n";
}

double KNLpAutTestFunctional::funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol)
{
  if (first)
  {
    init(col, par, sol);
    first = false;
  }

  col.jotf_x(AHAT.getA11(), par, ZZ);
  col.jotf_mB(mB, par, ZZ);

  col.interpolateOnMesh(solMSHData, sol);
  col.jotf_trivialKernelOnMesh(phi, par, solMSHData);

  col.star(AHAT.getA31(0), phi);
  AHAT.getA13(0) = mB * phi;

  // continuing the non-trivial kernel
  AHAT.solve(2, vv2, gg2, rhs, one2);
  AHAT.solveTr(2, uu2, hh2, rhs, one2);

  const double nrm_v = (1.0 / sqrt(vv2 * vv2 + gg2(0) * gg2(0)));
  const double nrm_u = (1.0 / sqrt(uu2 * uu2 + hh2(0) * hh2(0)));
  AHAT.getA31(1) = nrm_v * vv2;
  AHAT.getA13(1) = nrm_u * uu2;
  AHAT.getA33(1, 0) = nrm_v * gg2(0);
  AHAT.getA33(0, 1) = nrm_u * hh2(0);

  // for subsequent use
  col.interpolate(phiData, phi);
  col.interpolate(vv2Data, vv2);

//  std::cout<<"TF: "<<gg2(1)<<", "<<hh2(1)<<":"<<gg2(0)<<", "<<hh2(0)<<"\n";
//  if( gg2(1) > 0.0 ) std::cout<<"\t+++\n";
//  else               std::cout<<"\t---\n";
  return gg2(1);
}

double KNLpAutTestFunctional::funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/, size_t alpha)
{
  col.jotf_x_p(A_p, par, vv2Data, ZZ, alpha);
  col.jotf_mB_p(mB_p, par, phiData, ZZ, alpha);

  col.jotf_trivialKernelOnMesh_p(DpPhi, par, solMSHData, alpha);
  // check
//  col.jotf_x_p( temp, par, phiData, ZZ, alpha );
//  std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
//  temp = AHAT.getA11() * DpPhi;
//  std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
//  std::cout<<"d: "<<alpha<<", "<<DpPhi*DpPhi<<"\n";
  // end check
  temp = mB * DpPhi;
  temp += mB_p;

  const double gg_p = (uu2 * A_p) + gg2(0) * (uu2 * temp) + hh2(0) * col.integrate(DpPhi, vv2);
//  std::cout<<"GP "<<alpha<<":"<<gg_p;
  return gg_p;
}

void   KNLpAutTestFunctional::funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  col.jotf_x_x(A_x, par, vv2Data, ZZ);
  col.jotf_mB_x(mB_x, par, phiData, ZZ);
  func = !A_x * uu2;

  temp = !mB * uu2;
  col.jotf_trivialKernelOnMesh_x<true>(DxPhi, par, solMSHData, temp);
  func += gg2(0) * DxPhi;
  temp = !mB_x * uu2;
  func += gg2(0) * temp;

  col.jotf_trivialKernelOnMesh_x<true>(DxPhi, par, solMSHData, vv2);
  col.star(temp, DxPhi);
  func += hh2(0) * temp;
}

void   KNLpAutTestFunctional::kernel(KNVector& phi)
{
  phi = vv2;
}

/// -----------------------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems with SIMMETRY
/// -----------------------------------------------------------------------

KNLpAutRotTestFunctional::KNLpAutRotTestFunctional(KNDdeBvpCollocation& col, KNArray1D<size_t> CRe, KNArray1D<size_t> CIm, double Z) :
    ZZ(Z),
    Re(CRe), Im(CIm),
    AHAT(NDIM*(NDEG*NINT + 1), 0, 3, NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    A_p(NDIM*(NDEG*NINT + 1)),
    A_x('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    mB('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    mB_p(NDIM*(NDEG*NINT + 1)),
    mB_x('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    phi(NDIM*(NDEG*NINT + 1)),
    DpPhi(NDIM*(NDEG*NINT + 1)),
    DxPhi(NDIM*(NDEG*NINT + 1)),
    LAM(NDIM*(NDEG*NINT + 1)),
    DxLAM(NDIM*(NDEG*NINT + 1)),
    uu3(NDIM*(NDEG*NINT + 1)),
    vv3(NDIM*(NDEG*NINT + 1)),
    gg3((size_t)3),
    hh3((size_t)3),
    rhs3(NDIM*(NDEG*NINT + 1)),
    one3((size_t)3),
    temp(NDIM*(NDEG*NINT + 1)),
    phiData(NDIM, 2*NTAU + 1, NDEG*NINT),
    LAMData(NDIM, 2*NTAU + 1, NDEG*NINT),
    vv3Data(NDIM, 2*NTAU + 1, NDEG*NINT),
    solMSHData(NDIM, 2*NTAU, NDEG*NINT + 1)
{
  first = true;
}

KNLpAutRotTestFunctional::~KNLpAutRotTestFunctional()
{}

double KNLpAutRotTestFunctional::initStep()
{
  return 0.0;
}

void KNLpAutRotTestFunctional::init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol)
{
  // creating the matrix
  col.jotf_x(AHAT.getA11(), par, ZZ);
  col.jotf_mB(mB, par, ZZ);

  col.interpolateOnMesh(solMSHData, sol);
  col.jotf_trivialKernelOnMesh(phi, par, solMSHData);

  temp = AHAT.getA11() * phi;
//   std::cout << "temp: " << temp*temp << " phi: " << phi*phi << "\n";
  col.star(AHAT.getA31(0), phi);
  AHAT.getA13(0) = mB * phi;

  rotbord<false>(LAM, col, sol, Re, Im);
  col.star(AHAT.getA31(1), LAM);
  AHAT.getA13(1) = mB * LAM;

  AHAT.getA13(2).random();
  AHAT.getA31(2).random();
  // norming the borders
  AHAT.getA31(2) /= sqrt(AHAT.getA31(2) * AHAT.getA31(2));
  AHAT.getA13(2) /= sqrt(AHAT.getA13(2) * AHAT.getA13(2));
  AHAT.getA33().clear(); // 2x2 matrix

  // generating the non-trivial kernel
  one3(0) = 0.0;
  one3(1) = 0.0;
  one3(2) = 1.0;
  rhs3.clear();

  KNVector v3diff(vv3), u3diff(uu3);
  double g30diff, h30diff, g31diff, h31diff;
  double diffnorm = 1.0;
  size_t it = 0;
  do
  {
    AHAT.solve(3, vv3, gg3, rhs3, one3);
    AHAT.solveTr(3, uu3, hh3, rhs3, one3);
    u3diff = AHAT.getA13(2);
    u3diff -= uu3;
    h30diff = AHAT.getA33(0, 2) - hh3(0);
    h31diff = AHAT.getA33(1, 2) - hh3(1);
    v3diff = AHAT.getA31(2);
    v3diff -= vv3;
    g30diff = AHAT.getA33(2, 0) - gg3(0);
    g31diff = AHAT.getA33(2, 1) - gg3(1);
    diffnorm = std::max<double>(sqrt(u3diff * u3diff + h30diff * h30diff + h31diff * h31diff),
                                sqrt(v3diff * v3diff + g30diff * g30diff + g31diff * g31diff));
//     std::cout << "std::max(u3diff, v3diff) " << diffnorm << "\n";
    const double nrm_v = (1.0 / sqrt(vv3 * vv3 + gg3(0) * gg3(0) + gg3(1) * gg3(1)));
    const double nrm_u = (1.0 / sqrt(uu3 * uu3 + hh3(0) * hh3(0) + hh3(1) * hh3(1)));
    AHAT.getA31(2)   = nrm_v * vv3;
    AHAT.getA33(2, 0) = nrm_v * gg3(0);
    AHAT.getA33(2, 1) = nrm_v * gg3(1);
    AHAT.getA13(2)   = nrm_u * uu3;
    AHAT.getA33(0, 2) = nrm_u * hh3(0);
    AHAT.getA33(1, 2) = nrm_u * hh3(1);
  }
  while ((++it < kernIter) && (diffnorm > kernEps));
  if (diffnorm > kernEps) std::cout << "KNLpAutRotTestFunctional::Init: warning: No convergence in finding the singular vector. Residual = " << diffnorm << "\n";
//   std::cout << "TF: " << gg3(2) << ", " << hh3(2) << "\n";
}

double KNLpAutRotTestFunctional::funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol)
{
  if (first)
  {
    init(col, par, sol);
    first = false;
  }

  col.jotf_x(AHAT.getA11(), par, ZZ);
  col.jotf_mB(mB, par, ZZ);

  col.interpolateOnMesh(solMSHData, sol);
  col.jotf_trivialKernelOnMesh(phi, par, solMSHData);

//  temp = AHAT.getA11() * phi;
//  std::cout<<"temp: "<<temp*temp<<" phi: "<<phi*phi<<"\n";

  col.star(AHAT.getA31(0), phi);
  AHAT.getA13(0) = mB * phi;

  rotbord<false>(LAM, col, sol, Re, Im);
  col.star(AHAT.getA31(1), LAM);
  AHAT.getA13(1) = mB * LAM;

//  temp = AHAT.getA11() * LAM;
//  std::cout<<"temp: "<<temp*temp<<" LAM: "<<LAM*LAM<<"\n";

  AHAT.solve(3, vv3, gg3, rhs3, one3);
  AHAT.solveTr(3, uu3, hh3, rhs3, one3);
  const double nrm_v = (1.0 / sqrt(vv3 * vv3 + gg3(0) * gg3(0) + gg3(1) * gg3(1)));
  const double nrm_u = (1.0 / sqrt(uu3 * uu3 + hh3(0) * hh3(0) + hh3(1) * hh3(1)));
  AHAT.getA31(2)   = nrm_v * vv3;
  AHAT.getA33(2, 0) = nrm_v * gg3(0);
  AHAT.getA33(2, 1) = nrm_v * gg3(1);
  AHAT.getA13(2)   = nrm_u * uu3;
  AHAT.getA33(0, 2) = nrm_u * hh3(0);
  AHAT.getA33(1, 2) = nrm_u * hh3(1);

  // for subsequent use
  col.interpolate(phiData, phi);
  col.interpolate(LAMData, LAM);
  col.interpolate(vv3Data, vv3);
 std::cout<<" TF1: "<<gg3(0)<<", "<<hh3(0)<<" TF2: "<<gg3(1)<<", "<<hh3(1)<<" TF3: "<<gg3(2)<<", "<<hh3(2)<<"\n";

//   if (gg3(2) > 0.0) std::cout << "\t+++\n";
//   else               std::cout << "\t---\n";
  return gg3(2);
}

double KNLpAutRotTestFunctional::funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/, size_t alpha)
{
  col.jotf_x_p(A_p, par, vv3Data, ZZ, alpha);
  col.jotf_trivialKernelOnMesh_p(DpPhi, par, solMSHData, alpha);
  double gg_p = (uu3 * A_p);
  // check
//  col.jotf_x_p( temp, par, phiData, ZZ, alpha );
//  std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
//  temp = AHAT.getA11() * DpPhi;
//  std::cout<<"t: "<<alpha<<", "<<temp*temp<<"\n";
//  std::cout<<"d: "<<alpha<<", "<<DpPhi*DpPhi<<"\n";
  // end check
  // first phi
  temp = mB * DpPhi;
  col.jotf_mB_p(mB_p, par, phiData, ZZ, alpha);
  temp += mB_p;
  gg_p += gg3(0) * (uu3 * temp) + hh3(0) * col.integrate(DpPhi, vv3);
  // second LAM ( mB * LAM -> only mB_p is nonzero )
  col.jotf_mB_p(mB_p, par, LAMData, ZZ, alpha);
  gg_p += gg3(1) * (uu3 * mB_p);

//  if( alpha==0 ) 
//   gg_p *= 22.0;
//   std::cout << "\nGP " << alpha << ":" << gg_p;
  return gg_p;
}

void   KNLpAutRotTestFunctional::funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  col.interpolate(vv3Data, vv3);
  col.jotf_x_x(A_x, par, vv3Data, ZZ);
  func = !A_x * uu3;

  temp = !mB * uu3;
  col.jotf_trivialKernelOnMesh_x<true>(DxPhi, par, solMSHData, temp);
  func += gg3(0) * DxPhi;
  col.jotf_mB_x(mB_x, par, phiData, ZZ);
  temp = !mB_x * uu3;
  func += gg3(0) * temp;

  col.jotf_mB_x(mB_x, par, LAMData, ZZ);
  temp = !mB_x * uu3;
  func += gg3(1) * temp;

  rotbord<true>(DxLAM, col, uu3, Re, Im);
  func += gg3(1) * DxLAM;

  col.jotf_trivialKernelOnMesh_x<true>(DxPhi, par, solMSHData, vv3);
  col.star(temp, DxPhi);
  func += hh3(0) * temp;

  rotbord<true>(DxLAM, col, vv3, Re, Im);
  col.star(temp, DxLAM);
  func += hh3(1) * temp;
//   func *= 2;
}

void   KNLpAutRotTestFunctional::kernel(KNVector& phi)
{
  phi = vv3;
}

/// -----------------------------------------------------------------------
/// test function for FOLD BIFURCATIONS in autonomous systems with SYMMETRY
/// EXTENDED
/// -----------------------------------------------------------------------

KNLpAutRotTestFunctional2::KNLpAutRotTestFunctional2(KNDdeBvpCollocation& col, KNArray1D<size_t> CRe, KNArray1D<size_t> CIm, double Z) :
    ZZ(Z),
    Re(CRe), Im(CIm),
    AHAT(NDIM*(NDEG*NINT + 1), 0, 3, NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    A_p(NDIM*(NDEG*NINT + 1)),
    A_x('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    mB('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    mB_p(NDIM*(NDEG*NINT + 1)),
    mB_x('R', NDIM*(NDEG*NINT + 1), NDIM*(NDEG*NINT + 1)*NTAU*NDIM*(NDEG + 1)),
    phi(NDIM*(NDEG*NINT + 1)),
    DpPhi(NDIM*(NDEG*NINT + 1)),
    DxPhi(NDIM*(NDEG*NINT + 1)),
    LAM(NDIM*(NDEG*NINT + 1)),
    DxLAM(NDIM*(NDEG*NINT + 1)),
    uu3(NDIM*(NDEG*NINT + 1)),
    vv3(NDIM*(NDEG*NINT + 1)),
    gg3((size_t)3),
    hh3((size_t)3),
    one3((size_t)3),
    uu1(NDIM*(NDEG*NINT + 1)),
    vv1(NDIM*(NDEG*NINT + 1)),
    gg1((size_t)2),
    hh1((size_t)2),
    one1((size_t)2),
    uu2(NDIM*(NDEG*NINT + 1)),
    vv2(NDIM*(NDEG*NINT + 1)),
    gg2((size_t)2),
    hh2((size_t)2),
    one2((size_t)2),
    rhs(NDIM*(NDEG*NINT + 1)),
    temp(NDIM*(NDEG*NINT + 1)),
    vv1Data(NDIM, 2*NTAU + 1, NDEG*NINT),
    vv2Data(NDIM, 2*NTAU + 1, NDEG*NINT),
    vv3Data(NDIM, 2*NTAU + 1, NDEG*NINT),
    solMSHData(NDIM, 2*NTAU + 1, NDEG*NINT + 1)
{
  first = true;
}

KNLpAutRotTestFunctional2::~KNLpAutRotTestFunctional2()
{}

void KNLpAutRotTestFunctional2::init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  // creating the matrix
  col.jotf_x(AHAT.getA11(), par, ZZ);
  col.jotf_mB(mB, par, ZZ);

  vv1.random();
  uu1.random();
  vv2.random();
  uu2.random();

  vv1 /= sqrt(vv1 * vv1);
  uu1 /= sqrt(uu1 * uu1);
  vv2 -= (vv1 * vv2) * vv1;
  uu2 -= (uu1 * uu2) * uu1;
  vv2 /= sqrt(vv2 * vv2);
  uu2 /= sqrt(uu2 * uu2);

  // looking for the singular vectors (geometric multiplicity)
  AHAT.getA13(0) = uu1;
  AHAT.getA13(1) = uu2;
  AHAT.getA31(0) = vv1;
  AHAT.getA31(1) = vv2;

  AHAT.getA33().clear();
  rhs.clear();
  one1(0) = 1.0;
  one1(1) = 0.0;
  one2(0) = 0.0;
  one2(1) = 1.0;
  for (size_t i = 0; i < kernIter; i++)
  {
    AHAT.solve(2, vv1, gg1, rhs, one1);
    AHAT.solve(2, vv2, gg2, rhs, one2);
    AHAT.solveTr(2, uu1, hh1, rhs, one1);
    AHAT.solveTr(2, uu2, hh2, rhs, one2);
    const double nrm_v1 = (1.0 / sqrt(vv1 * vv1));
    const double nrm_v2 = (1.0 / sqrt(vv2 * vv2));
    const double nrm_u1 = (1.0 / sqrt(uu1 * uu1));
    const double nrm_u2 = (1.0 / sqrt(uu2 * uu2));
    AHAT.getA13(0) = nrm_u1 * uu1;
    AHAT.getA13(1) = nrm_u2 * uu2;
    AHAT.getA31(0) = nrm_v1 * vv1;
    AHAT.getA31(1) = nrm_v2 * vv2;
  }

  AHAT.getA13(0) = mB * vv1;
  AHAT.getA13(1) = mB * vv2;

  vv3.random();
  uu3.random();
  vv3 /= sqrt(vv1 * vv1);
  uu3 /= sqrt(uu1 * uu1);
  AHAT.getA13(2) = uu3;
  AHAT.getA31(2) = vv3;
  // generating the non-trivial kernel
  one3(0) = 0.0;
  one3(1) = 0.0;
  one3(2) = 1.0;
  for (size_t i = 0; i < kernIter; i++)
  {
    AHAT.solve(3, vv3, gg3, rhs, one3);
    AHAT.solveTr(3, uu3, hh3, rhs, one3);
    const double nrm_v = (1.0 / sqrt(vv3 * vv3 + gg3(0) * gg3(0) + gg3(1) * gg3(1)));
    const double nrm_u = (1.0 / sqrt(uu3 * uu3 + hh3(0) * hh3(0) + hh3(1) * hh3(1)));
    AHAT.getA31(2)   = nrm_v * vv3;
    AHAT.getA33(2, 0) = nrm_v * gg3(0);
    AHAT.getA33(2, 1) = nrm_v * gg3(1);
    AHAT.getA13(2)   = nrm_u * uu3;
    AHAT.getA33(0, 2) = nrm_u * hh3(0);
    AHAT.getA33(1, 2) = nrm_u * hh3(1);
  }
  std::cout << "TF: " << gg3(2) << ", " << hh3(2) << "\n";
}

double KNLpAutRotTestFunctional2::funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol)
{
  if (first)
  {
    init(col, par, sol);
    first = false;
  }

  col.jotf_x(AHAT.getA11(), par, ZZ);
  col.jotf_mB(mB, par, ZZ);

  AHAT.getA13(0) = uu1;
  AHAT.getA13(1) = uu2;

  AHAT.solve(2, vv1, gg1, rhs, one1);
  AHAT.solve(2, vv2, gg2, rhs, one2);
  AHAT.solveTr(2, uu1, hh1, rhs, one1);
  AHAT.solveTr(2, uu2, hh2, rhs, one2);
  const double nrm_v1 = (1.0 / sqrt(vv1 * vv1));
  const double nrm_v2 = (1.0 / sqrt(vv2 * vv2));
  const double nrm_u1 = (1.0 / sqrt(uu1 * uu1));
  const double nrm_u2 = (1.0 / sqrt(uu2 * uu2));
  AHAT.getA13(0) = nrm_u1 * uu1;
  AHAT.getA13(1) = nrm_u2 * uu2;
  AHAT.getA31(0) = nrm_v1 * vv1;
  AHAT.getA31(1) = nrm_v2 * vv2;

  AHAT.getA13(0) = mB * vv1;
  AHAT.getA13(1) = mB * vv2;

//  temp = AHAT.getA11() * LAM;
//  std::cout<<"temp: "<<temp*temp<<" LAM: "<<LAM*LAM<<"\n";

  AHAT.solve(3, vv3, gg3, rhs, one3);
  AHAT.solveTr(3, uu3, hh3, rhs, one3);
  const double nrm_v = (1.0 / sqrt(vv3 * vv3 + gg3(0) * gg3(0) + gg3(1) * gg3(1)));
  const double nrm_u = (1.0 / sqrt(uu3 * uu3 + hh3(0) * hh3(0) + hh3(1) * hh3(1)));
  AHAT.getA31(2)   = nrm_v * vv3;
  AHAT.getA33(2, 0) = nrm_v * gg3(0);
  AHAT.getA33(2, 1) = nrm_v * gg3(1);
  AHAT.getA13(2)   = nrm_u * uu3;
  AHAT.getA33(0, 2) = nrm_u * hh3(0);
  AHAT.getA33(1, 2) = nrm_u * hh3(1);

  // for subsequent use
  col.interpolate(vv1Data, vv1);
  col.interpolate(vv2Data, vv2);
  col.interpolate(vv3Data, vv3);

  if (gg3(2) > 0.0) std::cout << "\t+++\n";
  else               std::cout << "\t---\n";
  return gg3(2);
}

double KNLpAutRotTestFunctional2::funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/, size_t alpha)
{
  col.jotf_x_p(A_p, par, vv3Data, ZZ, alpha);
  col.jotf_mB_p(mB_p, par, vv1Data, ZZ, alpha);
  temp = gg3(0) * mB_p;
  col.jotf_mB_p(mB_p, par, vv2Data, ZZ, alpha);
  temp += gg3(1) * mB_p;

  double gg_p = (uu3 * A_p) + (uu3 * temp);
  std::cout << "\nGP " << alpha << ":" << gg_p;
  return gg_p;
}

void   KNLpAutRotTestFunctional2::funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& /*sol*/)
{
  col.jotf_x_x(A_x, par, vv3Data, ZZ);
  func = !A_x * uu3;
  col.jotf_mB_x(mB_x, par, vv1Data, ZZ);
  temp = !mB_x * uu3;
  func += gg3(0) * temp;
  col.jotf_mB_x(mB_x, par, vv2Data, ZZ);
  temp = !mB_x * uu3;
  func += gg3(1) * temp;
}

void   KNLpAutRotTestFunctional2::kernel(KNVector& phi)
{
  phi = vv3;
}
