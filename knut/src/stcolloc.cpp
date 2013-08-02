// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2013 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "stcolloc.h"
#include "exprsystem.h"
#include "spmatrix.h"

#define NDIM ndim
#define NPAR npar
#define NTAU ntau

KNSteadyStateJacobian::KNSteadyStateJacobian(KNSystem& _sys) : sys(&_sys),
    ndim(_sys.ndim()), npar(_sys.npar()), ntau(_sys.ntau()),
    time((size_t)1), solData(ndim, 2*ntau + 1, 1),
    p_fx(ndim, 1), p_dfx(ndim,ndim, 1), p_dfp(ndim,1, 1),
    p_fxRe(p_fx), p_fxIm(ndim, 1),
    p_dfxRe(p_dfx), p_dfxIm(ndim,ndim, 1),
    p_dummy(0,0,1) 
{ }

void KNSteadyStateJacobian::init(const KNVector& sol, const KNVector& /*par*/)
{
  for (size_t k=0; k < ntau; k++)
  {
    for (size_t p=0; p < ndim; p++)
    {
      solData(p, k, 0) = sol( p );
    }
  }
  time(0) = 0.0;
}

void KNSteadyStateJacobian::rightHandSide(KNVector& rhs, const KNVector& par, const KNVector& /*sol*/)
{
  sys->p_rhs(p_fx, time, solData, par, 0);
  for (size_t q=0; q < ndim; q++) rhs(q) = p_fx(q,0);
}

void KNSteadyStateJacobian::rightHandSide_p(KNVector& rhs, const KNVector& par, const KNVector& /*sol*/, size_t alpha)
{
  if (alpha == 0)
  {
//    sys->p_rhs(p_fx, time, solData, par, 0);
//    for (size_t q=0; q < ndim; q++) rhs(q) = p_fx(q,0);
  } else
  {
    size_t nx=0, vx, np=1, vp = alpha;
    sys->p_deri(p_dfp, time, solData, par, 0, nx, &vx, np, &vp, p_dummy);
    for (size_t q=0; q < ndim; q++) rhs(q) = -p_dfp(q,0,0);
  }
}

void KNSteadyStateJacobian::rightHandSide_x(KNSparseMatrix& A, const KNVector& par, const KNVector& sol)
{
  A.clear('R');
  for (size_t r = 0; r < ndim; r++) A.newLine(NDIM);
  for (size_t k = 0; k < NTAU; k++)
  {
    size_t nx = 1, vx = k, np = 0, vp;
    sys->p_deri(p_dfx, time, solData, par, 0, nx, &vx, np, &vp, p_dummy);
    for (size_t p=0; p < ndim; p++)
    {
      for (size_t q=0; q < ndim; q++)
      {
        A.writeIndex(p, q) = q;
        A.writeData(p, q) -= p_dfx(p,q,0);
//        std::cerr << "[" << p << "," << q << "]\n";
      }
    }
  }
//  par.print();
//  sol.print();
//  A.print();   
}
