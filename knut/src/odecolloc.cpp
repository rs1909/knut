// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "knerror.h"
#include "odecolloc.h"
#include "exprsystem.h"
#include "matrix.h"
#include "polynomial.h"

#include <cmath>

#ifdef DEBUG
#include <iostream>
#include <fstream>
#endif /*DEBUG*/

#define NDIM ndim
#define NPAR npar
#define NINT nint
#define NDEG ndeg

KNOdeBvpCollocation::KNOdeBvpCollocation(KNExprSystem& _sys, const size_t _nint, const size_t _ndeg)
    : KNAbstractBvpCollocation(_sys, _nint, _ndeg),
    tt(2, ndeg + 1, ndeg*nint),
    ttMSH(2, ndeg + 1, ndeg*nint + 1),
    solData(ndim, 2 + 1, ndeg*nint),
    p_fx(ndim, nint*ndeg), p_dfx(ndim,ndim, nint*ndeg), p_dfp(ndim,1, nint*ndeg),
    p_fxRe(p_fx), p_fxIm(ndim, nint*ndeg),
    p_dfxRe(p_dfx), p_dfxIm(ndim,ndim, nint*ndeg),
    p_fxMSH(ndim, nint*ndeg+1), p_dfxMSH(ndim,ndim, nint*ndeg+1), p_dfpMSH(ndim,1, nint*ndeg+1),
    p_dummy(0,0,nint*ndeg+1)
{ }

/// Determines the structure of the sparse matrices
/// Input:
///    mesh    : the mesh [0,1]
///    meshINT : internal mesh of the collocation subinterval (repr. points) [0,1]
///    col     : the collocation points in [0,1]
/// Output:
///    time(i)  : (double) all the collocation point
///    ttMSH(j,l,i) : j == 0, interp the solution
///                   j == 1, interp derivative
///    tt(j,l,i) : j == 0, interp the solution
///                j == 1, interp derivative
///
void KNOdeBvpCollocation::init(const KNVector& sol, const KNVector& /*par*/)
{
  for (size_t i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (size_t j = 0; j < NDEG; j++)
    {
      const size_t idx = j + i * NDEG;
      const double h0 = mesh(i + 1) - mesh(i);
      time(idx) = mesh(i) + h0 * col(j);
      timeMSH(idx) = mesh(i) + h0 * meshINT(j);
    }
  }
  timeMSH(NDEG*NINT) = 1.0;

  for (size_t i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (size_t j = 0; j < NDEG; j++)
    {
      const size_t idx = j + i * NDEG;
      const double h0 = mesh(i + 1) - mesh(i);

      poly_lgr(meshINT, out, col(j));
      for (size_t l = 0; l < NDEG + 1; l++)
      {
        tt(0, l, idx) = out(l);
      }
      poly_dlg(meshINT, out, col(j));
      for (size_t k = 0; k < NDEG + 1; k++)
      {
        tt(1, k, idx) = out(k) / h0;
      }
      poly_lgr(meshINT, out, meshINT(j));
      for (size_t l = 0; l < NDEG + 1; l++)
      {
        ttMSH(0, l, idx) = out(l);
      }
      poly_dlg(meshINT, out, meshINT(j));
      for (size_t l = 0; l < NDEG + 1; l++)
      {
        ttMSH(1, l, idx) = out(l) / h0;
      }
    }
  }

  // making periodic
  for (size_t l = 0; l < NDEG + 1; l++)
  {
    ttMSH(0, l, NDEG*NINT) = ttMSH(0, l, 0);
    ttMSH(1, l, NDEG*NINT) = ttMSH(1, l, 0);
  }
  // evaluate the solution and its derivatives at the collocation points
  interpolate(solData, sol);
}

void KNOdeBvpCollocation::interpolate(KNArray3D<double>& solData, const KNVector& sol)
{
  for (size_t i = 0; i < NINT; ++i)
  {
    for (size_t j = 0; j < NDEG; ++j)
    {
      const size_t idx = j + i * NDEG;
      for (size_t k = 0; k < NDIM; k++)
      {
        solData(k, 0, idx) = 0.0;
        solData(k, 1, idx) = 0.0;
        for (size_t l = 0; l < NDEG + 1; l++)
        {
          solData(k, 0, idx) += sol(k + NDIM * (l + i * NDEG)) * tt(0, l, idx);
          solData(k, 1, idx) += sol(k + NDIM * (l + i * NDEG)) * tt(1, l, idx);
        }
      }
    }
  }
}

void KNOdeBvpCollocation::interpolateOnMesh(KNArray3D<double>& solData, const KNVector& sol)
{
  for (size_t i = 0; i < NINT; ++i)
  {
    for (size_t j = 0; j < NDEG; ++j)
    {
      const size_t idx = j + i * NDEG;
      for (size_t k = 0; k < NDIM; k++)
      {
        // std::cout<<"InterR: "<<idx<<", "<<out(idx).row()<<", "<<out(idx).col()<<", "<<k<<", "<<p<<"\n";
        solData(k, 0, idx) = 0.0;
        solData(k, 1, idx) = 0.0;
        for (size_t l = 0; l < NDEG + 1; l++)
        {
          solData(k, 0, idx) += sol(k + NDIM * (l + i * NDEG)) * ttMSH(0, l, idx);
          solData(k, 1, idx) += sol(k + NDIM * (l + i * NDEG)) * ttMSH(1, l, idx);
        }
      }
    }
  }
}

// in complex form on 2*i th places are the reals and on 2*i+1 th places the imaginary parts

void KNOdeBvpCollocation::interpolateComplex(KNArray3D<double>& solDataRe, KNArray3D<double>& solDataIm, const KNVector& sol)
{
  for (size_t i = 0; i < NINT; ++i)
  {
    for (size_t j = 0; j < NDEG; ++j)
    {
      const size_t idx = j + i * NDEG;
      for (size_t k = 0; k < NDIM; k++)
      {
        solDataRe(k, 0, idx) = 0.0;
        solDataIm(k, 0, idx) = 0.0;
        solDataRe(k, 1, idx) = 0.0;
        solDataIm(k, 1, idx) = 0.0;
        for (size_t l = 0; l < NDEG + 1; l++)
        {
          solDataRe(k, 0, idx) += sol(2 * (k + NDIM * (l + i * NDEG))) * tt(0, l, idx);
          solDataIm(k, 0, idx) += sol(2 * (k + NDIM * (l + i * NDEG)) + 1) * tt(0, l, idx);
          solDataRe(k, 1, idx) += sol(2 * (k + NDIM * (l + i * NDEG))) * tt(1, l, idx);
          solDataIm(k, 1, idx) += sol(2 * (k + NDIM * (l + i * NDEG)) + 1) * tt(1, l, idx);
        }
      }
    }
  }
}

void KNOdeBvpCollocation::rightHandSide(KNVector& rhs, const KNVector& par, const KNVector& sol)
{
  // boundary conditions
  for (size_t r = 0; r < NDIM; r++)
  {
    rhs(r) = sol(r + NDIM * NDEG * NINT) - sol(r);
  }

  sys->p_rhs(p_fx, time, solData, par, 0);
  for (size_t k = 0; k < NDIM; k++)
  {
    for (size_t idx = 0; idx < NDEG*NINT; ++idx)
    {
      rhs(NDIM + k + NDIM*idx) = par(0) * p_fx(k, idx) - solData(k, 1, idx);
    }
  }
}

void KNOdeBvpCollocation::rightHandSide_p(KNVector& rhs, const KNVector& par, const KNVector& /*sol*/, size_t alpha)
{
  // boundary conditions
  for (size_t r = 0; r < NDIM; r++)
  {
    rhs(r) = 0.0;
  }

  if (alpha == 0)
  {
    sys->p_rhs(p_fx, time, solData, par, 0);

    for (size_t p = 0; p < NDIM; p++)
    {
      for (size_t idx = 0; idx < NDEG*NINT; ++idx)
      {
        rhs(NDIM + p + NDIM*idx) = -p_fx(p, idx);
      }
    }
  }
  else
  {
    size_t nx=0, vx=0, np=1, vp=alpha;
    sys->p_deri(p_dfp, time, solData, par, 0, nx, &vx, np, &vp, p_dummy);

    for (size_t k = 0; k < NDIM; k++)
    {
      for (size_t idx = 0; idx < NDEG*NINT; ++idx)
      {
        rhs(NDIM + k + NDIM*idx) = -par(0) * p_dfp(k, 0, idx);
      }
    }
  }
}

void KNOdeBvpCollocation::rightHandSide_x(KNSparseMatrix& A, const KNVector& par, const KNVector& /*sol*/)
{
  RHS_jacobian<true>(A, par);
}

void KNOdeBvpCollocation::jacobianOfStability(KNSparseMatrix& A, const KNVector& par)
{
  RHS_jacobian<false>(A, par);
}

template <bool periodic>
void KNOdeBvpCollocation::RHS_jacobian(KNSparseMatrix& A, const KNVector& par )
{
  A.clear('R');

  // boundary conditions
  for (size_t r = 0; r < NDIM; r++)
  {
    if (periodic) A.newLine(2);
    else          A.newLine(1);
    A.writeIndex(r, 0) = static_cast<int>( r );
    A.writeData(r, 0) = 1.0;
    if (periodic)
    {
      A.writeIndex(r, 1) = static_cast<int>(r + NDIM * NDEG * NINT);
      A.writeData(r, 1) = -1.0;
    }
  }

  // MUST PRESERVE THE ORDER. DO NOT CHANGE THE LOOPS
  for (size_t idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (size_t r = 0; r < NDIM; r++)
    {
      A.newLine(NDIM*(NDEG + 1));
    }
  }

  size_t nx = 1, vx = 0, np = 0, vp = 0;
  sys->p_deri(p_dfx, time, solData, par, 0, nx, &vx, np, &vp, p_dummy);
  for (size_t i = 0; i < NINT; ++i)
  {
    for (size_t j = 0; j < NDEG; ++j)
    {
      const size_t idx = j + i * NDEG;
      for (size_t l = 0; l < NDEG + 1; l++) // degree
      {
        for (size_t p = 0; p < NDIM; p++)   // row
        {
          for (size_t q = 0; q < NDIM; q++)   //column
          {
            WRIDX(A, idx, p, l, q) = static_cast<int>(q + NDIM * (l + NDEG * i));
            if (p == q) 
              WRDAT(A, idx, p, l, q) += tt(1, l, idx) - par(0) * p_dfx(p, q, idx) * tt(0, l, idx);
            else 
              WRDAT(A, idx, p, l, q) += 0.0 - par(0) * p_dfx(p, q, idx) * tt(0, l, idx);
          }
        }
      }
    }
  }
}
