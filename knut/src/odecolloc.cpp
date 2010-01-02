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
#include "system.h"
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

ODEColloc::ODEColloc(System& _sys, const int _nint, const int _ndeg)
    : PerSolColloc(_sys, _nint, _ndeg),
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
void ODEColloc::init(const Vector& sol, const Vector& par)
{
  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;
      const double h0 = mesh(i + 1) - mesh(i);
      time(idx) = mesh(i) + h0 * col(j);
      timeMSH(idx) = mesh(i) + h0 * meshINT(j);
    }
  }
  timeMSH(NDEG*NINT) = 1.0;

  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;
      const double h0 = mesh(i + 1) - mesh(i);

      poly_lgr(meshINT, out, col(j));
      for (int l = 0; l < NDEG + 1; l++)
      {
        tt(0, l, idx) = out(l);
      }
      poly_dlg(meshINT, out, col(j));
      for (int k = 0; k < NDEG + 1; k++)
      {
        tt(1, k, idx) = out(k) / h0;
      }
      poly_lgr(meshINT, out, meshINT(j));
      for (int l = 0; l < NDEG + 1; l++)
      {
        ttMSH(0, l, idx) = out(l);
      }
      poly_dlg(meshINT, out, meshINT(j));
      for (int l = 0; l < NDEG + 1; l++)
      {
        ttMSH(1, l, idx) = out(l) / h0;
      }
    }
  }

  // making periodic
  for (int l = 0; l < NDEG + 1; l++)
  {
    ttMSH(0, l, NDEG*NINT) = ttMSH(0, l, 0);
    ttMSH(1, l, NDEG*NINT) = ttMSH(1, l, 0);
  }
  // evaluate the solution and its derivatives at the collocation points
  interpolate(solData, sol);
}

void ODEColloc::interpolate(Array3D<double>& solData, const Vector& sol)
{
  for (int i = 0; i < NINT; ++i)
  {
    for (int j = 0; j < NDEG; ++j)
    {
      const int idx = j + i * NDEG;
      for (int k = 0; k < NDIM; k++)
      {
        solData(k, 0, idx) = 0.0;
        solData(k, 1, idx) = 0.0;
        for (int l = 0; l < NDEG + 1; l++)
        {
          solData(k, 0, idx) += sol(k + NDIM * (l + i * NDEG)) * tt(0, l, idx);
          solData(k, 1, idx) += sol(k + NDIM * (l + i * NDEG)) * tt(1, l, idx);
        }
      }
    }
  }
}

void ODEColloc::interpolateOnMesh(Array3D<double>& solData, const Vector& sol)
{
  for (int i = 0; i < NINT; ++i)
  {
    for (int j = 0; j < NDEG; ++j)
    {
      const int idx = j + i * NDEG;
      for (int k = 0; k < NDIM; k++)
      {
        // std::cout<<"InterR: "<<idx<<", "<<out(idx).row()<<", "<<out(idx).col()<<", "<<k<<", "<<p<<"\n";
        solData(k, 0, idx) = 0.0;
        solData(k, 1, idx) = 0.0;
        for (int l = 0; l < NDEG + 1; l++)
        {
          solData(k, 0, idx) += sol(k + NDIM * (l + i * NDEG)) * ttMSH(0, l, idx);
          solData(k, 1, idx) += sol(k + NDIM * (l + i * NDEG)) * ttMSH(1, l, idx);
        }
      }
    }
  }
}

// in complex form on 2*i th places are the reals and on 2*i+1 th places the imaginary parts

void ODEColloc::interpolateComplex(Array3D<double>& solDataRe, Array3D<double>& solDataIm, const Vector& sol)
{
  for (int i = 0; i < NINT; ++i)
  {
    for (int j = 0; j < NDEG; ++j)
    {
      const int idx = j + i * NDEG;
      for (int k = 0; k < NDIM; k++)
      {
        solDataRe(k, 0, idx) = 0.0;
        solDataIm(k, 0, idx) = 0.0;
        solDataRe(k, 1, idx) = 0.0;
        solDataIm(k, 1, idx) = 0.0;
        for (int l = 0; l < NDEG + 1; l++)
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

void ODEColloc::RHS(Vector& rhs, const Vector& par, const Vector& sol)
{
  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    rhs(r) = sol(r + NDIM * NDEG * NINT) - sol(r);
  }

  sys->p_rhs(p_fx, time, solData, par, 0);
  for (int k = 0; k < NDIM; k++)
  {
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      rhs(NDIM + k + NDIM*idx) = par(0) * p_fx(k, idx) - solData(k, 1, idx);
    }
  }
}

void ODEColloc::RHS_p(Vector& rhs, const Vector& par, const Vector& /*sol*/, int alpha)
{
  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    rhs(r) = 0.0;
  }

  if (alpha == 0)
  {
    sys->p_rhs(p_fx, time, solData, par, 0);

    for (int p = 0; p < NDIM; p++)
    {
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        rhs(NDIM + p + NDIM*idx) = -p_fx(p, idx);
      }
    }
  }
  else
  {
    int nx=0, vx=0, np=1, vp=alpha;
    sys->p_deri(p_dfp, time, solData, par, 0, nx, &vx, np, &vp, p_dummy);

    for (int k = 0; k < NDIM; k++)
    {
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        rhs(NDIM + k + NDIM*idx) = -par(0) * p_dfp(k, 0, idx);
      }
    }
  }
}

void ODEColloc::RHS_x(SpMatrix& A, const Vector& par, const Vector& /*sol*/)
{
  RHS_jacobian<true>(A, par);
}

void ODEColloc::StabJac(SpMatrix& A, const Vector& par)
{
  RHS_jacobian<false>(A, par);
}

template <bool periodic>
void ODEColloc::RHS_jacobian(SpMatrix& A, const Vector& par )
{
  A.clear('R');

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    if (periodic) A.NewL(2);
    else          A.NewL(1);
    A.WrLi(r, 0) = r;
    A.WrLx(r, 0) = 1.0;
    if (periodic)
    {
      A.WrLi(r, 1) = r + NDIM * NDEG * NINT;
      A.WrLx(r, 1) = -1.0;
    }
  }

  // MUST PRESERVE THE ORDER. DO NOT CHANGE THE LOOPS
  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      A.NewL(NDIM*(NDEG + 1));
    }
  }

  int nx = 1, vx = 0, np = 0, vp;
  sys->p_deri(p_dfx, time, solData, par, 0, nx, &vx, np, &vp, p_dummy);
  for (int i = 0; i < NINT; ++i)
  {
    for (int j = 0; j < NDEG; ++j)
    {
      const int idx = j + i * NDEG;
      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
          {
            WRIDX(A, idx, p, l, q) = q + NDIM * (l + NDEG * i);
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
