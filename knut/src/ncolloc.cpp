// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "knerror.h"
#include "ncolloc.h"
#include "system.h"
#include "matrix.h"
#include "polynomial.h"

#include <cmath>

#ifdef DEBUG
#include <iostream>
#include <fstream>
#endif /*DEBUG*/

// #define MADD

// support routines


//
// the NColloc class
//

#define NDIM ndim
#define NPAR npar
#define NTAU ntau
#define NINT nint
#define NDEG ndeg
#define NMAT nmat

NColloc::NColloc(System& _sys, const int _nint, const int _ndeg, int _nmat) :
    PerSolColloc(_sys, _nint, _ndeg),
    ntau(_sys.ntau()), nmat(_nmat),

    kk(ntau + 1, nint*ndeg), ee(ntau + 1, nint*ndeg),
    dd(ntau + 1, nint*ndeg), rr(ntau + 1, nint*ndeg),

    kkS(ntau + 1, nint*ndeg), eeS(ntau + 1, nint*ndeg),
    rrS(ntau + 1, nint*ndeg), ddS(ntau + 1, nint*ndeg),
    mmS(ntau + 1, nint*ndeg), szS(2, nint*ndeg),

    kkI(ntau + 1, nint*ndeg), eeI(ntau + 1, nint*ndeg),
    rrI(ntau + 1, nint*ndeg), ddI(ntau + 1, nint*ndeg),
    mmI(ntau + 1, nint*ndeg), szI(nmat + 1, nint*ndeg),

    kkMSH(ntau, nint*ndeg + 1),tt(2*ntau + 1, ndeg + 1, ndeg*nint),
    ttMSH(2*ntau, ndeg + 1, ndeg*nint + 1),
    solData(ndim, 2*ntau + 1, ndeg*nint),
    p_tau(ntau, nint*ndeg), p_dtau(ntau, nint*ndeg),
    p_fx(ndim, nint*ndeg), p_dfx(ndim,ndim, nint*ndeg), p_dfp(ndim,1, nint*ndeg),
    p_fxRe(p_fx), p_fxIm(ndim, nint*ndeg),
    p_dfxRe(p_dfx), p_dfxIm(ndim,ndim, nint*ndeg),
    p_tauMSH(ntau, nint*ndeg+1), p_dtauMSH(ntau, nint*ndeg+1),
    p_fxMSH(ndim, nint*ndeg+1), p_dfxMSH(ndim,ndim, nint*ndeg+1), p_dfpMSH(ndim,1, nint*ndeg+1),
    p_dummy(0,0,nint*ndeg+1)
{
#ifdef DEBUG
  count_reset();
#endif
}


/// Determines the structure of the sparse matrices
/// Input:
///    mesh    : the mesh [0,1]
///    meshINT : internal mesh of the collocation subinterval (repr. points) [0,1]
///    col     : the collocation points in [0,1]
/// Output:
///    time(i)  : (double) all the collocation point
///    kk(j,i)  : (int) the collocation interval into which the delayed time falls
///               this must be in [0,NDEG*NINT)
///               j == 0   for the derivative d/dt x(t)
///               j == k+1 for x(t - \tau_k) (most of the time \tau_0 = 0)
///    ee(j,i)  : kk(ee(j,i),i) is in increasing order
///    dd(j,i)  : if there are neighbouring subintervals, then this increases by one
///    rr(j,i)  : counts the real number of intervals. If i & i+1 overlap they have the same rr number.
///    kkS(j,i) : the kk for two separate matrices, e.g., when calculating the; NOT used
///    eeS(j,i) : see ee
///    ddS(j,i) : see dd
///    rrS(j,i) : see rr
///    mmS(j,i) : if == 0 then matrix A, if == 1 the \Sum_k B_k; NOT used
///    szS(0,i) : size of A
///    szS(1,i) : size of \Sum B_k ; used in CharJac_mB
///    kkI(j,i) : kk without the mod operation
///    eeI(j,i) : sorting kkI
///    ddI(j,i) : see dd
///    rrI(j,i) : see rr
///    mmI(j,i) : the exponent of \mu in (A - \Sum_k \mu^k B_k)
///    szI(j,i) : in Stability
///    kkMSH(j,i) : j = k, the delay. Here, we don't shift as in others
///                 this is same as kk, but at the representation points
///    tt(j,l,i) : j == 0, interp derivative
///                j == k+1, interp at the k-th delay;
///                j == k+NDEG+1; interp of derivative at the k-th delay
///    ttMSH(j,l,i) : j = k, the delay. Here, we don't shift as in others
///                   it contains only the middle of tt j=1...NDEG at indices j=0...NDEG-1
///
void NColloc::Init(const Vector& sol, const Vector& par)
{
  Array1D<double> t(NTAU+1);
  Array1D<double> tMSH(NTAU);

  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;
      const double h0 = mesh(i + 1) - mesh(i);
      time(idx) = mesh(i) + h0 * col(j);
//      if (!std::isfinite(mesh(i)) || !std::isfinite(col(j))) std::cout << mesh(i) << ", " << col(j) << "\n";
      timeMSH(idx) = mesh(i) + h0 * meshINT(j);
    }
  }
  timeMSH(NDEG*NINT) = 1.0;

  sys->p_tau(p_tau, time, par); // this will be rescaled, so have to call it again

  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;
      const double h0 = mesh(i + 1) - mesh(i);

      /* first, the derivative: x'(t)*/
      t(0) = time(idx);
      kk(0, idx) = i;
      kkS(0, idx) = i;
      kkI(0, idx) = i;
      poly_dlg(meshINT, out, col(j));
      // out.Print(); std::cout<<"SS\n";
      for (int k = 0; k < NDEG + 1; k++)
      {
        tt(0, k, idx) = out(k) / h0;
      }

      /* second, the right-hand side */
      for (int k = 0; k < NTAU; k++)
      {
        p_tau(k,idx) /= par(0);
        P_ERROR_X3(p_tau(k,idx) <= NMAT, "The scaled delay became greater then NMAT times the period. k=", k, ".");
        P_ERROR_X3(p_tau(k,idx) >= 0.0, "Either the delay or the period became negative. k=", k, ".");
        t(1+k) = (t(0) - p_tau(k,idx)) - floor(t(0) - p_tau(k,idx));  // nem szetvalasztott

        // binary search for in which interval is t-p_tau(k,idx)
        const int low = meshlookup(mesh, t(1+k));
        kk(k + 1, idx) = low;

        if (t(0) - p_tau(k,idx) >= 0) kkS(k + 1, idx) = low;
        else kkS(k + 1, idx) = low - NINT;

        kkI(k + 1, idx) = low + NINT * static_cast<int>(floor(t(0) - p_tau(k,idx)));

//        std::cout << "low " << low << " t(0) " << t(0) << " p_tau(k,idx) " << p_tau(k,idx);
        const double hk = mesh(low + 1) - mesh(low);
        // x(t-\tau_i)
        poly_lgr(meshINT, out, (t(1+k) - mesh(low)) / hk);
        for (int l = 0; l < NDEG + 1; l++)
        {
          tt(1 + k, l, idx) = out(l);
        }
        // x'(t-\tau_i)
        poly_dlg(meshINT, out, (t(1+k) - mesh(low)) / hk);
        for (int l = 0; l < NDEG + 1; l++)
        {
          tt(NTAU + 1 + k, l, idx) = out(l) / hk;
        }
        // creating interpolation at the representation points
        tMSH(k) = mesh(i) + h0 * meshINT(j) - p_tau(k,idx) - floor(mesh(i) + h0 * meshINT(j) - p_tau(k,idx));
        const int lowMSH = meshlookup(mesh, tMSH(k));
        kkMSH(k, idx) = lowMSH;
        const double hkMSH = mesh(lowMSH + 1) - mesh(lowMSH);
        poly_lgr(meshINT, out, (tMSH(k) - mesh(lowMSH)) / hkMSH);
        for (int l = 0; l < NDEG + 1; l++)
        {
          ttMSH(k, l, idx) = out(l);
        }
        // x'(t-\tau_i)
        poly_dlg(meshINT, out, (tMSH(k) - mesh(lowMSH)) / hkMSH);
        for (int l = 0; l < NDEG + 1; l++)
        {
          ttMSH(NTAU + k, l, idx) = out(l) / hkMSH;
        }
        // REPR ENDS
      }

      // sorting kk

      for (int r = 0; r < NTAU + 1; r++)
      {
        ee(r, idx) = r;
        eeS(r, idx) = r;
        eeI(r, idx) = r;
      }

      int tmp;
      // probably there is a better algorithm. It is n^2, but NTAU is generally small
      for (int r = 1; r < NTAU + 1; r++)
      {
        for (int s = r; s < NTAU + 1; s++)
        {
          if (kk(ee(r - 1, idx), idx) > kk(ee(s, idx), idx))
          {
            tmp = ee(r - 1, idx);
            ee(r - 1, idx) = ee(s, idx);
            ee(s, idx) = tmp;
          }
          if (kkS(eeS(r - 1, idx), idx) > kkS(eeS(s, idx), idx))
          {
            tmp = eeS(r - 1, idx);
            eeS(r - 1, idx) = eeS(s, idx);
            eeS(s, idx) = tmp;
          }
          if (kkI(eeI(r - 1, idx), idx) > kkI(eeI(s, idx), idx))
          {
            tmp = eeI(r - 1, idx);
            eeI(r - 1, idx) = eeI(s, idx);
            eeI(s, idx) = tmp;
          }
        }
      }

      // filter out elements at the same position and find adjacent elements
      // for one matrix
      int idp = 0, del = 0;
      rr(ee(0, idx), idx) = 0;
      dd(ee(0, idx), idx) = 0;

      // for two separate matrices
      for (int r = 0; r < 2; r++) szS(r, idx) = 0;
      for (int r = 0; r < NTAU + 1; r++) mmS(r, idx) = (-kkS(r, idx) + NINT - 1) / NINT;
      int idpS = 0, delS = 0;
      rrS(eeS(0, idx), idx) = 0;
      ddS(eeS(0, idx), idx) = 0;
      szS(mmS(eeS(0, idx), idx), idx) = NDEG + 1;

      // for many separate matrices
      for (int r = 0; r < NMAT; r++) szI(r, idx) = 0;
      for (int r = 0; r < NTAU + 1; r++) mmI(r, idx) = (-kkI(r, idx) + NINT - 1) / NINT;
      int idpI = 0, delI = 0;
      rrI(eeI(0, idx), idx) = 0;
      ddI(eeI(0, idx), idx) = 0;
//      std::cout << " mmI " << mmI(eeI(0, idx),idx) << " eeI " << eeI(0, idx) << " kkI " << kkI(eeI(0, idx),idx) << "\n";
      szI(mmI(eeI(0, idx), idx), idx) = NDEG + 1;

      for (int r = 1; r < NTAU + 1; r++)
      {
        // for one matrix
        if (kk(ee(r - 1, idx), idx) != kk(ee(r, idx), idx)) idp++;
        rr(ee(r, idx), idx) = idp;
        if (kk(ee(r, idx), idx) - kk(ee(r - 1, idx), idx) == 1) del++;
        dd(ee(r, idx), idx) = del;

        // for two separate matrices
        if (mmS(eeS(r - 1, idx), idx) != mmS(eeS(r, idx), idx))
        {
          idpS = 0, delS = 0;
          szS(mmS(eeS(r, idx), idx), idx) = NDEG + 1;
        }
        else
        {
          if (kkS(eeS(r - 1, idx), idx) != kkS(eeS(r, idx), idx))
          {
            idpS++;
            szS(mmS(eeS(r, idx), idx), idx) += NDEG + 1;
          }
          if (kkS(eeS(r, idx), idx) - kkS(eeS(r - 1, idx), idx) == 1)
          {
            delS++;
            szS(mmS(eeS(r, idx), idx), idx) -= 1;
          }
        }
        rrS(eeS(r, idx), idx) = idpS;
        ddS(eeS(r, idx), idx) = delS;

        // for many separate matrices
        if (mmI(eeI(r - 1, idx), idx) != mmI(eeI(r, idx), idx))
        {
          idpI = 0, delI = 0;
          szI(mmI(eeI(r, idx), idx), idx) = NDEG + 1;
        }
        else
        {
          if (kkI(eeI(r - 1, idx), idx) != kkI(eeI(r, idx), idx))
          {
            idpI++;
            szI(mmI(eeI(r, idx), idx), idx) += NDEG + 1;
          }
          if (kkI(eeI(r, idx), idx) - kkI(eeI(r - 1, idx), idx) == 1)
          {
            delI++;
            szI(mmI(eeI(r, idx), idx), idx) -= 1;
          }
        }
        rrI(eeI(r, idx), idx) = idpI;
        ddI(eeI(r, idx), idx) = delI;

      }
    }
  }
  // making periodic
  for (int k = 0; k < NTAU; k++)
  {
    kkMSH(k, NDEG*NINT) = kkMSH(k, 0);
    for (int l = 0; l < NDEG + 1; l++)
    {
      ttMSH(k, l, NDEG*NINT) = ttMSH(k, l, 0);
      ttMSH(NTAU + k, l, NDEG*NINT) = ttMSH(NTAU + k, l, 0);
    }
  }
  // all the delays. the previous was rescaled and we need unscaled
  sys->p_tau(p_tau, time, par);
  sys->p_tau(p_tauMSH, timeMSH, par);
  // evaluate the solution and its derivatives at the collocation points
  Interpolate(solData, sol);
  // diagnostics
#ifdef DEBUG
  count_print();
  count_reset();
#endif
}

void NColloc::Interpolate(Array3D<double>& solData, const Vector& sol)
{
  for (int idx = 0; idx < NDEG*NINT; ++idx)   // i: interval; j: which collocation point
  {
    for (int p = 0; p < NTAU; p++)
    {
      for (int k = 0; k < NDIM; k++)
      {
        solData(k, p, idx) = 0.0;
        for (int l = 0; l < NDEG + 1; l++)
        {
          solData(k, p, idx) += sol(k + NDIM * (l + kk(p + 1, idx) * NDEG)) * tt(p + 1, l, idx);
        }
      }
    }
    // x'(t)
    for (int k = 0; k < NDIM; k++)
    {
      solData(k, NTAU, idx) = 0.0;
      for (int l = 0; l < NDEG + 1; l++)
      {
        solData(k, NTAU, idx) += sol(k + NDIM * (l + kk(0, idx) * NDEG)) * tt(0, l, idx);
      }
    }
    // x'(t-\tau_i)
    for (int p = 0; p < NTAU; p++)
    {
      for (int k = 0; k < NDIM; k++)
      {
        solData(k, NTAU + 1 + p, idx) = 0.0;
        for (int l = 0; l < NDEG + 1; l++)
        {
          solData(k, NTAU + 1 + p, idx) += sol(k + NDIM * (l + kk(p + 1, idx) * NDEG)) * tt(NTAU + 1 + p, l, idx);
        }
      }
    }
  }
}

void NColloc::InterpolateMSH(Array3D<double>& solData, const Vector& sol)
{
  for (int idx = 0; idx < NDEG*NINT+1; ++idx)
  {
    for (int p = 0; p < NTAU; p++)
    {
      for (int k = 0; k < NDIM; k++)
      {
        // std::cout<<"InterR: "<<idx<<", "<<out(idx).Row()<<", "<<out(idx).Col()<<", "<<k<<", "<<p<<"\n";
        solData(k, p, idx)      = 0.0;
        solData(k, NTAU + p, idx) = 0.0;
        for (int l = 0; l < NDEG + 1; l++)
        {
          solData(k, p, idx)      += sol(k + NDIM * (l + kkMSH(p, idx) * NDEG)) * ttMSH(p, l, idx);
          solData(k, NTAU + p, idx) += sol(k + NDIM * (l + kkMSH(p, idx) * NDEG)) * ttMSH(NTAU + p, l, idx);
        }
      }
    }
  }
}

// in complex form on 2*i th places are the reals and on 2*i+1 th places the imaginary parts

void NColloc::InterpolateCPLX(Array3D<double>& solDataRe, Array3D<double>& solDataIm, const Vector& sol)
{
  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int p = 0; p < NTAU; p++)
    {
      for (int k = 0; k < NDIM; k++)
      {
        solDataRe(k, p, idx) = 0.0;
        solDataIm(k, p, idx) = 0.0;
        for (int l = 0; l < NDEG + 1; l++)
        {
          solDataRe(k, p, idx) += sol(2 * (k + NDIM * (l + kk(p + 1, idx) * NDEG))) * tt(p + 1, l, idx);
          solDataIm(k, p, idx) += sol(2 * (k + NDIM * (l + kk(p + 1, idx) * NDEG)) + 1) * tt(p + 1, l, idx);
        }
      }
    }
    for (int k = 0; k < NDIM; k++)
    {
      solDataRe(k, NTAU, idx) = 0.0;
      solDataIm(k, NTAU, idx) = 0.0;
      for (int l = 0; l < NDEG + 1; l++)
      {
        solDataRe(k, NTAU, idx) += sol(2 * (k + NDIM * (l + kk(0, idx) * NDEG))) * tt(0, l, idx);
        solDataIm(k, NTAU, idx) += sol(2 * (k + NDIM * (l + kk(0, idx) * NDEG)) + 1) * tt(0, l, idx);
      }
    }
  }
}

void NColloc::RHS(Vector& rhs, const Vector& par, const Vector& sol)
{
#ifdef DEBUG
  count_RHS++;
#endif
  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    rhs(r) = sol(r + NDIM * NDEG * NINT) - sol(r);
  }

  sys->p_rhs(p_fx, time, solData, par);
  for (int k = 0; k < NDIM; k++)
  {
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      rhs(NDIM + k + NDIM*idx) = par(0) * p_fx(k, idx) - solData(k, NTAU, idx);
    }
  }
}

void NColloc::RHS_p(Vector& rhs, const Vector& par, const Vector& /*sol*/, int alpha)
{
#ifdef DEBUG
  count_RHS_p++;
#endif
  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    rhs(r) = 0.0;
  }

  if (alpha == 0)
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    sys->p_rhs(p_fx, time, solData, par);

    for (int p = 0; p < NDIM; p++)
    {
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        rhs(NDIM + p + NDIM*idx) = -p_fx(p, idx);
      }
    }

    for (int r = 0; r < NTAU; r++)
    {
      int nx=1, vx=r, np=0, vp;

      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
      for (int p = 0; p < NDIM; p++)
      {
        for (int q = 0; q < NDIM; q++)
        {
          for (int idx = 0; idx < NDEG*NINT; ++idx)
          {
            const double d = (p_dtau(r,idx) - p_tau(r,idx) / par(0));
            rhs(NDIM + p + NDIM*idx) += d * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
          }
        }
      }
    }
  }
  else
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    int nx, vx, np, vp;
    nx = 0, np = 1;
    vp = alpha;

    sys->p_deri(p_dfp, time, solData, par, nx, &vx, np, &vp, p_dummy);

    for (int k = 0; k < NDIM; k++)
    {
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        rhs(NDIM + k + NDIM*idx) = -par(0) * p_dfp(k, 0, idx);
      }
    }

    for (int r = 0; r < NTAU; r++)
    {
      nx = 1; np = 0; vx = r;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
      for (int p = 0; p < NDIM; p++)
      {
        for (int q = 0; q < NDIM; q++)
        {
          for (int idx = 0; idx < NDEG*NINT; ++idx)
          {
            rhs(NDIM + p + NDIM*idx) += p_dtau(r, idx) * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
          }
        }
      }
    }
  }
}


void NColloc::RHS_x(SpMatrix& A, const Vector& par, const Vector& /*sol*/)
{
#ifdef DEBUG
  count_RHS_x++;
#endif
  A.Clear('R');

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    A.NewL(2);
    A.WrLi(r, 0) = r;
    A.WrLi(r, 1) = r + NDIM * NDEG * NINT;
    A.WrLx(r, 0) = 1.0;
    A.WrLx(r, 1) = -1.0;
  }

  // MUST PRESERVE THE ORDER. DO NOT CHANGE THE LOOPS
  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      A.NewL(NDIM*((NDEG + 1)*(rr(ee(NTAU,idx), idx) + 1) - dd(ee(NTAU,idx), idx))); // check the line
    }
  }

  for (int k = 0; k < NTAU + 1; k++)
  {
    int nx = 1, vx, np = 0, vp;
    if (k != 0)
    {
      vx = k - 1;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
    }
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
          {
            WRIDX(A, idx, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
            if (k == 0)
            {
              if (p == q) WRDAT(A, idx, p, k, l, q) += tt(0, l, idx);
            }
            else
            {
              WRDAT(A, idx, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
            }
          }
        }
      }
    }
  }

}


//! its very different from all of them
void NColloc::StabJac(StabMatrix& AB, const Vector& par)
{
#ifdef DEBUG
  count_StabJac++;
#endif
  AB.getA0().Clear('R');
  for (int s = 0; s < NMAT; s++) AB.getAI(s).Clear('R');

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    AB.getA0().NewL(1);        // L
    AB.getA0().WrLi(r, 0) = r;
    AB.getA0().WrLx(r, 0) = 1.0;

    for (int s = 0; s < NMAT; s++)
    {
      if (s == 0)
      {
        AB.getAI(s).NewL(1);        // M
        AB.getAI(s).WrLi(r, 0) = r + NDIM * NDEG * NINT;
        AB.getAI(s).WrLx(r, 0) = 1.0;
      }
      else
      {
        AB.getAI(s).NewL(0);
      }
    }
  }

  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      AB.getA0().NewL(NDIM*szI(0, idx));
    }
    for (int s = 0; s < NMAT; s++)
    {
      for (int r = 0; r < NDIM; r++)
      {
        AB.getAI(s).NewL(NDIM*szI(s + 1, idx));
      }
    }
  }

  for (int k = 0; k < NTAU + 1; k++)
  {
    int nx = 1, vx, np = 0, vp;
    if (k != 0)
    {
      vx = k - 1;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
    }
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      for (int l = 0; l < NDEG + 1; l++)
      {
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            if (mmI(k, idx) == 0)  // A matrix -- kkS(eeS(k,idx),idx) >= 0
            {
              WRIDXI(AB.getA0(), idx, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              if (k == 0)
              {
                if (p == q) WRDATI(AB.getA0(), idx, p, k, l, q) += tt(0, l, idx);
              }
              else
              {
                WRDATI(AB.getA0(), idx, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
            }
            else // B matrices
            {
//                 std::cout<<"m"<<mmI(k, idx)<<" ";
              WRIDXI(AB.getAI(mmI(k, idx) - 1), idx, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              if (k == 0)
              {
//          if( p == q ) WRDATSB(B,idx, i,j,p, k,l,q) -= tt(0,l,idx);
              }
              else
              {
                WRDATI(AB.getAI(mmI(k, idx) - 1), idx, p, k, l, q) += par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
            }
          }
        }
      }
    }
  }
}

// similar to RHS_x but with one boundary condition only and multiplied by Z
void NColloc::CharJac_x(SpMatrix& A, const Vector& par, double Z)
{
#ifdef DEBUG
  count_CharJac_x++;
#endif
  A.Clear('R');

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    A.NewL(2);
    A.WrLi(r, 0) = r;
    A.WrLx(r, 0) = 1.0;
    // just for the test function
    A.WrLi(r, 1) = NDIM * NINT * NDEG + r;
    A.WrLx(r, 1) = - 1.0 * Z;
  }

  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      A.NewL(NDIM*((NDEG + 1)*(rr(ee(NTAU, idx), idx) + 1) - dd(ee(NTAU, idx), idx))); // check the line
    }
  }

  for (int k = 0; k < NTAU + 1; k++)
  {
    if (k != 0)
    {
      int nx = 1, vx, np = 0, vp;
      vx = k - 1;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
    }
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;
      const double ZP = pow(Z, zpow);

      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
          {
            WRIDX(A, idx, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
            if (k == 0)             // A
            {
              if (p == q) WRDAT(A, idx, p, k, l, q) += tt(0, l, idx);
            }
            else
            {
              if (zpow > 0)   // -zB0 - z^2*B1 ... - z^n*BN
              {
                WRDAT(A, idx, p, k, l, q) -= ZP * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
              else                         // A
              {
                WRDAT(A, idx, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
            }
          }
        }
      }
    }
  }
}


// this has to be changed only to packed complex.
void NColloc::CharJac_x(SpMatrix& A, const Vector& par, double Re, double Im)
{
#ifdef DEBUG
  count_CharJac_x++;
#endif
  A.Clear('R');

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    A.NewL(3);      // Re
    // L
    A.WrLi(2*r, 0) = 2 * r;
    A.WrLx(2*r, 0) = 1.0;
    // -z M
    A.WrLi(2*r, 1) = 2 * NDIM * NINT * NDEG + 2 * r;
    A.WrLi(2*r, 2) = 2 * NDIM * NINT * NDEG + 2 * r + 1;
    A.WrLx(2*r, 1) = -1.0 * Re;
    A.WrLx(2*r, 2) = 1.0 * Im;
    A.NewL(3);      // Im
    // L
    A.WrLi(2*r + 1, 0) = 2 * r + 1;
    A.WrLx(2*r + 1, 0) = 1.0;
    // -z M
    A.WrLi(2*r + 1, 1) = 2 * NDIM * NINT * NDEG + 2 * r;
    A.WrLi(2*r + 1, 2) = 2 * NDIM * NINT * NDEG + 2 * r + 1;
    A.WrLx(2*r + 1, 1) = -1.0 * Im;
    A.WrLx(2*r + 1, 2) = -1.0 * Re;
  }

  // computing the powers of the multiplier
  Vector ZReP(NMAT + 1), ZImP(NMAT + 1);
  ZReP(0) = 1.0;
  ZImP(0) = 0.0;
  ZReP(1) = Re;
  ZImP(1) = Im;
  for (int r = 2; r < NMAT + 1; r++)
  {
    ZReP(r) = Re * ZReP(r - 1) - Im * ZImP(r - 1);
    ZImP(r) = Re * ZImP(r - 1) + Im * ZReP(r - 1);
  }

  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      A.NewL(2*NDIM*((NDEG + 1)*(rr(ee(NTAU, idx), idx) + 1) - dd(ee(NTAU, idx), idx))); // Re
      A.NewL(2*NDIM*((NDEG + 1)*(rr(ee(NTAU, idx), idx) + 1) - dd(ee(NTAU, idx), idx))); // Im
    }
  }

  for (int k = 0; k < NTAU + 1; k++)
  {
    if (k != 0)
    {
      int nx = 1, vx, np = 0, vp;
      vx = k - 1;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
    }
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;
      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
          {
            WRIDXCPLX(A, idx, p, 0, k, l, q, 0) = 2 * (q + NDIM * (l + NDEG * kk(k, idx)));
            WRIDXCPLX(A, idx, p, 0, k, l, q, 1) = 2 * (q + NDIM * (l + NDEG * kk(k, idx))) + 1;
            WRIDXCPLX(A, idx, p, 1, k, l, q, 0) = 2 * (q + NDIM * (l + NDEG * kk(k, idx)));
            WRIDXCPLX(A, idx, p, 1, k, l, q, 1) = 2 * (q + NDIM * (l + NDEG * kk(k, idx))) + 1;
            if (k == 0)             // A
            {
              if (p == q)
              {
                WRDATCPLX(A, idx, p, 0, k, l, q, 0) += tt(0, l, idx);
                WRDATCPLX(A, idx, p, 1, k, l, q, 1) += tt(0, l, idx);
              }
            }
            else
            { //kkS(ee(k,idx),idx) < 0
              if (zpow > 0)   // -zB0 - z^2*B1 ... - z^n*BN
              {
                WRDATCPLX(A, idx, p, 0, k, l, q, 0) -= ZReP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                WRDATCPLX(A, idx, p, 0, k, l, q, 1) += ZImP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                WRDATCPLX(A, idx, p, 1, k, l, q, 0) -= ZImP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                WRDATCPLX(A, idx, p, 1, k, l, q, 1) -= ZReP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
              else                         // A
              {
                WRDATCPLX(A, idx, p, 0, k, l, q, 0) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                WRDATCPLX(A, idx, p, 1, k, l, q, 1) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
            }
          }
        }
      }
    }
  }
}

// REAL

void NColloc::CharJac_x_p(Vector& V, const Vector& par, const Array3D<double>& phiData, double Z, int alpha)
{
#ifdef DEBUG
  count_CharJac_x_p++;
#endif
  V.Clear();

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    V(r) = 0.0;
  }

  if (alpha == 0)   // a periodusido szerint deriv...
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    for (int k = 0; k < NTAU; k++)
    {
      int nx = 1, vx[2], np = 1, vp = alpha;
      nx = 1;
      np = 0;
      vx[0] = k;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, p_dummy);
      p_fx.Clear(); //OK this is cleared, but why here?
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            p_fx(p, idx) -= p_dfx(p, q, idx) * phiData(q, k, idx);
          }
        }
      }

      nx = 2;
      np = 0;
      for (int r = 0; r < NTAU; r++)
      {
        vx[1] = r;
        vx[0] = k;
        sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, phiData);
        for (int idx = 0; idx < NDEG*NINT; ++idx)
        {
          for (int p = 0; p < NDIM; p++)
          {
            for (int q = 0; q < NDIM; q++)
            {
              p_fx(p, idx) += (p_dtau(r,idx) - p_tau(r,idx) / par(0)) * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
            }
          }
        }
      }
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;
        const double ZP = pow(Z, zpow);

        for (int p = 0; p < NDIM; p++)
        {
          if (zpow > 0)
          {
            V(NDIM + p + NDIM*(idx)) += ZP * p_fx(p,idx);
          }
          else
          {
            V(NDIM + p + NDIM*(idx)) += p_fx(p,idx);
          }
        }
      }
    }
  }
  else
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    for (int k = 0; k < NTAU; k++)
    {
      int nx = 1, vx[2], np = 1, vp = alpha;
      nx = 1;
      np = 1;
      vx[0] = k;
      vp = alpha; // though, this last is never changed...
      sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, p_dummy);
      p_fx.Clear();
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            p_fx(p,idx) -= par(0) * p_dfx(p, q, idx) * phiData(q, k, idx);
          }
        }
      }

      nx = 2;
      np = 0;
      for (int r = 0; r < NTAU; r++)
      {
        vx[1] = r;
        vx[0] = k; // CHANGE THIS to 0, 1
        sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, phiData);
        for (int idx = 0; idx < NDEG*NINT; ++idx)
        {
          for (int p = 0; p < NDIM; p++)
          {
            for (int q = 0; q < NDIM; q++)
            {
              p_fx(p, idx) += p_dtau(r, idx) * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
            }
          }
        }
      }
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;
        const double ZP = pow(Z, zpow);

        for (int p = 0; p < NDIM; p++)
        {
          if (zpow > 0)
          {
            V(NDIM + p + NDIM*(idx)) += ZP * p_fx(p,idx);
          }
          else
          {
            V(NDIM + p + NDIM*(idx)) += p_fx(p,idx);
          }
        }
      }
    }
  }
}

// COMPLEX

void NColloc::CharJac_x_p(Vector& V, const Vector& par,
                          const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm,
                          double Re, double Im, int alpha)
{
#ifdef DEBUG
  count_CharJac_x_p++;
#endif
  V.Clear();

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    V(2*r) = 0.0;
    V(2*r + 1) = 0.0;
  }

  // computing the powers of the multiplier
  Vector ZReP(NMAT + 1), ZImP(NMAT + 1);
  ZReP(0) = 1.0;
  ZImP(0) = 0.0;
  ZReP(1) = Re;
  ZImP(1) = Im;
  for (int r = 2; r < NMAT + 1; r++)
  {
    ZReP(r) = Re * ZReP(r - 1) - Im * ZImP(r - 1);
    ZImP(r) = Re * ZImP(r - 1) + Im * ZReP(r - 1);
  }

  if (alpha == 0)   // a periodusido szerint deriv...
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    for (int k = 0; k < NTAU; k++)
    {
      int nx = 1, vx[2], np = 1, vp = alpha;
      nx = 1;
      np = 0;
      vx[0] = k;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, p_dummy);
      p_fxRe.Clear(); //OK this is cleared, but why here?
      p_fxIm.Clear();
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            p_fxRe(p, idx) -= p_dfx(p, q, idx) * phiDataRe(q, k, idx);
            p_fxIm(p, idx) -= p_dfx(p, q, idx) * phiDataIm(q, k, idx);
          }
        }
      }

      nx = 2;
      np = 0;
      for (int r = 0; r < NTAU; r++)
      {
        vx[1] = r;
        vx[0] = k;
        sys->p_deri(p_dfxRe, time, solData, par, nx, &vx[0], np, &vp, phiDataRe);
        sys->p_deri(p_dfxIm, time, solData, par, nx, &vx[0], np, &vp, phiDataIm);
        for (int idx = 0; idx < NDEG*NINT; ++idx)
        {
          for (int p = 0; p < NDIM; p++)
          {
            for (int q = 0; q < NDIM; q++)
            {
              p_fxRe(p, idx) += (p_dtau(r,idx) - p_tau(r,idx) / par(0)) * p_dfxRe(p, q, idx) * solData(q, NTAU + 1 + r, idx);
              p_fxIm(p, idx) += (p_dtau(r,idx) - p_tau(r,idx) / par(0)) * p_dfxIm(p, q, idx) * solData(q, NTAU + 1 + r, idx);
            }
          }
        }
      }
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

        for (int p = 0; p < NDIM; p++)
        {
          if (zpow > 0)
          {
            V(2*(NDIM + p + NDIM*(idx)))     += ZReP(zpow) * p_fxRe(p,idx) - ZImP(zpow) * p_fxIm(p,idx);
            V(2*(NDIM + p + NDIM*(idx)) + 1) += ZImP(zpow) * p_fxRe(p,idx) + ZReP(zpow) * p_fxIm(p,idx);
          }
          else
          {
            V(2*(NDIM + p + NDIM*(idx)))     += p_fxRe(p,idx);
            V(2*(NDIM + p + NDIM*(idx)) + 1) += p_fxIm(p,idx);
          }
        }
      }
    }
  }
  else
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    for (int k = 0; k < NTAU; k++)
    {
      int nx = 1, vx[2], np = 1, vp = alpha;
      nx = 1;
      np = 1;
      vx[0] = k;
      vp = alpha; // though, this last is never changed...
      sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, p_dummy);
      p_fxRe.Clear();
      p_fxIm.Clear();
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            p_fxRe(p,idx) -= par(0) * p_dfx(p, q, idx) * phiDataRe(q, k, idx);
            p_fxIm(p,idx) -= par(0) * p_dfx(p, q, idx) * phiDataIm(q, k, idx);
          }
        }
      }

      nx = 2;
      np = 0;
      for (int r = 0; r < NTAU; r++)
      {
        vx[1] = r;
        vx[0] = k; // CHANGE THIS to 0, 1
        sys->p_deri(p_dfxRe, time, solData, par, nx, &vx[0], np, &vp, phiDataRe);
        sys->p_deri(p_dfxIm, time, solData, par, nx, &vx[0], np, &vp, phiDataIm);
        for (int idx = 0; idx < NDEG*NINT; ++idx)
        {
          for (int p = 0; p < NDIM; p++)
          {
            for (int q = 0; q < NDIM; q++)
            {
              p_fxRe(p, idx) += p_dtau(r, idx) * p_dfxRe(p, q, idx) * solData(q, NTAU + 1 + r, idx);
              p_fxIm(p, idx) += p_dtau(r, idx) * p_dfxIm(p, q, idx) * solData(q, NTAU + 1 + r, idx);
            }
          }
        }
      }
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

        for (int p = 0; p < NDIM; p++)
        {
          if (zpow > 0)
          {
            V(2*(NDIM + p + NDIM*(idx)))     += ZReP(zpow) * p_fxRe(p,idx) - ZImP(zpow) * p_fxIm(p,idx);
            V(2*(NDIM + p + NDIM*(idx)) + 1) += ZImP(zpow) * p_fxRe(p,idx) + ZReP(zpow) * p_fxIm(p,idx);
          }
          else
          {
            V(2*(NDIM + p + NDIM*(idx)))     += p_fxRe(p,idx);
            V(2*(NDIM + p + NDIM*(idx)) + 1) += p_fxIm(p,idx);
          }
        }
      }
    }
  }
}

void NColloc::CharJac_x_x(SpMatrix& A, const Vector& par, const Array3D<double>& phiData, double Z)
{
#ifdef DEBUG
  count_CharJac_x_x++;
#endif
  Array3D<double> p_t_dfx(NDIM, NDIM, NDEG*NINT);

  A.Clear('R');

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    A.NewL(1);
    A.WrLi(r, 0) = r;
    A.WrLx(r, 0) = 0.0;
  }

  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      A.NewL(NDIM*((NDEG + 1)*(rr(ee(NTAU, idx), idx) + 1) - dd(ee(NTAU, idx), idx))); // check the line
    }
  }

  for (int k = 1; k < NTAU + 1; k++)
  {
    int nx = 2, vx[2], np = 0, vp;
    vx[1] = k - 1;  // CHANGE THIS to 1

    p_dfx.Clear();
    p_t_dfx.Clear();
    for (int r = 0; r < NTAU; r++)
    {
      vx[0] = r;        // CHANGE THIS to 0
      sys->p_deri(p_t_dfx, time, solData, par, nx, &vx[0], np, &vp, phiData);
      //std::cout<<"t_dfx "; t_dfx.Print();
      for (int ra = 0; ra < NDIM; ra++)
      {
        for (int rb = 0; rb < NDIM; rb++)
        {
          for (int idx = 0; idx < NDEG*NINT; ++idx)
          {
            p_dfx(ra, rb, idx) += p_t_dfx(ra, rb, idx);
          }
        }
      }
    }

    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;
      const double ZP = pow(Z, zpow);

      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
          {
            WRIDX(A, idx, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
            // This is valid only for k != 0
            if (zpow > 0)   // -zB
            {
              WRDAT(A, idx, p, k, l, q) -= ZP * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
            }
            else                         // A
            {
              WRDAT(A, idx, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
            }
          }
        }
      }
    }
  }
}

void NColloc::CharJac_x_x(SpMatrix& A, const Vector& par,
                          const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm,
                          double Re, double Im)
{
#ifdef DEBUG
  count_CharJac_x_x++;
#endif
  Array3D<double> p_t_dfxRe(NDIM, NDIM, NDEG*NINT);
  Array3D<double> p_t_dfxIm(NDIM, NDIM, NDEG*NINT);

  A.Clear('R');

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    A.NewL(1);      // Re
    A.WrLi(2*r, 0) = 2 * r;
    A.WrLx(2*r, 0) = 0.0;
    A.NewL(1);      // Im
    A.WrLi(2*r + 1, 0) = 2 * r + 1;
    A.WrLx(2*r + 1, 0) = 0.0;
  }

  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      A.NewL(NDIM*((NDEG + 1)*(rr(ee(NTAU, idx), idx) + 1) - dd(ee(NTAU, idx), idx))); // Re
      A.NewL(NDIM*((NDEG + 1)*(rr(ee(NTAU, idx), idx) + 1) - dd(ee(NTAU, idx), idx))); // Im
    }
  }

  // computing the powers of the multiplier
  Vector ZReP(NMAT + 1), ZImP(NMAT + 1);
  ZReP(0) = 1.0;
  ZImP(0) = 0.0;
  ZReP(1) = Re;
  ZImP(1) = Im;
  for (int r = 2; r < NMAT + 1; r++)
  {
    ZReP(r) = Re * ZReP(r - 1) - Im * ZImP(r - 1);
    ZImP(r) = Re * ZImP(r - 1) + Im * ZReP(r - 1);
  }

  for (int k = 0; k < NTAU + 1; k++)
  {
    int nx = 2, vx[2], np = 0, vp;
    if (k != 0)
    {
      vx[1] = k - 1;  // CHANGE THIS to 1

      p_dfxRe.Clear();
      p_t_dfxRe.Clear();
      p_dfxIm.Clear();
      p_t_dfxIm.Clear();
      for (int r = 0; r < NTAU; r++)
      {
        vx[0] = r;        // CHANGE THIS to 0
        sys->p_deri(p_t_dfxRe, time, solData, par, nx, &vx[0], np, &vp, phiDataRe);
        sys->p_deri(p_t_dfxIm, time, solData, par, nx, &vx[0], np, &vp, phiDataIm);
        //std::cout<<"t_dfx "; t_dfx.Print();
        for (int ra = 0; ra < NDIM; ra++)
        {
          for (int rb = 0; rb < NDIM; rb++)
          {
            for (int idx = 0; idx < NDEG*NINT; ++idx)
            {
              p_dfxRe(ra, rb, idx) += p_t_dfxRe(ra, rb, idx);
              p_dfxIm(ra, rb, idx) += p_t_dfxIm(ra, rb, idx);
            }
          }
        }
      }
      //std::cout<<"dfx "; dfx.Print();
    }
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;

      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
          {
            WRIDXCPLXM(A, idx, p, 0, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
            WRIDXCPLXM(A, idx, p, 1, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
            if (k == 0)             // A
            {
              // WRDATCPLX(A,idx, i,j,p,0, k,l,q,0) += tt(0,l,idx); // no derivative !!!!!!
              // WRDATCPLX(A,idx, i,j,p,1, k,l,q,1) += tt(0,l,idx);
            }
            else
            {
              if (zpow > 0)   // -zB
              {
                WRDATCPLXM(A, idx, p, 0, k, l, q) -=
                  ZReP(zpow) * par(0) * p_dfxRe(p, q, idx) * tt(k, l, idx) - ZImP(zpow) * par(0) * p_dfxIm(p, q, idx) * tt(k, l, idx);
                WRDATCPLXM(A, idx, p, 1, k, l, q) -=
                  ZImP(zpow) * par(0) * p_dfxRe(p, q, idx) * tt(k, l, idx) + ZReP(zpow) * par(0) * p_dfxIm(p, q, idx) * tt(k, l, idx);
              }
              else                         // A
              {
                WRDATCPLXM(A, idx, p, 0, k, l, q) -= par(0) * p_dfxRe(p, q, idx) * tt(k, l, idx);
                WRDATCPLXM(A, idx, p, 1, k, l, q) -= par(0) * p_dfxIm(p, q, idx) * tt(k, l, idx);
              }
            }
          }
        }
      }
    }
  }
}

// this is for CharmatCPLX

void NColloc::CharJac_x_z(Vector& V, const Vector& par, const Vector& phi,
                          const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm, double Re, double Im)
{
#ifdef DEBUG
  count_CharJac_x_z++;
#endif
  V.Clear();

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    V(2*r)   = -phi(2 * (r + NDIM * NDEG * NINT));
    V(2*r + 1) = -phi(2 * (r + NDIM * NDEG * NINT) + 1);
  }

  // computing the powers of the multiplier
  Vector ZReP(NMAT + 1), ZImP(NMAT + 1);
  ZReP(0) = 1.0;
  ZImP(0) = 0.0;
  ZReP(1) = Re;
  ZImP(1) = Im;
  for (int r = 2; r < NMAT + 1; r++)
  {
    ZReP(r) = Re * ZReP(r - 1) - Im * ZImP(r - 1);
    ZImP(r) = Re * ZImP(r - 1) + Im * ZReP(r - 1);
  }

  for (int k = 0; k < NTAU; k++)
  {
    int nx = 1, vx, np = 0, vp;
    vx = k;
    sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);

    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

      for (int p = 0; p < NDIM; p++)
      {
        for (int q = 0; q < NDIM; q++)
        {
          if (zpow > 0)
          {
            V(2*(NDIM + p + NDIM*(idx)))     -= zpow * (ZReP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataRe(q, k, idx) -
                                                                ZImP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataIm(q, k, idx));
            V(2*(NDIM + p + NDIM*(idx)) + 1) -= zpow * (ZImP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataRe(q, k, idx) +
                                                                ZReP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataIm(q, k, idx));
          }
        }
      }
    }
  }
}

//!
//! from now CharmatLPAUT
//!

void NColloc::CharJac_mB(SpMatrix& B, const Vector& par, double Z)
{
#ifdef DEBUG
  count_CharJac_mB++;
#endif
  B.Clear('R');

  Vector ZP(NMAT + 1);
  ZP(0) = 1.0;
  for (int r = 1; r < NMAT + 1; r++)
  {
    ZP(r) = Z * ZP(r - 1);
  }

  // boundary conditions: no boundary condition
  for (int r = 0; r < NDIM; r++)
  {
    B.NewL(1);        // M
    B.WrLi(r, 0) = r + NDIM * NDEG * NINT;
    B.WrLx(r, 0) = -1.0;
  }

  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      B.NewL(NDIM*szS(1, idx));
    }
  }

  for (int k = 1; k < NTAU + 1; k++)
  {
    int nx = 1, vx, np = 0, vp;
    vx = k - 1;
    sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;

      for (int l = 0; l < NDEG + 1; l++)
      {
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            if (zpow > 0)  // B matrix -- kkS(eeS(k,idx),idx) < 0
            {
              WRIDXS(B, idx, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              WRDATS(B, idx, p, k, l, q) -= zpow * ZP(zpow - 1) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
            }
          }
        }
      }
    }
  }
}

// same as CharJac_x_p, but only writes the B part

void NColloc::CharJac_mB_p(Vector& V, const Vector& par, const Array3D<double>& phiData, double Z, int alpha)
{
#ifdef DEBUG
  count_CharJac_mB_p++;
#endif
  Matrix t_dfx(NDIM, NDIM);

  V.Clear();

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    V(r) = 0.0;
  }

  Vector ZP(NMAT + 1);
  ZP(0) = 1.0;
  for (int r = 1; r < NMAT + 1; r++)
  {
    ZP(r) = Z * ZP(r - 1);
  }

  if (alpha == 0)   // a periodusido szerint deriv...
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    for (int k = 0; k < NTAU; k++)
    {
      // START verbatim from CharJac_x_p
      int nx = 1, vx[2], np = 1, vp = alpha;
      nx = 1;
      np = 0;
      vx[0] = k;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, p_dummy);
      p_fx.Clear(); //OK this is cleared, but why here?
      for (int p = 0; p < NDIM; p++)
      {
        for (int q = 0; q < NDIM; q++)
        {
          for (int idx = 0; idx < NDEG*NINT; ++idx)
          {
            p_fx(p, idx) -= p_dfx(p, q, idx) * phiData(q, k, idx);
          }
        }
      }

      nx = 2;
      np = 0;
      for (int r = 0; r < NTAU; r++)
      {
        vx[1] = r;
        vx[0] = k;
        sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, phiData);
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            for (int idx = 0; idx < NDEG*NINT; ++idx)
            {
              p_fx(p, idx) += (p_dtau(r,idx) - p_tau(r,idx) / par(0)) * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
            }
          }
        }
      }
      // END verbatim from CharJac_x_p
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

        for (int p = 0; p < NDIM; p++)
        {
          if (zpow > 0)
          {
            V(NDIM + p + NDIM*(idx)) += zpow * ZP(zpow - 1) * p_fx(p, idx);
          }
        }
      }
    }
  }
  else
  {
    sys->p_dtau(p_dtau, time, par, alpha);
    for (int k = 0; k < NTAU; k++)
    {
      // START verbatim from CharJac_x_p
      int nx = 1, vx[2], np = 1, vp = alpha;
      nx = 1;
      np = 1;
      vx[0] = k;
      vp = alpha; // though, this last is never changed...
      sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, p_dummy);
      p_fx.Clear();
      for (int p = 0; p < NDIM; p++)
      {
        for (int q = 0; q < NDIM; q++)
        {
          for (int idx = 0; idx < NDEG*NINT; ++idx)
          {
            p_fx(p,idx) -= par(0) * p_dfx(p, q, idx) * phiData(q, k, idx);
          }
        }
      }

      nx = 2;
      np = 0;
      for (int r = 0; r < NTAU; r++)
      {
        vx[1] = r;
        vx[0] = k; // CHANGE THIS to 0, 1
        sys->p_deri(p_dfx, time, solData, par, nx, &vx[0], np, &vp, phiData);
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            for (int idx = 0; idx < NDEG*NINT; ++idx)
            {
              p_fx(p, idx) += p_dtau(r, idx) * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
            }
          }
        }
      }
      // END verbatim from CharJac_x_p
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

        for (int p = 0; p < NDIM; p++)
        {
          if (zpow > 0)
          {
            V(NDIM + p + NDIM*(idx)) += zpow * ZP(zpow - 1) * p_fx(p, idx);
          }
        }
      }
    }
  }
}

// like x_x but write bpart only
void NColloc::CharJac_mB_x(SpMatrix& B, const Vector& par, const Array3D<double>& phiData, double Z)
{
#ifdef DEBUG
  count_CharJac_mB_x++;
#endif
  Array3D<double> p_t_dfx(NDIM, NDIM, NDEG*NINT);

  B.Clear('R');

  // boundary conditions: no boundary condition
  for (int r = 0; r < NDIM; r++)
  {
    B.NewL(0);        // M
//   B.WrLi(r,0) = r+NDIM*NDEG*NINT;
//   B.WrLx(r,0) = 0.0;
  }

  for (int idx = 0; idx < NDEG*NINT; ++idx)
  {
    for (int r = 0; r < NDIM; r++)
    {
      B.NewL(NDIM*szS(1, idx));
    }
  }

  Vector ZP(NMAT + 1);
  ZP(0) = 1.0;
  for (int r = 1; r < NMAT + 1; r++)
  {
    ZP(r) = Z * ZP(r - 1);
  }

  for (int k = 1; k < NTAU + 1; k++)
  {
    // START verbatim from CharJac_x_x
    int nx = 2, vx[2], np = 0, vp;
    vx[1] = k - 1;  // CHANGE THIS to 1

    p_dfx.Clear();
    p_t_dfx.Clear();
    for (int r = 0; r < NTAU; r++)
    {
      vx[0] = r;        // CHANGE THIS to 0
      sys->p_deri(p_t_dfx, time, solData, par, nx, &vx[0], np, &vp, phiData);
      //std::cout<<"t_dfx "; t_dfx.Print();
      for (int ra = 0; ra < NDIM; ra++)
      {
        for (int rb = 0; rb < NDIM; rb++)
        {
          for (int idx = 0; idx < NDEG*NINT; ++idx)
          {
            p_dfx(ra, rb, idx) += p_t_dfx(ra, rb, idx);
          }
        }
      }
    }
    // END verbatim from CharJac_x_x
    for (int idx = 0; idx < NDEG*NINT; ++idx)
    {
      const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;

      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
          {
            if (zpow > 0)  // A matrix -- kkS(eeS(k,idx),idx) >= 0
            {
              WRIDXS(B, idx, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              WRDATS(B, idx, p, k, l, q) -= zpow * ZP(zpow - 1) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
            }
          }
        }
      }
    }
  }
}

//!
//! until this CharmatLPAUT
//!

//!
//! this is also for CharmatLPAUT: for computing q_0 and its derivatives: q_0, D_x q_0, D_p q_0
//!

void NColloc::CharJac_MSHphi(Vector& V, const Vector& par, const Array3D<double>& solMSHData)
{
#ifdef DEBUG
  count_CharJac_MSHphi++;
#endif
  sys->p_rhs(p_fxMSH, timeMSH, solData, par);

  for (int idx = 0; idx < NDEG*NINT+1; ++idx)   // i: interval; j: which collocation point
  {
    for (int k = 0; k < NDIM; k++) V(k + NDIM*idx) = - p_fxMSH(k, idx);
  }
}

void NColloc::CharJac_MSHphi_p(Vector& V, const Vector& par, const Array3D<double>& solMSHData, int alpha)
{
#ifdef DEBUG
  count_CharJac_MSHphi_p++;
#endif
  V.Clear(); /// it is not cleared otherwise!!!!

  // boundary conditions
  if (alpha == 0)
  {
    sys->p_dtau(p_dtauMSH, timeMSH, par, alpha);
    for (int r = 0; r < NTAU; r++)
    {
      int nx, vx, np, vp;
      nx = 1;
      np = 0;
      vx = r;
      sys->p_deri(p_dfxMSH, timeMSH, solData, par, nx, &vx, np, &vp, p_dummy);

      for (int p = 0; p < NDIM; p++)
      {
        for (int q = 0; q < NDIM; q++)
        {
          for (int idx = 0; idx < NDEG*NINT+1; ++idx)   // i: interval; j: which collocation point
          {
            V(p + NDIM*idx) += ((p_dtauMSH(r,idx) - p_tauMSH(r,idx) / par(0)) / par(0)) * p_dfxMSH(p, q, idx) * solData(q, NTAU + r, idx);
          }
        }
      }
    }
  }
  else
  {
    sys->p_dtau(p_dtauMSH, timeMSH, par, alpha);
    for (int r = 0; r < NTAU; r++)
    {
      int nx, vx, np, vp;
      nx = 0, np = 1;
      vp = alpha;

      sys->p_deri(p_dfpMSH, timeMSH, solData, par, nx, &vx, np, &vp, p_dummy);

      nx = 1, np = 0;
      vx = r;
      sys->p_deri(p_dfxMSH, timeMSH, solData, par, nx, &vx, np, &vp, p_dummy);

      for (int k = 0; k < NDIM; k++)
      {
        for (int idx = 0; idx < NDEG*NINT+1; ++idx)   // i: interval; j: which collocation point
        {
          V(k + NDIM*idx) = - p_dfpMSH(k,0,idx);
        }
      }
      for (int p = 0; p < NDIM; p++)
      {
        for (int q = 0; q < NDIM; q++)
        {
          for (int idx = 0; idx < NDEG*NINT+1; ++idx)   // i: interval; j: which collocation point
          {
            V(p + NDIM*idx) += (p_dtauMSH(r,idx) / par(0)) * p_dfxMSH(p, q, idx) * solData(q, NTAU + r, idx);
          }
        }
      }
    }
  }
}


//!
//! until this CharmatLPAUT
//!
