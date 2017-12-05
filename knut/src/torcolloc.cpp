// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "exprsystem.h"
#include "matrix.h"
#include "spmatrix.h"
#include "hypermatrix.h"
#include "plot.h"
#include "polynomial.h"
#include <cmath>
#include <fstream>
#include <functional>
#include <algorithm>

struct comp : public std::binary_function<size_t, size_t, bool>
{
  const size_t* V;
  comp(size_t *V_) : V(V_) {}
  bool operator()(size_t i, size_t j)
  {
    return V[i] < V[j];
  }
};

// kk[.]     : i1,i2,j1,j2 -> idx2
// kk[ee[.]] :      s      -> idx2
// ee[.]     : i1,i2,j1,j2 ->  s
// rr[.]     : i1,i2,j1,j2 ->  index of Ai

//-----------------------------------------------------------------------------------------------
//
//    KNDdeTorusCollocation CLASS
//
//-----------------------------------------------------------------------------------------------

#define NTAU ntau
#define NDIM ndim
#define NPAR npar

#define RHO VarToIndex(VarRot,NPAR)

#include "torcolloc.h"
#include "pointtype.h"

KNDdeTorusCollocation::KNDdeTorusCollocation(KNExprSystem& sys_, size_t ndeg1_, size_t ndeg2_, size_t nint1_, size_t nint2_) :
    sys(&sys_),
    ndim(sys_.ndim()), ntau(sys_.ntau()), npar(sys_.npar()),
    ndeg1(ndeg1_), ndeg2(ndeg2_), nint1(nint1_), nint2(nint2_),
    col1(ndeg1_),    col2(ndeg2_),
    mesh1(ndeg1_ + 1), mesh2(ndeg2_ + 1),
    lgr1(ndeg1_ + 1, ndeg1_ + 1), lgr2(ndeg2_ + 1, ndeg2_ + 1),
    dlg1(ndeg1_ + 1, ndeg1_ + 1), dlg2(ndeg2_ + 1, ndeg2_ + 1),
    I1((ndeg1_ + 1), (ndeg1_ + 1)),
    ID1((ndeg1_ + 1), (ndeg1_ + 1)),
    I2((ndeg2_ + 1), (ndeg2_ + 1)),
    ID2((ndeg2_ + 1), (ndeg2_ + 1)),
    mlg1((ndeg1_ + 1)*(ndeg1_ + 1)),
    mlg2((ndeg2_ + 1)*(ndeg2_ + 1)),
    mlgd1((ndeg1_ + 1)*(ndeg1_ + 1)),
    mlgd2((ndeg2_ + 1)*(ndeg2_ + 1)),
    ilg1((ndeg1_ + 1)*(ndeg1_ + 1) + 1),
    ilg2((ndeg2_ + 1)*(ndeg2_ + 1) + 1),
    ilgd1((ndeg1_ + 1)*(ndeg1_ + 1) + 1),
    ilgd2((ndeg2_ + 1)*(ndeg2_ + 1) + 1),
    time1(ndeg1*ndeg2*nint1*nint2), time2(ndeg1*ndeg2*nint1*nint2),
    kk((ntau+1)*(ndeg1+1)*(ndeg2+1), ndeg1*ndeg2*nint1*nint2),
    ee((ntau+1)*(ndeg1+1)*(ndeg2+1), ndeg1*ndeg2*nint1*nint2),
    rr((ntau+1)*(ndeg1+1)*(ndeg2+1), ndeg1*ndeg2*nint1*nint2),
    p_tau(ntau, ndeg1*ndeg2*nint1*nint2), p_dtau(ntau, ndeg1*ndeg2*nint1*nint2),
    p_xx(ndim, ntau+2*(ntau+1), ndeg1*ndeg2*nint1*nint2),
    p_fx(ndim, ndeg1*ndeg2*nint1*nint2),
    p_dfp(ndim, 1, ndeg1*ndeg2*nint1*nint2),
    p_dfx(ndim, ndim, ndeg1*ndeg2*nint1*nint2),
    p_dummy(0, 0, ndeg1*ndeg2*nint1*nint2)
{
  lobatto(mesh1);
  lobatto(mesh2);
  gauss(col1);
  gauss(col2);
  for (size_t i = 0; i < mesh1.size(); i++)
  {
    poly_coeff_lgr(lgr1(i), mesh1, i);
    poly_coeff_diff(dlg1(i), lgr1(i));
  }
  for (size_t i = 0; i < mesh2.size(); i++)
  {
    poly_coeff_lgr(lgr2(i), mesh2, i);
    poly_coeff_diff(dlg2(i), lgr2(i));
  }

  // here comes the phase condition for par(OMEGA0) and par(OMEGA1)
  // the integration in the bottom border
  // construct the diffint matrix
  for (size_t i = 0; i < ndeg1 + 1; i++)
  {
    for (size_t k = 0; k < ndeg1 + 1; k++)
    {
      mlg1.clear();
      mlgd1.clear();
      ilg1.clear();
      ilgd1.clear();
      poly_coeff_mul(mlg1,  lgr1(i), lgr1(k));
      poly_coeff_mul(mlgd1, dlg1(i), lgr1(k));
      poly_coeff_int(ilg1,  mlg1);
      poly_coeff_int(ilgd1,  mlgd1);
      I1(i, k) = (poly_eval(ilg1, 1.0)-poly_eval(ilg1, 0.0)) / nint1;
      ID1(i, k) = poly_eval(ilgd1, 1.0);
    }
  }
  for (size_t j = 0; j < ndeg2 + 1; j++)
  {
    for (size_t l = 0; l < ndeg2 + 1; l++)
    {
      mlg2.clear();
      mlgd2.clear();
      ilg2.clear();
      ilgd2.clear();
      poly_coeff_mul(mlg2,  lgr2(j), lgr2(l));
      poly_coeff_mul(mlgd2, dlg2(j), lgr2(l));
      poly_coeff_int(ilg2,  mlg2);
      poly_coeff_int(ilgd2,  mlgd2);
      I2(j, l) = (poly_eval(ilg2, 1.0) - poly_eval(ilg2, 0.0))/ nint2;
      ID2(j, l) = poly_eval(ilgd2, 1.0);
    }
  }
}

void KNDdeTorusCollocation::init(const KNVector& sol, const KNVector& par)
{
  std::vector<double> t1(NTAU+1), t2(NTAU+1);
  
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1; j1++)
        {
          const size_t idx = idxmap(j1, j2, i1, i2); // this is a unique

          time1(idx) = ((double)i1 + col1(j1)) / nint1;
          time2(idx) = ((double)i2 + col2(j2)) / nint2;
        }
      }
    }
  }

  sys->p_tau(p_tau, time1, par);

  p_xx.clear();
  for (size_t idx = 0; idx < time1.size(); ++idx)
  {
    t1[0] = time1(idx);
    t2[0] = time2(idx);
    for (size_t k = 0; k < NTAU; k++)
    {
      t1[1+k] = (t1[0] - p_tau(k,idx) / par(0)) - floor(t1[0] - p_tau(k,idx) / par(0));
      t2[1+k] = (t2[0] - par(RHO) * p_tau(k,idx) / par(0)) - floor(t2[0] - par(RHO) * p_tau(k,idx) / par(0));
    }

    for (size_t k = 0; k < NTAU + 1; k++)
    {
      auto i1 = static_cast<size_t>(floor(nint1 * t1[k]));
      auto i2 = static_cast<size_t>(floor(nint2 * t2[k]));
      // some error checking
      if ((t1[k] > 1.0) || (t2[k] > 1.0) ||
          (t1[k] < 0.0) || (t2[k] < 0.0)) std::cout << "Er ";
      if ((i1 >= nint1) || (i2 >= nint2)) std::cout << "Ei ";
      // end error checking
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          const size_t idxKK = idxkk(j1, j2, k);
          kk(idxKK,idx) = idxmap(j1, j2, i1, i2);
          ee(idxKK,idx) = idxKK;
        }
      }
    }

    // sorting
    comp aa(kk.pointer(0,idx));
    std::sort(ee.pointer(0,idx), ee.pointer(0,idx) + (NTAU + 1)*(ndeg1 + 1)*(ndeg2 + 1), aa);

    // filtering same indices
    rr(ee(0,idx),idx) = 0;
    size_t nz = 0;
    for (size_t i = 1; i < (NTAU + 1)*(ndeg1 + 1)*(ndeg2 + 1); i++)
    {
      if (kk(ee(i-1,idx),idx) != kk(ee(i,idx),idx)) nz++;
      rr(ee(i,idx),idx) = nz;
    }

    // interpolate the solution
    for (size_t k = 1; k < NTAU + 1; k++)
    {
      const double c1 = nint1 * t1[k] - floor(nint1 * t1[k]);
      const double c2 = nint2 * t2[k] - floor(nint2 * t2[k]);
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          const size_t idxKK = idxkk(j1, j2, k);
          const double cf = poly_lgr_eval(mesh1, j1, c1) * poly_lgr_eval(mesh2, j2, c2);
          for (size_t p = 0; p < NDIM; p++)
          {
            p_xx(p, k - 1, idx) += cf * sol(p + NDIM * kk(idxKK, idx));
          }
        }
      }
    }
    // derivatives at all delays
    for (size_t k = 0; k < NTAU + 1; k++)
    {
      const double c1 = nint1 * t1[k] - floor(nint1 * t1[k]);
      const double c2 = nint2 * t2[k] - floor(nint2 * t2[k]);
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          const size_t idxKK = idxkk(j1, j2, k);
          const double cf1 = poly_dlg_eval(mesh1, j1, c1) * nint1 * poly_lgr_eval(mesh2, j2, c2);
          const double cf2 = poly_lgr_eval(mesh1, j1, c1) * poly_dlg_eval(mesh2, j2, c2) * nint2;
          for (size_t p = 0; p < NDIM; p++)
          {
            p_xx(p, NTAU + 2*k, idx)   += cf1 * sol(p + NDIM * kk(idxKK, idx));
            p_xx(p, NTAU + 2*k + 1, idx) += cf2 * sol(p + NDIM * kk(idxKK, idx));
          }
        }
      }
    }
  }
//   for (int idx = 0; idx < time1.size(); ++idx) std::cout<<p_xx(0, NTAU, idx)<<"\t";
//   std::cout<<"\np_xx(0,...) was\n";
}

//***************************************************************************
//
// NOT Just for periodic equations
//
//  D[u[x,y],x] + \rho*D[u[x,y],y] = T * f[t, u(x-tau(T)/T, y-rho*tau(T)/T), ... ]
//
//  ParRot is RHO
//
//***************************************************************************

void KNDdeTorusCollocation::jacobian(KNSparseMatrix& A, KNArray1D< KNVector* > Avar, KNVector& rhs, KNVector& par, KNVector& /*sol*/, KNArray1D<size_t>& var)
{
  A.clear();
  rhs.clear();
  for (size_t r = 0; r < var.size(); r++) Avar(r)->clear();
  // rhs doesn't need to be cleared

  // creates kk, ee, rr & interpolates p_xx & gets p_tau
//   init(sol, par);
  // builds up the structure of the sparse matrix
  for (size_t idx = 0; idx < time1.size(); ++idx)
  {
    const size_t lend = NDIM * (rr(ee((NTAU+1)*(ndeg1+1)*(ndeg2+1)-1,idx),idx) + 1);
    for (size_t p = 0; p < NDIM; p++) A.newLine(lend);
  }

  // the right-hand side
  sys->p_rhs(p_fx, time1, p_xx, par, 0);
  for (size_t idx = 0; idx < time1.size(); ++idx)
  {
    for (size_t p = 0; p < NDIM; p++) 
    {
      rhs(p + NDIM*idx) = par(0) * p_fx(p,idx) - p_xx(p, NTAU, idx) - par(RHO) * p_xx(p, NTAU + 1, idx);
    }
  }

  // derivatives w.r.t the parameters
  for (size_t r = 0; r < var.size(); r++)
  {
    KNVector& deri = *Avar(r);
    deri.clear();
//!!!!!!!!!!!!!!!!!
//!BEGIN w.r.t the period
//!!!!!!!!!!!!!!!!!
    if (var(r) == 0)
    {
      for (size_t idx = 0; idx < time1.size(); ++idx)
      {
        for (size_t p = 0; p < NDIM; p++) deri(p + NDIM*idx) = -p_fx(p,idx);
      }
      sys->p_dtau(p_dtau, time1, par, 0);
      for (size_t k = 0; k < NTAU; k++)
      {
        size_t nx = 1, vx = k, np = 0, vp = 0;
        sys->p_deri(p_dfx, time1, p_xx, par, 0, nx, &vx, np, &vp, p_dummy);
        for (size_t idx = 0; idx < time1.size(); ++idx)
        {
          for (size_t p = 0; p < NDIM; p++)
          {
            for (size_t q = 0; q < NDIM; q++)
            {
              const double d = p_tau(k,idx) / par(0) - p_dtau(k,idx);
              deri(p + NDIM*idx) -= p_dfx(p, q, idx) * d * (p_xx(q, NTAU + 2 * (k + 1), idx) + par(RHO) * p_xx(q, NTAU + 2 * (k + 1) + 1, idx));
            }
          }
        }
      }
    }
    else
//!!!!!!!!!!!!!!!!!
//!END w.r.t the period
//!!!!!!!!!!!!!!!!!
//
//!!!!!!!!!!!!!!!!!
//!BEGIN w.r.t the ordinary parameters
//!!!!!!!!!!!!!!!!!
    // derivatives w.r.t. the real parameters
    if (var(r) < NPAR)
    {
      size_t nx = 0, vx = 0, np = 1, vp = var(r);
      sys->p_deri(p_dfp, time1, p_xx, par, 0, nx, &vx, np, &vp, p_dummy);
      for (size_t idx = 0; idx < time1.size(); ++idx)
      {
        for (size_t p = 0; p < NDIM; p++) deri(p + NDIM*idx) = - par(0) * p_dfp(p, 0, idx);
      }
    }
    else
//!!!!!!!!!!!!!!!!!
//!END w.r.t the ordinary parameters
//!!!!!!!!!!!!!!!!!
//
//!!!!!!!!!!!!!!!!!
//!BEGIN w.r.t the rotation number
//!!!!!!!!!!!!!!!!!
    if (var(r) == RHO)
    {
      for (size_t idx = 0; idx < time1.size(); ++idx)
      {
        for (size_t p = 0; p < NDIM; p++) deri(p + NDIM*idx) = p_xx(p, NTAU + 1, idx);
      }
      sys->p_dtau(p_dtau, time1, par, 0);
      for (size_t k = 0; k < NTAU; k++)
      {
        const size_t nx = 1, vx = k, np = 0, vp = 0;
        sys->p_deri(p_dfx, time1, p_xx, par, 0, nx, &vx, np, &vp, p_dummy);
        for (size_t idx = 0; idx < time1.size(); ++idx)
        {
          const double d = -p_tau(k, idx);
          for (size_t p = 0; p < NDIM; p++)
          {
            for (size_t q = 0; q < NDIM; q++)
            {
              deri(p + NDIM*idx) -= p_dfx(p, q, idx) * d * p_xx(q, NTAU + 2 * (k + 1) + 1, idx);
            }
          }
        }
      }
    }
//!!!!!!!!!!!!!!!!!
//!END w.r.t the rotation number
//!!!!!!!!!!!!!!!!!
    else std::cout << "Jac:NNN\n";
  }
  // fill the matrix
  for (size_t k = 0; k < NTAU + 1; k++)
  {
    if (k != 0)
    {
      // evaluate the solution
      const size_t nx = 1, vx = k - 1, np = 0, vp = 0;
      sys->p_deri(p_dfx, time1, p_xx, par, 0, nx, &vx, np, &vp, p_dummy);
    }

    for (size_t idx = 0; idx < time1.size(); ++idx)
    {
      double tk1, tk2;
      if (k != 0)
      {
        tk1 = time1(idx) - p_tau(k-1,idx) / par(0);
        tk2 = time2(idx) - par(RHO) * p_tau(k-1,idx) / par(0);
      } else
      {
        tk1 = time1(idx);
        tk2 = time2(idx);
      }
      const double c1 = nint1 * tk1 - floor(nint1 * tk1);
      const double c2 = nint2 * tk2 - floor(nint2 * tk2);

      // the real jacobian
      for (size_t l2 = 0; l2 < ndeg2 + 1; l2++)
      {
        for (size_t l1 = 0; l1 < ndeg1 + 1; l1++)
        {
          const size_t idxK = idxkk(l1, l2, k);
          if (k != 0)
          {
            const double cf = poly_lgr_eval(mesh1, l1, c1) * poly_lgr_eval(mesh2, l2, c2);
            for (size_t p = 0; p < NDIM; p++)
            {
              for (size_t q = 0; q < NDIM; q++)
              {
                // jacobian of rhs
                A.writeIndex(p + NDIM*idx, q + NDIM*rr(idxK, idx)) = static_cast<int>(q + NDIM * kk(idxK, idx));
                A.writeData(p + NDIM*idx, q + NDIM*rr(idxK, idx)) -= cf * par(0) * p_dfx(p, q, idx);
              }
            }
          }
          else
          {
            const double cf1 = poly_dlg_eval(mesh1, l1, c1) * nint1 * poly_lgr_eval(mesh2, l2, c2);
            const double cf2 = poly_lgr_eval(mesh1, l1, c1) * poly_dlg_eval(mesh2, l2, c2) * nint2;
            for (size_t p = 0; p < NDIM; p++)
            {
              for (size_t q = 0; q < NDIM; q++)
              {
                // derivative part of the jacobian
                A.writeIndex(p + NDIM*idx, q + NDIM*rr(idxK, idx)) = static_cast<int>(q + NDIM * kk(idxK, idx));
                if (p == q)
                  A.writeData(p + NDIM*idx, q + NDIM*rr(idxK, idx)) += cf1 + par(RHO) * cf2;
              }
            }
          }
          // error check
          if (kk(idxK, idx) > (ndeg1*nint1)*(ndeg2*nint2)) std::cout << "D" << kk(idxK, idx);
        }
      }
    }
  }
}

//***************************************************************************************
//
// End of the periodic case
//
//***************************************************************************************

void KNDdeTorusCollocation::PhaseONE(KNVector& ph, KNVector& presol)
{
  ph.clear();
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          size_t idx1 = idxmap(j1, j2, i1, i2);
          // matrix multiplication
          for (size_t l2 = 0; l2 < ndeg2 + 1; l2++)
          {
            for (size_t l1 = 0; l1 < ndeg1 + 1; l1++)
            {
              const size_t idx2 = idxmap(l1, l2, i1, i2);
              for (size_t p = 0; p < NDIM; p++)
              {
                ph(p + NDIM*idx1) += presol(p + NDIM * idx2) * I1(l1, j1) * ID2(l2, j2);
              }
            }
          }
        }
      }
    }
  }
}

void KNDdeTorusCollocation::PhaseBOTH(KNVector& ph0, KNVector& ph1, KNVector& presol)
{
  ph0.clear();
  ph1.clear();
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          size_t idx1 = idxmap(j1, j2, i1, i2);
          // matrix multiplication
          for (size_t l2 = 0; l2 < ndeg2 + 1; l2++)
          {
            for (size_t l1 = 0; l1 < ndeg1 + 1; l1++)
            {
              const size_t idx2 = idxmap(l1, l2, i1, i2);
              for (size_t p = 0; p < NDIM; p++)
              {
                ph0(p + NDIM*idx1) += presol(p + NDIM * idx2) * ID1(l1, j1) * I2(l2, j2);
                ph1(p + NDIM*idx1) += presol(p + NDIM * idx2) * I1(l1, j1) * ID2(l2, j2);
              }
            }
          }
        }
      }
    }
  }
}

double KNDdeTorusCollocation::integrate(const KNVector& ph1, const KNVector& ph2)
{
  double res = 0.0;
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          size_t idx1 = idxmap(j1, j2, i1, i2);
          // matrix multiplication
          for (size_t l2 = 0; l2 < ndeg2 + 1; l2++)
          {
            for (size_t l1 = 0; l1 < ndeg1 + 1; l1++)
            {
              const size_t idx2 = idxmap(l1, l2, i1, i2);
              for (size_t p = 0; p < NDIM; p++)
              {
                res += ph1(p + NDIM * idx1) * I1(j1, l1) * I2(j2, l2) * ph2(p + NDIM * idx2);
              }
            }
          }
        }
      }
    }
  }
  return res;
}

void KNDdeTorusCollocation::star(KNVector& ph1, const KNVector& ph2)
{
  ph1.clear();
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          size_t idx1 = idxmap(j1, j2, i1, i2);
          // matrix multiplication
          for (size_t l2 = 0; l2 < ndeg2 + 1; l2++)
          {
            for (size_t l1 = 0; l1 < ndeg1 + 1; l1++)
            {
              const size_t idx2 = idxmap(l1, l2, i1, i2);
              for (size_t p = 0; p < NDIM; p++)
              {
                ph1(p + NDIM*idx2) += I1(l1, j1) * I2(l2, j2) * ph2(p + NDIM * idx1);
              }
            }
          }
        }
      }
    }
  }
}

double KNDdeTorusCollocation::IntegrateDIFF(KNVector& ph1, KNVector& ph2, KNVector& ph3)
{
  double res = 0.0;
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2 + 1; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1 + 1; j1++)
        {
          size_t idx1 = idxmap(j1, j2, i1, i2);
          // matrix multiplication
          for (size_t l2 = 0; l2 < ndeg2 + 1; l2++)
          {
            for (size_t l1 = 0; l1 < ndeg1 + 1; l1++)
            {
              const size_t idx2 = idxmap(l1, l2, i1, i2);
              for (size_t p = 0; p < NDIM; p++)
              {
                res += (ph1(p + NDIM * idx1) - ph2(p + NDIM * idx1)) * I1(j1, l1) * I2(j2, l2) * ph3(p + NDIM * idx2);
              }
            }
          }
        }
      }
    }
  }
  return res;
}

void KNDdeTorusCollocation::importSolution(KNVector& out, KNVector& in)
{
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1; j1++)
        {
          const size_t idx1 = idxmap(j1, j2, i1, i2);
          for (size_t p = 0; p < NDIM; p++)
          {
            out(p + NDIM*idx1) = in(p + NDIM * (j1 + i1 * ndeg1));
          }
        }
      }
    }
  }
}

void KNDdeTorusCollocation::importTangent(KNVector& out, KNVector& Re, KNVector& Im, double alpha)
{
  // alpha /= 2.0;
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t i1 = 0; i1 < nint1; i1++)
    {
      for (size_t j2 = 0; j2 < ndeg2; j2++)
      {
        for (size_t j1 = 0; j1 < ndeg1; j1++)
        {
          const size_t idx1 = idxmap(j1, j2, i1, i2);
          const double t1 = (mesh1(j1) + i1) / ((double)nint1);
          const double t2 = (mesh2(j2) + i2) / ((double)nint2);
          for (size_t p = 0; p < NDIM; p++)
          {
            out(p + NDIM*idx1) = cos(-alpha * t1 + 2 * M_PI * t2) * Re(p + NDIM * (j1 + i1 * ndeg1)) + sin(-alpha * t1 + 2 * M_PI * t2) * Im(p + NDIM * (j1 + i1 * ndeg1));
          }
        }
      }
    }
  }
  double norm = sqrt(integrate(out, out));
  out /= norm;
}

void KNDdeTorusCollocation::Save(const char* dat, const char* idx, const KNVector& in)
{
  // writing to file for plotting
  std::ofstream ff(dat);
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t j2 = 0; j2 < ndeg2; j2++)
    {
      for (size_t i1 = 0; i1 < nint1; i1++)
      {
        for (size_t j1 = 0; j1 < ndeg1; j1++)
        {
          const size_t idx1 = idxmap(j1, j2, i1, i2);
          for (size_t p = 0; p < NDIM; p++)
          {
            ff << in(p + NDIM*idx1) << "\t";
          }
        }
      }
      ff << "\n";
    }
  }

  std::ofstream gg(idx);
  for (size_t i1 = 0; i1 < nint1; i1++)
  {
    for (size_t j1 = 0; j1 < ndeg1; j1++)
    {
      const double t1 = (mesh1(j1) + i1) / ((double)nint1);
      gg << t1 << "\t";
    }
  }
  gg << "\n";
  for (size_t i2 = 0; i2 < nint2; i2++)
  {
    for (size_t j2 = 0; j2 < ndeg2; j2++)
    {
      const double t2 = (mesh2(j2) + i2) / ((double)nint2);
      gg << t2 << "\t";
    }
  }
  gg << "\n";
}

#undef NTAU
#undef NDIM

