// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include "pderror.h"
#include "ncolloc.h"
#include "system.h"
#include "matrix.h"
#include "polynomial.h"

#include <cmath>


// #define MADD

// support routines

static void equidist(Vector& mesh)
{
  const int m = mesh.Size();
  for (int i = 0; i < m; i++) mesh(i) = (double)i / ((double)m - 1);
}

static void poly_gau(Vector& roots)
{
  const int m = roots.Size();

  Matrix a(m, m);
  /* construct the matrix */
  switch (m)
  {
    case 1:
      a(0, 0) = 0.0;
      break;
    case 2:
      a(0, 0) = -1.0 / 3;
      a(0, 1) = 0.0;
      break;
    case 3:
      a(0, 0) = 0.0;
      a(0, 1) = -3.0 / 5;
      a(0, 2) = 0.0;
      break;
    case 4:
      a(0, 0) = 3.0 / 35;
      a(0, 1) = 0.0;
      a(0, 2) = -6.0 / 7;
      a(0, 3) = 0.0;
      break;
    case 5:
      a(0, 0) = 0.0;
      a(0, 1) = 5.0 / 21;
      a(0, 2) = 0.0;
      a(0, 3) = -10.0 / 9;
      a(0, 4) = 0.0;
      break;
    case 6:
      a(0, 0) = -5.0 / 231;
      a(0, 1) = 0.0;
      a(0, 2) = 5.0 / 11;
      a(0, 3) = 0.0;
      a(0, 4) = -15.0 / 11;
      a(0, 5) = 0.0;
      break;
    case 7:
      a(0, 0) = 0.0;
      a(0, 1) = -35.0 / 429;
      a(0, 2) = 0.0;
      a(0, 3) = 105.0 / 143;
      a(0, 4) = 0.0;
      a(0, 5) = -21.0 / 13;
      a(0, 6) = 0;
      break;
    default:
      P_MESSAGE("Unsupported degree of collocation polinomial is selected.");
      return;
      break;
  }

  for (int i = 0; i < m; i++) a(m - 1, i) = -a(0, i);
  for (int i = 0; i < m - 1; i++)
  {
    for (int j = 0; j < m; j++)
    {
      if (i + 1 == j) a(i, j) = 1.0;
      else a(i, j) = 0.0;
    }
  }

  Vector wi(m);
  a.Eigval(roots, wi);

  /* scaling */
  for (int i = 0; i < m; i++) roots(i) = (roots(i) + 1.0) / 2.0;

  /* sorting */
  double tmp;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < (m - 1); j++)
    {
      if (roots(j) > roots(j + 1))
      {
        tmp = roots(j);
        roots(j) = roots(j + 1);
        roots(j + 1) = tmp;
      }
    }
  }
}

inline static void repr_mesh(Vector& V)
{
  equidist(V);
}

inline static void col_mesh(Vector& V)
{
  poly_gau(V);
}

static void poly_lgr(const Vector& t, Vector &out, double c)
{

  P_ASSERT_X(t.Size() == out.Size(), "poly_lgr: wrong dimensions");
  for (int i = 0; i < t.Size(); i++)
  {
    out(i) = 1.0;
    for (int j = 0; j < t.Size(); j++)
    {
      if (i != j)
      {
        out(i) *= (c - t(j)) / (t(i) - t(j));
      }
    }
  }
}

static void poly_dlg(const Vector& t, Vector& out, double c)
{
  int j, k, l;
  double f;

  P_ASSERT_X(t.Size() == out.Size(), "poly_dlg: wrong dimensions");

  for (j = 0; j < t.Size(); j++)
  {
    out(j) = 0.0;
    for (k = 0; k < t.Size(); k++)
    {
      if (k != j)
      {
        f = 1.0;
        for (l = 0; l < t.Size(); l++)
        {
          if ((l != k) && (l != j)) f *= (c - t(l)) / (t(j) - t(l));
        }
        out(j) += f / (t(j) - t(k));
      }
    }
  }
}

static void poly_d2lg(const Vector& t, Vector& out, double c)
{
  int i, j, k, l;
  double f;

  P_ASSERT_X(t.Size() == out.Size(), "poly_dlg: wrong dimensions");

  for (i = 0; i < t.Size(); i++)
  {
    out(i) = 0.0;
    for (l = 0; l < t.Size(); l++)
    {
      if (l==i) continue;
      for (k = 0; k < t.Size(); k++)
      {
        if (k==l) continue;
        if (k==i) continue;
        f = 1.0;
        for (j = 0; j < t.Size(); j++)
        {
          if (j==k) continue;
          if (j==l) continue;
          if (j==i) continue;
          f *= (c - t(j)) / (t(i) - t(j));
        }
        out(i) += f / (t(i) - t(k)) / (t(i) - t(l));
      }
    }
  }
}


static inline void poly_mul(Vector& pp, double bb, double aa)
{
  // pp * ( aa + bb*x )
  Vector tmp(pp.Size());
  for (int i = 0; i < pp.Size(); i++)
  {
    tmp(i) = pp(i) * bb;
    pp(i) *= aa;
  }
  for (int i = 1; i < pp.Size(); i++) pp(i) += tmp(i - 1);
}

static void poly_int(Matrix& out, const Vector& t)
{
  int i, j, k;
  Vector poly(2*t.Size());

  for (i = 0; i < t.Size(); i++)
  {
    for (j = 0; j < t.Size(); j++)
    {
      poly(0) = 1.0;
      for (k = 1; k < poly.Size(); k++) poly(k) = 0.0;
      //      poly.Print();
      // i,j az out matrix indexe
      //      cout<<"in:poly_mul\n";
      for (k = 0; k < t.Size(); k++)
      {
        if (k != i)
          poly_mul(poly, 1.0 / (t(i) - t(k)), -t(k) / (t(i) - t(k)));
        if (k != j)
          poly_mul(poly, 1.0 / (t(j) - t(k)), -t(k) / (t(j) - t(k)));
      }
      //      cout<<"out:poly_mul\n";
      //      t.Print();
      //      poly.Print();
      // integrate
      for (k = 0; k < poly.Size(); k++) poly(k) /= k + 1.0;
      out(i, j) = 0.0;
      // evaluate at x = 0..1
      for (k = 0; k < poly.Size(); k++) out(i, j) += poly(k);
    }
  }
}

static void poly_diff_int(Matrix& out, const Vector& t)
{
  Vector poly(2*t.Size());
  Vector poly_fin(2*t.Size());

  for (int i = 0; i < t.Size(); i++)
  {
    for (int j = 0; j < t.Size(); j++)
    {
      poly_fin.Clear();
      for (int s = 0; s < t.Size(); s++)
      {
        if (s != i)
        {
          poly(0) = 1.0;
          for (int r = 1; r < poly.Size(); r++) poly(r) = 0.0;

          for (int r = 0; r < t.Size(); r++)
          {
            if ((i != r) && (i != s)) poly_mul(poly, 1.0 / (t(i) - t(r)), -t(r) / (t(i) - t(r)));
            if (j != r)             poly_mul(poly, 1.0 / (t(j) - t(r)), -t(r) / (t(j) - t(r)));
          }
          // adding
          for (int r = 0; r < poly.Size(); r++) poly_fin(r) += poly(r) / (t(i) - t(s));
        }
      }
      // integrate
      for (int k = 0; k < poly_fin.Size(); k++) poly_fin(k) /= k + 1.0;
      out(i, j) = 0.0;
      // evaluate at x = 0..1
      for (int k = 0; k < poly_fin.Size(); k++) out(i, j) += poly_fin(k);
    }
  }
}

//
// the NColloc class
//

static inline int meshlookup(const Vector& mesh, double t)
{
  // binary search for in which interval is t-tau(k)
  int mid, low = 0, up = mesh.Size() - 1;
  while (up - low > 1)
  {
    mid = low + (up - low) / 2;
    if ((mesh(low) <= t) && (mesh(mid) > t)) up = mid;
    else low = mid;
  }
  return low;
}

#define NDIM ndim
#define NPAR npar
#define NTAU ntau
#define NINT nint
#define NDEG ndeg
#define NMAT nmat

NColloc::NColloc(System& _sys, const int _nint, const int _ndeg, int _nmat)
    :
    ndim(_sys.ndim()), npar(_sys.npar()), ntau(_sys.ntau()),
    nint(_nint), ndeg(_ndeg), nmat(_nmat),
    mesh(nint + 1), time(nint*ndeg),

    kk(ntau + 1, nint*ndeg), ee(ntau + 1, nint*ndeg),
    dd(ntau + 1, nint*ndeg), rr(ntau + 1, nint*ndeg),

    kkS(ntau + 1, nint*ndeg), eeS(ntau + 1, nint*ndeg),
    rrS(ntau + 1, nint*ndeg), ddS(ntau + 1, nint*ndeg),
    mmS(ntau + 1, nint*ndeg), szS(2, nint*ndeg),

    kkI(ntau + 1, nint*ndeg), eeI(ntau + 1, nint*ndeg),
    rrI(ntau + 1, nint*ndeg), ddI(ntau + 1, nint*ndeg),
    mmI(ntau + 1, nint*ndeg), szI(nmat + 1, nint*ndeg),

    tt(2*ntau + 1, ndeg + 1, ndeg*nint),
    timeMSH(ndeg*nint + 1),
    ttMSH(2*ntau, ndeg + 1, ndeg*nint + 1),
    kkMSH(ntau, nint*ndeg + 1),
    metric(ndeg + 1, ndeg + 1),
    metricPhase(ndeg + 1, ndeg + 1),
    col(ndeg),
    out(ndeg + 1),
    meshINT(ndeg + 1),
    lgr(ndeg+1, ndeg+1),
    p_tau(ntau, nint*ndeg), p_dtau(ntau, nint*ndeg),
    p_fx(ndim, nint*ndeg), p_dfx(ndim,ndim, nint*ndeg), p_dfp(ndim,1, nint*ndeg),
    p_fxRe(p_fx), p_fxIm(ndim, nint*ndeg),
    p_dfxRe(p_dfx), p_dfxIm(ndim,ndim, nint*ndeg),
    p_tauMSH(ntau, nint*ndeg+1), p_dtauMSH(ntau, nint*ndeg+1),
    p_fxMSH(ndim, nint*ndeg+1), p_dfxMSH(ndim,ndim, nint*ndeg+1), p_dfpMSH(ndim,1, nint*ndeg+1),
    p_dummy(0,0,nint*ndeg+1)
{
  sys = &_sys;
  for (int i = 0; i < NINT + 1; i++) mesh(i) = i * 1.0 / NINT;
  repr_mesh(meshINT);

  // computes the largrange coefficients
  for (int i = 0; i < ndeg+1; i++)
  {
    poly_coeff_lgr(lgr(i), meshINT, i);
  }

  col_mesh(col);
  poly_int(metric, meshINT);   // works for any meshes
  poly_diff_int(metricPhase, meshINT);
}

static void meshConstruct(Vector& newmesh, const Vector& oldmesh, const Vector& eqf)
{
//   for (int i = 1; i < eqf.Size()-1; i++) if (isnan(eqf(i))) std::cout<<i<<": nan ";
//   std::cout<<"first "<<eqf(1)<<" end "<<eqf(NINT)<<" ratio "<< eqf(1)/eqf(NINT)<<"\n";
  // now computing the new mesh
  const int nint = oldmesh.Size()-1;
  newmesh(0) = 0.0;
  for (int i = 1; i < newmesh.Size()-1; i++)
  {
    const double t = eqf(nint)*i/(newmesh.Size()-1);
    const int idx = meshlookup( eqf, t );
    const double d = (t - eqf(idx))/(eqf(idx+1)-eqf(idx));
//     std::cout<<t<<":"<<d<<":"<<i<<":"<<idx<<":"<<mesh(idx) + d*(mesh(idx+1)-mesh(idx))<<" : "<<mesh(idx+1)-mesh(idx)<<"\n";
    newmesh(i) = oldmesh(idx) + d*(oldmesh(idx+1)-oldmesh(idx));
//     if (eqf(i) < eqf(i-1)) std::cout<<"bad "<<eqf(i-1)<<", "<<eqf(i)<<"\n";
//     if (newmesh(i) < newmesh(i-1)) std::cout<<"very bad "<<newmesh(i-1)<<", "<<newmesh(i)<<"\n";
  }
  newmesh(newmesh.Size()-1) = 1.0;
}

void NColloc::meshAdapt( Vector& newmesh, const Vector& profile )
{
  // compute the coeff of the highest degree term in each interval
  bool small_deri = true;
  const double hmach = 1e-6;
  Matrix hd(NINT+1, NDIM);
  for (int i = 0; i < NINT; i++)
  {
    for (int p = 0; p < NDIM; p++)
    {
      hd(i, p) = 0.0;
      for (int j = 0; j < NDEG+1; j++)
      {
        hd(i,p) += lgr(j)(NDEG)*profile( p + NDIM*( j + NDEG*i ) );
      }
      // adjust by the mesh interval
      hd(i,p) /= pow(mesh(i+1)-mesh(i),NDEG);
      if (fabs(hd(i,p)) > hmach) small_deri = false;
    }
  }
//   if (small_deri) std::cout<<"small derivatives\n";
  // takes care of periodicity
  // this has to be changed when other boundary condition is used
  for (int p = 0; p < NDIM; p++)
  {
    hd(NINT,p) = hd(0,p);
  }
  // computes the (m+1)-th derivative.
  // The mesh modulo need not to be changed, when not periodic BC is used.
  for (int i = 0; i < NINT; i++)
  {
    double dtav;
    if ( i+2 < NINT ) dtav = 0.5*(mesh(i+2)-mesh(i));
    else dtav = 0.5*(1.0+mesh((i+2)-NINT)-mesh(i));
    if( dtav < 0.0 ) std::cout<<"dtav<0\n";
    for (int p = 0; p < NDIM; p++)
    {
      hd(i,p) = (hd(i+1,p) - hd(i,p))/dtav;
    }
  }
  // takes care of periodicity
  // this has to be changed when other boundary condition is used
  for (int p = 0; p < NDIM; p++)
  {
    hd(NINT,p) = hd(0,p);
  }
  // eqf contains the integral which has to be equidistributed
  Vector eqf(NINT+1);
  // when the derivatives are too small;
  if (small_deri) for (int i = 0; i < NINT+1; ++i) eqf(i) = i;
  else eqf(0) = 0.0;
  // computing eqf
  const double pwr=1.0/(NDEG+1.0);
  for (int j=0; j < NINT; ++j)
  {
    double EP=0;
    if (j == 0)
    {
      for (int i = 0; i < NDIM; ++i)
      {
        EP+=pow(fabs(hd(NINT,i)),pwr);
      }
    } else
    {
      for (int i = 0; i < NDIM; ++i)
      {
        EP+=pow(fabs(hd(j-1,i)),pwr);
      }
    }
    double E=0;
    for (int i = 0; i < NDIM; ++i)
    {
      E+=pow(fabs(hd(j,i)),pwr);
    }
    eqf(j+1)=eqf(j)+0.5*(mesh(j+1)-mesh(j))*(E+EP);
  }
  meshConstruct(newmesh, mesh, eqf);
}

static void profileConvert(Vector& newprofile, const Vector& newmesh, const Vector& profile, const Vector& mesh,
                           const Array1D< Array1D<double> >& old_lgr, const int ndim)
{
  const int old_nint = mesh.Size()-1;
  const int old_ndeg = (profile.Size()/ndim - 1)/old_nint;
  const int new_nint = newmesh.Size()-1;
  const int new_ndeg = (newprofile.Size()/ndim - 1)/new_nint;
  // creating the new profile with the same old meshINT
  for (int i = 0; i < new_nint; i++)
  {
    for (int j = 0; j < new_ndeg; j++)
    {
      const double t = newmesh(i) + j*(newmesh(i+1)-newmesh(i))/new_ndeg;
      int idx = meshlookup( mesh, t );
      const double d = (t - mesh(idx))/(mesh(idx+1)-mesh(idx));
      for (int p = 0; p < ndim; p++)
        newprofile(p+ndim*(j+i*new_ndeg)) = 0.0;
      for ( int k = 0; k < old_ndeg+1; k++)
      {
        double c = poly_eval(old_lgr(k), d);
        for (int p = 0; p < ndim; p++)
        {
          newprofile(p+ndim*(j+i*new_ndeg)) += c*profile(p+ndim*(k+idx*old_ndeg));
        }
      }
    }
  }
  for (int p = 0; p < ndim; p++)
    newprofile(p+ndim*(new_ndeg*(newmesh.Size()-1))) = profile(p+ndim*(old_ndeg*(mesh.Size()-1)));
}

void NColloc::profileAdapt(Vector& newprofile, const Vector& newmesh, const Vector& profile)
{
  profileConvert(newprofile, newmesh, profile, mesh, lgr, NDIM);
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
void NColloc::Init(const Vector& par, const Vector& /*sol*/)
{
  double *t = new double[NTAU+1];
  double *tMSH = new double[NTAU];

  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      int idx = j + i * NDEG;
      const double h0 = mesh(i + 1) - mesh(i);

      /* first, the derivative: x'(t)*/
      t[0] = mesh(i) + h0 * col(j); // xdot(c_i,j)
      time(idx) = t[0];
      timeMSH(idx) = mesh(i) + h0 * meshINT(j);
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
      sys->p_tau(p_tau, time, par);
      for (int k = 0; k < NTAU; k++)
      {
        p_tau(k,idx) /= par(0);
        P_ERROR_X2(p_tau(k,idx) <= NMAT, "\n NColloc::Init: DELAY > NMAT*PERIOD ", k);
        P_ERROR_X2(p_tau(k,idx) >= 0.0, "\n NColloc::Init: Negative DELAY ", k);
        t[1+k] = (t[0] - p_tau(k,idx)) - floor(t[0] - p_tau(k,idx));  // nem szetvalasztott

        // binary search for in which interval is t-p_tau(k,idx)
        const int low = meshlookup(mesh, t[1+k]);
        kk(k + 1, idx) = low;

        if (t[0] - p_tau(k,idx) >= 0) kkS(k + 1, idx) = low;
        else kkS(k + 1, idx) = low - NINT;

        kkI(k + 1, idx) = low + NINT * static_cast<int>(floor(t[0] - p_tau(k,idx)));

        const double hk = mesh(low + 1) - mesh(low);
        // x(t-\tau_i)
        poly_lgr(meshINT, out, (t[1+k] - mesh(low)) / hk);
        for (int l = 0; l < NDEG + 1; l++)
        {
          tt(1 + k, l, idx) = out(l);
        }
        // x'(t-\tau_i)
        poly_dlg(meshINT, out, (t[1+k] - mesh(low)) / hk);
        for (int l = 0; l < NDEG + 1; l++)
        {
          tt(NTAU + 1 + k, l, idx) = out(l) / hk;
        }
        // creating interpolation at the representation points
        tMSH[k] = mesh(i) + h0 * meshINT(j) - p_tau(k,idx) - floor(mesh(i) + h0 * meshINT(j) - p_tau(k,idx));
        const int lowMSH = meshlookup(mesh, tMSH[k]);
        kkMSH(k, idx) = lowMSH;
        const double hkMSH = mesh(lowMSH + 1) - mesh(lowMSH);
        poly_lgr(meshINT, out, (tMSH[k] - mesh(lowMSH)) / hkMSH);
        for (int l = 0; l < NDEG + 1; l++)
        {
          ttMSH(k, l, idx) = out(l);
        }
        // x'(t-\tau_i)
        poly_dlg(meshINT, out, (tMSH[k] - mesh(lowMSH)) / hkMSH);
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
  timeMSH(NDEG*NINT) = 1.0;
  for (int k = 0; k < NTAU; k++)
  {
    kkMSH(k, NDEG*NINT) = kkMSH(k, 0);
    for (int l = 0; l < NDEG + 1; l++)
    {
      ttMSH(k, l, NDEG*NINT) = ttMSH(k, l, 0);
      ttMSH(NTAU + k, l, NDEG*NINT) = ttMSH(NTAU + k, l, 0);
    }
  }
  delete[] tMSH;
  delete[] t;
}

void NColloc::Interpolate(Array3D<double>& solData, const Vector& sol)
{
  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;

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
}

void NColloc::InterpolateREAL(Array3D<double>& solData, const Vector& sol)
{
  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;

      for (int p = 0; p < NTAU; p++)
      {
        for (int k = 0; k < NDIM; k++)
        {
          // std::cout<<"InterR: "<<idx<<", "<<out(idx).Row()<<", "<<out(idx).Col()<<", "<<k<<", "<<p<<"\n";
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
    }
  }
}

void NColloc::InterpolateMSH(Array3D<double>& solData, const Vector& sol)
{
  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;

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
  for (int p = 0; p < NTAU; p++)
  {
    for (int k = 0; k < NDIM; k++)
    {
      // std::cout<<"InterR: "<<idx<<", "<<out(idx).Row()<<", "<<out(idx).Col()<<", "<<k<<", "<<p<<"\n";
      solData(k, p, NDEG*NINT)        = 0.0;
      solData(k, NTAU + p, NDEG*NINT) = 0.0;
      for (int l = 0; l < NDEG + 1; l++)
      {
        solData(k, p, NDEG*NINT)      += sol(k + NDIM * (l + kkMSH(p, NDEG * NINT) * NDEG)) * ttMSH(p, l, NDEG * NINT);
        solData(k, NTAU + p, NDEG*NINT) += sol(k + NDIM * (l + kkMSH(p, NDEG * NINT) * NDEG)) * ttMSH(NTAU + p, l, NDEG * NINT);
      }
    }
  }
}

// in complex form on 2*i th places are the reals and on 2*i+1 th places the imaginary parts

void NColloc::InterpolateCPLX(Array3D<double>& solDataRe, Array3D<double>& solDataIm, const Vector& sol)
{
  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;

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
}

void NColloc::RHS(Vector& rhs, const Vector& par, const Vector& sol, const Array3D<double>& solData)
{
  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    rhs(r) = sol(r + NDIM * NDEG * NINT) - sol(r);
  }

  sys->p_rhs(p_fx, time, solData, par);
  for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
  {
    for (int j = 0; j < NDEG; j++)
    {
      const int idx = j + i * NDEG;

      for (int k = 0; k < NDIM; k++)
        rhs(NDIM + k + NDIM*idx) = par(0) * p_fx(k, idx) - solData(k, NTAU, idx);
    }
  }
}

void NColloc::RHS_p(Vector& rhs, const Vector& par, const Vector& /*sol*/, const Array3D<double>& solData, int alpha)
{
  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    rhs(r) = 0.0;
  }

  if (alpha == 0)
  {
    sys->p_tau(p_tau, time, par);
    sys->p_dtau(p_dtau, time, par, alpha);
    sys->p_rhs(p_fx, time, solData, par);

    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;

        for (int p = 0; p < NDIM; p++)
        {
          rhs(NDIM + p + NDIM*idx) = -p_fx(p, idx);
        }
      }
    }

    for (int r = 0; r < NTAU; r++)
    {
      int nx=1, vx=r, np=0, vp;

      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;
          const double d = (p_dtau(r,idx) - p_tau(r,idx) / par(0));
          if (d != 0.0)
          {
            for (int p = 0; p < NDIM; p++)
            {
              for (int q = 0; q < NDIM; q++)
              {
                rhs(NDIM + p + NDIM*idx) += d * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
              }
            }
          }
        }
      }
    }
  }
  else
  {
    sys->p_tau(p_tau, time, par);
    sys->p_dtau(p_dtau, time, par, alpha);
    int nx, vx, np, vp;
    nx = 0, np = 1;
    vp = alpha;

    sys->p_deri(p_dfp, time, solData, par, nx, &vx, np, &vp, p_dummy);

    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;

        for (int k = 0; k < NDIM; k++) rhs(NDIM + k + NDIM*idx) = -par(0) * p_dfp(k, 0, idx);
      }
    }

    for (int r = 0; r < NTAU; r++)
    {
      nx = 1; np = 0; vx = r;
      sys->p_deri(p_dfx, time, solData, par, nx, &vx, np, &vp, p_dummy);
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;
          if (p_dtau(r, idx) != 0.0)
          {
            std::cout<<"P"<<alpha; // it is workng for the Glass and logistic eqns
            for (int p = 0; p < NDIM; p++)
            {
              for (int q = 0; q < NDIM; q++)
              {
                rhs(NDIM + p + NDIM*idx) += p_dtau(r, idx) * p_dfx(p, q, idx) * solData(q, NTAU + 1 + r, idx);
              }
            }
          }
        }
      }
    }
  }
}


void NColloc::RHS_x(SpMatrix& A, const Vector& par, const Vector& /*sol*/, const Array3D<double>& solData)
{

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
    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;

        for (int l = 0; l < NDEG + 1; l++) // degree
        {
          for (int p = 0; p < NDIM; p++)   // row
          {
            for (int q = 0; q < NDIM; q++)   //column
            {
              WRIDX(A, idx, i, j, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              if (k == 0)
              {
                if (p == q) WRDAT(A, idx, i, j, p, k, l, q) += tt(0, l, idx);
              }
              else
              {
                WRDAT(A, idx, i, j, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
            }
          }
        }
      }
    }
  }

}


//! its very different from all of them
void NColloc::StabJac(StabMatrix& AB, const Vector& par, const Array3D<double>& solData)
{

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
    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;

        for (int l = 0; l < NDEG + 1; l++)
        {
          for (int p = 0; p < NDIM; p++)
          {
            for (int q = 0; q < NDIM; q++)
            {
              if (mmI(k, idx) == 0)  // A matrix -- kkS(eeS(k,idx),idx) >= 0
              {
                WRIDXI(AB.getA0(), idx, i, j, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
                if (k == 0)
                {
                  if (p == q) WRDATI(AB.getA0(), idx, i, j, p, k, l, q) += tt(0, l, idx);
                }
                else
                {
                  WRDATI(AB.getA0(), idx, i, j, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                }
              }
              else // B matrices
              {
//                 std::cout<<"m"<<mmI(k, idx)<<" ";
                WRIDXI(AB.getAI(mmI(k, idx) - 1), idx, i, j, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
                if (k == 0)
                {
//          if( p == q ) WRDATSB(B,idx, i,j,p, k,l,q) -= tt(0,l,idx);
                }
                else
                {
                  WRDATI(AB.getAI(mmI(k, idx) - 1), idx, i, j, p, k, l, q) += par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                }
              }
            }
          }
        }
      }
    }
  }
}

// similar to RHS_x but with one boundary condition only and multiplied by Z
void NColloc::CharJac_x(SpMatrix& A, const Vector& par, const Array3D<double>& solData, double Z)
{

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
    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;

        const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;
        const double ZP = pow(Z, zpow);

        for (int l = 0; l < NDEG + 1; l++) // degree
        {
          for (int p = 0; p < NDIM; p++)   // row
          {
            for (int q = 0; q < NDIM; q++)   //column
            {
              WRIDX(A, idx, i, j, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              if (k == 0)             // A
              {
                if (p == q) WRDAT(A, idx, i, j, p, k, l, q) += tt(0, l, idx);
              }
              else
              {
                if (zpow > 0)   // -zB0 - z^2*B1 ... - z^n*BN
                {
                  WRDAT(A, idx, i, j, p, k, l, q) -= ZP * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                }
                else                         // A
                {
                  WRDAT(A, idx, i, j, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                }
              }
            }
          }
        }
      }
    }
  }
}


// this has to be changed only to packed complex.
void NColloc::CharJac_x(SpMatrix& A, const Vector& par, const Array3D<double>& solData, double Re, double Im)
{

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
    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;

        const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;
        for (int l = 0; l < NDEG + 1; l++) // degree
        {
          for (int p = 0; p < NDIM; p++)   // row
          {
            for (int q = 0; q < NDIM; q++)   //column
            {
              WRIDXCPLX(A, idx, i, j, p, 0, k, l, q, 0) = 2 * (q + NDIM * (l + NDEG * kk(k, idx)));
              WRIDXCPLX(A, idx, i, j, p, 0, k, l, q, 1) = 2 * (q + NDIM * (l + NDEG * kk(k, idx))) + 1;
              WRIDXCPLX(A, idx, i, j, p, 1, k, l, q, 0) = 2 * (q + NDIM * (l + NDEG * kk(k, idx)));
              WRIDXCPLX(A, idx, i, j, p, 1, k, l, q, 1) = 2 * (q + NDIM * (l + NDEG * kk(k, idx))) + 1;
              if (k == 0)             // A
              {
                if (p == q)
                {
                  WRDATCPLX(A, idx, i, j, p, 0, k, l, q, 0) += tt(0, l, idx);
                  WRDATCPLX(A, idx, i, j, p, 1, k, l, q, 1) += tt(0, l, idx);
                }
              }
              else
              { //kkS(ee(k,idx),idx) < 0
                if (zpow > 0)   // -zB0 - z^2*B1 ... - z^n*BN
                {
                  WRDATCPLX(A, idx, i, j, p, 0, k, l, q, 0) -= ZReP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                  WRDATCPLX(A, idx, i, j, p, 0, k, l, q, 1) += ZImP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                  WRDATCPLX(A, idx, i, j, p, 1, k, l, q, 0) -= ZImP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                  WRDATCPLX(A, idx, i, j, p, 1, k, l, q, 1) -= ZReP(zpow) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                }
                else                         // A
                {
                  WRDATCPLX(A, idx, i, j, p, 0, k, l, q, 0) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                  WRDATCPLX(A, idx, i, j, p, 1, k, l, q, 1) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
                }
              }
            }
          }
        }
      }
    }
  }
}

// REAL

void NColloc::CharJac_x_p(Vector& V, const Vector& par, const Array3D<double>& solData, const Array3D<double>& phiData, double Z, int alpha)
{

  V.Clear();

  // boundary conditions
  for (int r = 0; r < NDIM; r++)
  {
    V(r) = 0.0;
  }

  if (alpha == 0)   // a periodusido szerint deriv...
  {
    for (int k = 0; k < NTAU; k++)
    {
      sys->p_tau(p_tau, time, par);
      sys->p_dtau(p_dtau, time, par, alpha);
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
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;

          const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;
          const double ZP = pow(Z, zpow);

          for (int p = 0; p < NDIM; p++)
          {
            if (zpow > 0)
            {
              V(NDIM + p + NDIM*(j + NDEG*i)) += ZP * p_fx(p,idx);
            }
            else
            {
              V(NDIM + p + NDIM*(j + NDEG*i)) += p_fx(p,idx);
            }
          }
        }
      }
    }
  }
  else
  {
    for (int k = 0; k < NTAU; k++)
    {
      sys->p_tau(p_tau, time, par);
      sys->p_dtau(p_dtau, time, par, alpha);

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
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;
          const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;
          const double ZP = pow(Z, zpow);

          for (int p = 0; p < NDIM; p++)
          {
            if (zpow > 0)
            {
              V(NDIM + p + NDIM*(j + NDEG*i)) += ZP * p_fx(p,idx);
            }
            else
            {
              V(NDIM + p + NDIM*(j + NDEG*i)) += p_fx(p,idx);
            }
          }
        }
      }
    }
  }
}

// COMPLEX

void NColloc::CharJac_x_p(Vector& V, const Vector& par, const Array3D<double>& solData,
                          const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm,
                          double Re, double Im, int alpha)
{

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
    for (int k = 0; k < NTAU; k++)
    {
      sys->p_tau(p_tau, time, par);
      sys->p_dtau(p_dtau, time, par, alpha);
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
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;

          const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

          for (int p = 0; p < NDIM; p++)
          {
            if (zpow > 0)
            {
              V(2*(NDIM + p + NDIM*(j + NDEG*i)))     += ZReP(zpow) * p_fxRe(p,idx) - ZImP(zpow) * p_fxIm(p,idx);
              V(2*(NDIM + p + NDIM*(j + NDEG*i)) + 1) += ZImP(zpow) * p_fxRe(p,idx) + ZReP(zpow) * p_fxIm(p,idx);
            }
            else
            {
              V(2*(NDIM + p + NDIM*(j + NDEG*i)))     += p_fxRe(p,idx);
              V(2*(NDIM + p + NDIM*(j + NDEG*i)) + 1) += p_fxIm(p,idx);
            }
          }
        }
      }
    }
  }
  else
  {
    for (int k = 0; k < NTAU; k++)
    {
      sys->p_tau(p_tau, time, par);
      sys->p_dtau(p_dtau, time, par, alpha);

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
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;
          const int idxRe = 2 * (j + i * NDEG);
          const int idxIm = 2 * (j + i * NDEG) + 1;

          const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

          for (int p = 0; p < NDIM; p++)
          {
            if (zpow > 0)
            {
              V(2*(NDIM + p + NDIM*(j + NDEG*i)))     += ZReP(zpow) * p_fxRe(p,idx) - ZImP(zpow) * p_fxIm(p,idx);
              V(2*(NDIM + p + NDIM*(j + NDEG*i)) + 1) += ZImP(zpow) * p_fxRe(p,idx) + ZReP(zpow) * p_fxIm(p,idx);
            }
            else
            {
              V(2*(NDIM + p + NDIM*(j + NDEG*i)))     += p_fxRe(p,idx);
              V(2*(NDIM + p + NDIM*(j + NDEG*i)) + 1) += p_fxIm(p,idx);
            }
          }
        }
      }
    }
  }
}

void NColloc::CharJac_x_x(SpMatrix& A, const Vector& par, const Array3D<double>& solData, const Array3D<double>& phiData, double Z)
{
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
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        for (int ra = 0; ra < NDIM; ra++)
        {
          for (int rb = 0; rb < NDIM; rb++)
          {
            p_dfx(ra, rb, idx) += p_t_dfx(ra, rb, idx);
          }
        }
      }
    }

    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;
        const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;
        const double ZP = pow(Z, zpow);

        for (int l = 0; l < NDEG + 1; l++) // degree
        {
          for (int p = 0; p < NDIM; p++)   // row
          {
            for (int q = 0; q < NDIM; q++)   //column
            {
              WRIDX(A, idx, i, j, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              // This is valid only for k != 0
              if (zpow > 0)   // -zB
              {
                WRDAT(A, idx, i, j, p, k, l, q) -= ZP * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
              else                         // A
              {
                WRDAT(A, idx, i, j, p, k, l, q) -= par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
            }
          }
        }
      }
    }
  }
}

void NColloc::CharJac_x_x(SpMatrix& A, const Vector& par, const Array3D<double>& solData,
                          const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm,
                          double Re, double Im)
{
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
        for (int idx = 0; idx < NDEG*NINT; ++idx)
        {
          for (int ra = 0; ra < NDIM; ra++)
          {
            for (int rb = 0; rb < NDIM; rb++)
            {
              p_dfxRe(ra, rb, idx) += p_t_dfxRe(ra, rb, idx);
              p_dfxIm(ra, rb, idx) += p_t_dfxIm(ra, rb, idx);
            }
          }
        }
      }
      //std::cout<<"dfx "; dfx.Print();
    }
    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;
        const int idxRe = 2 * (j + i * NDEG);
        const int idxIm = 2 * (j + i * NDEG) + 1;
        const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;

        for (int l = 0; l < NDEG + 1; l++) // degree
        {
          for (int p = 0; p < NDIM; p++)   // row
          {
            for (int q = 0; q < NDIM; q++)   //column
            {
              WRIDXCPLXM(A, idx, i, j, p, 0, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              WRIDXCPLXM(A, idx, i, j, p, 1, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
              if (k == 0)             // A
              {
                // WRDATCPLX(A,idx, i,j,p,0, k,l,q,0) += tt(0,l,idx); // no derivative !!!!!!
                // WRDATCPLX(A,idx, i,j,p,1, k,l,q,1) += tt(0,l,idx);
              }
              else
              {
                if (zpow > 0)   // -zB
                {
                  WRDATCPLXM(A, idx, i, j, p, 0, k, l, q) -=
                    ZReP(zpow) * par(0) * p_dfxRe(p, q, idx) * tt(k, l, idx) - ZImP(zpow) * par(0) * p_dfxIm(p, q, idx) * tt(k, l, idx);
                  WRDATCPLXM(A, idx, i, j, p, 1, k, l, q) -=
                    ZImP(zpow) * par(0) * p_dfxRe(p, q, idx) * tt(k, l, idx) + ZReP(zpow) * par(0) * p_dfxIm(p, q, idx) * tt(k, l, idx);
                }
                else                         // A
                {
                  WRDATCPLXM(A, idx, i, j, p, 0, k, l, q) -= par(0) * p_dfxRe(p, q, idx) * tt(k, l, idx);
                  WRDATCPLXM(A, idx, i, j, p, 1, k, l, q) -= par(0) * p_dfxIm(p, q, idx) * tt(k, l, idx);
                }
              }
            }
          }
        }
      }
    }
  }
}

// this is for CharmatCPLX

void NColloc::CharJac_x_z(Vector& V, const Vector& par, const Array3D<double>& solData, const Vector& phi,
                          const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm, double Re, double Im)
{

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

    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;
        const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            if (zpow > 0)
            {
              V(2*(NDIM + p + NDIM*(j + NDEG*i)))     -= zpow * (ZReP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataRe(q, k, idx) -
                                                                 ZImP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataIm(q, k, idx));
              V(2*(NDIM + p + NDIM*(j + NDEG*i)) + 1) -= zpow * (ZImP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataRe(q, k, idx) +
                                                                 ZReP(zpow - 1) * par(0) * p_dfx(p, q, idx) * phiDataIm(q, k, idx));
            }
          }
        }
      }
    }
  }
}

//!
//! from now CharmatLPAUT
//!

void NColloc::CharJac_mB(SpMatrix& B, const Vector& par, const Array3D<double>& solData, double Z)
{

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
    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;

        const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;

        for (int l = 0; l < NDEG + 1; l++)
        {
          for (int p = 0; p < NDIM; p++)
          {
            for (int q = 0; q < NDIM; q++)
            {
              if (zpow > 0)  // B matrix -- kkS(eeS(k,idx),idx) < 0
              {
                WRIDXS(B, idx, i, j, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
                WRDATS(B, idx, i, j, p, k, l, q) -= zpow * ZP(zpow - 1) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
            }
          }
        }
      }
    }
  }
}

// same as CharJac_x_p, but only writes the B part

void NColloc::CharJac_mB_p(Vector& V, const Vector& par, const Array3D<double>& solData, const Array3D<double>& phiData, double Z, int alpha)
{
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
    for (int k = 0; k < NTAU; k++)
    {
      // START verbatim from CharJac_x_p
      sys->p_tau(p_tau, time, par);
      sys->p_dtau(p_dtau, time, par, alpha);
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
      // END verbatim from CharJac_x_p
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;
          const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

          for (int p = 0; p < NDIM; p++)
          {
            if (zpow > 0)
            {
              V(NDIM + p + NDIM*(j + NDEG*i)) += zpow * ZP(zpow - 1) * p_fx(p, idx);
            }
          }
        }
      }
    }
  }
  else
  {
    for (int k = 0; k < NTAU; k++)
    {
      // START verbatim from CharJac_x_p
      sys->p_tau(p_tau, time, par);
      sys->p_dtau(p_dtau, time, par, alpha);

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
      // END verbatim from CharJac_x_p
      for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
      {
        for (int j = 0; j < NDEG; j++)
        {
          const int idx = j + i * NDEG;
          const int zpow = (-kkI(k + 1, idx) + NINT - 1) / NINT;

          for (int p = 0; p < NDIM; p++)
          {
            if (zpow > 0)
            {
              V(NDIM + p + NDIM*(j + NDEG*i)) += zpow * ZP(zpow - 1) * p_fx(p, idx);
            }
          }
        }
      }
    }
  }
}

// like x_x but write bpart only
void NColloc::CharJac_mB_x(SpMatrix& B, const Vector& par, const Array3D<double>& solData, const Array3D<double>& phiData, double Z)
{
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
      for (int idx = 0; idx < NDEG*NINT; ++idx)
      {
        for (int ra = 0; ra < NDIM; ra++)
        {
          for (int rb = 0; rb < NDIM; rb++)
          {
            p_dfx(ra, rb, idx) += p_t_dfx(ra, rb, idx);
          }
        }
      }
    }
    // END verbatim from CharJac_x_x
    for (int i = 0; i < NINT; i++)   // i: interval; j: which collocation point
    {
      for (int j = 0; j < NDEG; j++)
      {
        const int idx = j + i * NDEG;
        const int zpow = (-kkI(k, idx) + NINT - 1) / NINT;

        for (int l = 0; l < NDEG + 1; l++) // degree
        {
          for (int p = 0; p < NDIM; p++)   // row
          {
            for (int q = 0; q < NDIM; q++)   //column
            {
              if (zpow > 0)  // A matrix -- kkS(eeS(k,idx),idx) >= 0
              {
                WRIDXS(B, idx, i, j, p, k, l, q) = q + NDIM * (l + NDEG * kk(k, idx));
                WRDATS(B, idx, i, j, p, k, l, q) -= zpow * ZP(zpow - 1) * par(0) * p_dfx(p, q, idx) * tt(k, l, idx);
              }
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

void NColloc::CharJac_MSHphi(Vector& V, const Vector& par, const Array3D<double>& solData)
{
  sys->p_rhs(p_fxMSH, timeMSH, solData, par);

  for (int idx = 0; idx < NDEG*NINT+1; ++idx)   // i: interval; j: which collocation point
  {
    for (int k = 0; k < NDIM; k++) V(k + NDIM*idx) = - p_fxMSH(k, idx);
  }
}

void NColloc::CharJac_MSHphi_p(Vector& V, const Vector& par, const Array3D<double>& solData, int alpha)
{
  V.Clear(); /// it is not cleared otherwise!!!!

  // boundary conditions
  if (alpha == 0)
  {
    for (int r = 0; r < NTAU; r++)
    {
      sys->p_tau(p_tauMSH, timeMSH, par);
      sys->p_dtau(p_dtauMSH, timeMSH, par, alpha);
      int nx, vx, np, vp;
      nx = 1;
      np = 0;
      vx = r;
      sys->p_deri(p_dfxMSH, timeMSH, solData, par, nx, &vx, np, &vp, p_dummy);

      for (int idx = 0; idx < NDEG*NINT+1; ++idx)   // i: interval; j: which collocation point
      {
        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
          {
            V(p + NDIM*idx) += ((p_dtauMSH(r,idx) - p_tauMSH(r,idx) / par(0)) / par(0)) * p_dfxMSH(p, q, idx) * solData(q, NTAU + r, idx);
          }
        }
      }
    }
  }
  else
  {
    for (int r = 0; r < NTAU; r++)
    {
      sys->p_tau(p_tauMSH, timeMSH, par);
      sys->p_dtau(p_dtauMSH, timeMSH, par, alpha);

      int nx, vx, np, vp;
      nx = 0, np = 1;
      vp = alpha;

      sys->p_deri(p_dfpMSH, timeMSH, solData, par, nx, &vx, np, &vp, p_dummy);

      nx = 1, np = 0;
      vx = r;
      sys->p_deri(p_dfxMSH, timeMSH, solData, par, nx, &vx, np, &vp, p_dummy);

      for (int idx = 0; idx < NDEG*NINT+1; ++idx)   // i: interval; j: which collocation point
      {
        for (int k = 0; k < NDIM; k++) V(k + NDIM*idx) = - p_dfpMSH(k,0,idx);

        for (int p = 0; p < NDIM; p++)
        {
          for (int q = 0; q < NDIM; q++)
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

//----------------------------------------------------------------------------
//
// here are the integration routines
//
//----------------------------------------------------------------------------

void NColloc::getMetric(Matrix& mt, const Vector& t)
{
  poly_int(mt, t);
}

void NColloc::getDiffMetric(Matrix& mt, const Vector& t)
{
  poly_int(mt, t);
}

void NColloc::star(Vector& out, const Vector& in, const Matrix& mt, const Vector& msh, int dim)
{
  const int t_deg = mt.Col() - 1;
  const int t_int = (out.Size() / dim - 1) / t_deg;
  P_ERROR(in.Size() == out.Size());
  P_ERROR(out.Size() == dim*(t_deg*t_int + 1));
  P_ERROR(msh.Size() == t_int + 1);
  out.Clear();
  for (int i = 0; i < t_int; ++i)
  {
    const double dx = msh(i + 1) - msh(i);
    // ez itt a matrixszorzas
    for (int k = 0; k < t_deg + 1; ++k)
    {
      for (int l = 0; l < t_deg + 1; ++l)
      {
        for (int j = 0; j < dim; ++j)
        {
          out(j + dim*(k + i*t_deg)) += dx * mt(k, l) * in(j + dim * (l + i * t_deg));
        }
      }
    }
  }
#ifdef MADD // whether we need to add the headpoint ?
  for (int j = 0; j < dim; j++)
  {
    out(t_int*t_deg*dim + j) += in(t_int * t_deg * dim + j);
  }
#endif
}

double NColloc::integrate(const Vector& v1, const Vector& v2, const Matrix& mt, const Vector& msh, int dim)
{
  double res = 0.0;
  const int t_deg = mt.Col() - 1;
  const int t_int = (v1.Size() / dim - 1) / t_deg;
  P_ERROR(v1.Size() == v2.Size());
  P_ERROR(v1.Size() == dim*(t_deg*t_int + 1));
  P_ERROR(msh.Size() == t_int + 1);
  for (int i = 0; i < t_int; ++i)
  {
    const double dx = msh(i + 1) - msh(i);
    // ez itt a matrixszorzas
    for (int k = 0; k < t_deg; ++k)
    {
      for (int l = 0; l < t_deg; ++l)
      {
        for (int j = 0; j < dim; ++j)
        {
          res += v1(j + dim * (k + i * t_deg)) * dx * mt(k, l) * v2(j + dim * (l + i * t_deg));
        }
      }
    }
  }
#ifdef MADD // whether we need to add the headpoint ?
  for (int j = 0; j < dim; j++)
  {
    res += v1(t_int * t_deg * dim + j) * v2(t_int * t_deg * dim + j);
  }
#endif
  return res;
}

void NColloc::Star(Vector& V1, const Vector& V2)
{
  V1.Clear();
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          V1(j + NDIM*(k + i*NDEG)) += dx * metric(k, l) * V2(j + NDIM * (l + i * NDEG));
        }
      }
    }
  }

#ifdef MADD // whether we need to add the headpoint ?

  // Now we add (M udot)^* M\phi + <udot, \phi > ...
  for (int j = 0; j < NDIM; j++)
  {
    V1(NINT*NDEG*NDIM + j) += V2(NINT * NDEG * NDIM + j);
  }

#endif
}

double NColloc::Integrate(const Vector& V1, const Vector& V2)
{
  double res = 0.0, head = 0.0;
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          res += dx * V1(j + NDIM * (k + i * NDEG)) * metric(k, l) * V2(j + NDIM * (l + i * NDEG));
        }
      }
    }
  }

#ifdef MADD
  for (int j = 0; j < NDIM; j++)
  {
    head += V1(NINT * NDEG * NDIM + j) * V2(NINT * NDEG * NDIM + j);
  }
#endif

  return res + head;
}

double NColloc::IntegrateCont(const Vector& V1, const Vector& V2, const Vector& V3)
{
  double res = 0.0, head = 0.0;
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          res += dx * V1(j + NDIM * (k + i * NDEG)) * metric(k, l) * (V2(j + NDIM * (l + i * NDEG)) - V3(j + NDIM * (l + i * NDEG)));
        }
      }
    }
  }
#ifdef MADD
  for (int j = 0; j < NDIM; j++)
  {
    head += V1(NINT * NDEG * NDIM + j) * (V2(NINT * NDEG * NDIM + j) - V3(NINT * NDEG * NDIM + j));
  }
#endif
  return res + head;
}

void NColloc::PhaseStar(Vector& V1, const Vector& V2)
{
  V1.Clear();
  for (int i = 0; i < NINT; i++)
  {
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          V1(j + NDIM*(k + i*NDEG)) += metricPhase(k, l) * V2(j + NDIM * (l + i * NDEG));
        }
      }
    }
  }
}

void NColloc::PhaseRotStar(Vector& V1, const Vector& V2, const Array1D<int>& Re, const Array1D<int>& Im)
{
  V1.Clear();
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < Re.Size(); j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          V1(Re(j) + NDIM*(k + i*NDEG)) -= dx * metric(k, l) * V2(Im(j) + NDIM * (l + i * NDEG));
          V1(Im(j) + NDIM*(k + i*NDEG)) += dx * metric(k, l) * V2(Re(j) + NDIM * (l + i * NDEG));
        }
      }
    }
  }
}

void NColloc::pdMeshConvert(Vector& newprofile, Vector& newtangent, const Vector& oldprofile, const Vector& oldtangent)
{
  Vector tmp_mesh(2*NINT+1);
  Vector tmp_profile(NDIM*(2*NINT*NDEG+1));
  Vector tmp_tangent(NDIM*(2*NINT*NDEG+1));
  for (int i = 0; i < NINT; ++i)
  {
    tmp_mesh(i) = 0.0 + 0.5*mesh(i);
    tmp_mesh(NINT+i) = 0.5 + 0.5*mesh(i);
    for (int j = 0; j < NDEG; ++j)
    {
      for (int p = 0; p < NDIM; ++p)
      {
        tmp_profile(p+NDIM*(j+i*NDEG))        = oldprofile(p+NDIM*(j+i*NDEG));
        tmp_profile(p+NDIM*(j+(i+NINT)*NDEG)) = oldprofile(p+NDIM*(j+i*NDEG));
        tmp_tangent(p+NDIM*(j+i*NDEG))        = oldtangent(p+NDIM*(j+i*NDEG));
        tmp_tangent(p+NDIM*(j+(i+NINT)*NDEG)) = -oldtangent(p+NDIM*(j+i*NDEG));
      }
    }
  }
  tmp_mesh(2*NINT) = 1.0;
  for (int p = 0; p < NDIM; ++p)
  {
    tmp_profile(p+NDIM*(2*NINT*NDEG)) = tmp_profile(p);
    tmp_tangent(p+NDIM*(2*NINT*NDEG)) = tmp_tangent(p);
  }
  // constructing the new mesh
  Vector eqf(tmp_mesh.Size());
  for (int i = 0; i < eqf.Size(); ++i) eqf(i) = i;
  meshConstruct(mesh, tmp_mesh, eqf);
  profileConvert(newprofile, mesh, tmp_profile, tmp_mesh, lgr, NDIM);
  profileConvert(newtangent, mesh, tmp_tangent, tmp_mesh, lgr, NDIM);
}

void NColloc::Import(Vector& newprofile, const Vector& oldprofile, const Vector& oldmesh, int old_ndeg)
{
  Vector old_meshINT(old_ndeg+1);
  repr_mesh(old_meshINT);

  Array1D< Array1D<double> > old_lgr(old_ndeg+1, old_ndeg+1);
  for (int i = 0; i < old_ndeg+1; i++)
  {
    poly_coeff_lgr(old_lgr(i), old_meshINT, i);
  }

  Vector eqf(oldmesh.Size());
  for (int i = 0; i < eqf.Size(); ++i) eqf(i) = i;
  meshConstruct(mesh, oldmesh, eqf);
  profileConvert(newprofile, mesh, oldprofile, oldmesh, old_lgr, NDIM);
}

// it exports for CollocTR and PointTR, so no last value is necessary
void   NColloc::Export(Vector& outs, const Vector& mshint, const Vector& mshdeg, const Vector& in)
{
  int nint_ = mshint.Size() - 1;
  int ndeg_ = mshdeg.Size() - 1;
  Vector in_mesh(NDEG + 1);
  Vector in_lgr(NDEG + 1);

  for (int i = 0; i < NDEG + 1; i++) in_mesh(i) = i * 1.0 / NDEG;

  for (int i = 0; i < nint_; i++)
  {
    for (int j = 0; j < ndeg_; j++)
    {
      double t = mshint(i) + mshdeg(j) / nint_;
      int k = meshlookup(mesh, t);
      // std::cout<<"int "<<i<<" "<<k<<"\n";
      double c = (t - mesh(k)) / (mesh(k + 1) - mesh(k));  // mesh is the interval mesh in the class

      poly_lgr(in_mesh, in_lgr, c);
      // in_lgr.Print();
      for (int p = 0; p < NDIM; p++)
      {
        outs(p + NDIM*(j + i*ndeg_)) = 0.0;
        for (int r = 0; r < NDEG + 1; r++)
        {
          outs(p + NDIM*(j + i*ndeg_)) += in(p + NDIM * (r + k * NDEG)) * in_lgr(r);
        }
      }
    }
  }
}
