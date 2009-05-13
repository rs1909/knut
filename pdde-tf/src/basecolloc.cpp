// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <cmath>
#include "basecolloc.h"
#include "polynomial.h"
#include "system.h"

#ifdef DEBUG
#include <iostream>
#include <fstream>
#endif /*DEBUG*/

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
      P_MESSAGE1("Unsupported degree of collocation polinomial is selected.");
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
      //      cout<<"in:poly_mul\\n";
      for (k = 0; k < t.Size(); k++)
      {
        if (k != i)
          poly_mul(poly, 1.0 / (t(i) - t(k)), -t(k) / (t(i) - t(k)));
        if (k != j)
          poly_mul(poly, 1.0 / (t(j) - t(k)), -t(k) / (t(j) - t(k)));
      }
      //      cout<<"out:poly_mul\\n";
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

#define NDIM ndim
#define NPAR npar
#define NINT nint
#define NDEG ndeg

PerSolColloc::PerSolColloc(System& _sys, const int _nint, const int _ndeg) :
    ndim(_sys.ndim()), npar(_sys.npar()),
    nint(_nint), ndeg(_ndeg),
    mesh(nint + 1), 
    time(nint*ndeg),
    timeMSH(ndeg*nint + 1),
    metric(ndeg + 1, ndeg + 1),
    metricPhase(ndeg + 1, ndeg + 1),
    col(ndeg),
    out(ndeg + 1),
    meshINT(ndeg + 1),
    lgr(ndeg+1, ndeg+1)
{
  sys = &_sys;
  for (int i = 0; i < nint + 1; i++) mesh(i) = i * 1.0 / nint;
  repr_mesh(meshINT);

  col_mesh(col);
  poly_int(metric, meshINT);   // works for any meshes
  poly_diff_int(metricPhase, meshINT);
  
  // computes the largrange coefficients
  for (int i = 0; i < ndeg+1; i++)
  {
    poly_coeff_lgr(lgr(i), meshINT, i);
  }
}

int PerSolColloc::meshlookup(const Vector& mesh, double t)
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

//----------------------------------------------------------------------------
//
// here are the mesh adaptation routines
//
//----------------------------------------------------------------------------

static void meshConstruct(Vector& newmesh, const Vector& oldmesh, const Vector& eqf)
{
//   for (int i = 1; i < eqf.Size()-1; i++) if (isnan(eqf(i))) std::cout<<i<<": nan ";
//   std::cout<<"first "<<eqf(1)<<" end "<<eqf(NINT)<<" ratio "<< eqf(1)/eqf(NINT)<<"\\n";
  // now computing the new mesh
  const int nint = oldmesh.Size()-1;
  newmesh(0) = 0.0;
  for (int i = 1; i < newmesh.Size()-1; i++)
  {
    const double t = eqf(nint)*i/(newmesh.Size()-1);
    const int idx = PerSolColloc::meshlookup( eqf, t );
    const double d = (t - eqf(idx))/(eqf(idx+1)-eqf(idx));
//     std::cout<<t<<":"<<d<<":"<<i<<":"<<idx<<":"<<mesh(idx) + d*(mesh(idx+1)-mesh(idx))<<" : "<<mesh(idx+1)-mesh(idx)<<"\\n";
    newmesh(i) = oldmesh(idx) + d*(oldmesh(idx+1)-oldmesh(idx));
//     if (eqf(i) < eqf(i-1)) std::cout<<"bad "<<eqf(i-1)<<", "<<eqf(i)<<"\\n";
//     if (newmesh(i) < newmesh(i-1)) std::cout<<"very bad "<<newmesh(i-1)<<", "<<newmesh(i)<<"\\n";
  }
  newmesh(newmesh.Size()-1) = 1.0;
}

static void meshAssess(Vector& eqf, const Vector& mesh, const Vector& profile, const Array1D< Array1D<double> >& lgr)
{
  int ndeg_ = lgr.Size() - 1;
  int nint_ = mesh.Size() - 1;
  P_ERROR_X1( profile.Size() % (ndeg_*nint_ + 1) == 0, "Wrong profile size.");
  int ndim_ = profile.Size() / (ndeg_*nint_ + 1);
  
  // compute the coeff of the highest degree term in each interval
  bool small_deri = true;
  const double hmach = 1e-6;
  Matrix hd(nint_+1, ndim_);
  for (int i = 0; i < nint_; i++)
  {
    for (int p = 0; p < ndim_; p++)
    {
      hd(i, p) = 0.0;
      for (int j = 0; j < ndeg_+1; j++)
      {
        hd(i,p) += lgr(j)(ndeg_)*profile( p + ndim_*( j + ndeg_*i ) );
      }
      // adjust by the mesh interval
      hd(i,p) /= pow(mesh(i+1)-mesh(i),ndeg_);
      if (fabs(hd(i,p)) > hmach) small_deri = false;
    }
  }
//   if (small_deri) std::cout<<"small derivatives\\n";
  // takes care of periodicity
  // this has to be changed when other boundary condition is used
  for (int p = 0; p < ndim_; p++)
  {
    hd(nint_,p) = hd(0,p);
  }
  // computes the (m+1)-th derivative.
  // The mesh modulo need not to be changed, when not periodic BC is used.
  for (int i = 0; i < nint_; i++)
  {
    double dtav;
    if ( i+2 < nint_ ) dtav = 0.5*(mesh(i+2)-mesh(i));
    else dtav = 0.5*(1.0+mesh((i+2)-nint_)-mesh(i));
    if( dtav < 0.0 ) std::cout<<"dtav<0\\n";
    for (int p = 0; p < ndim_; p++)
    {
      hd(i,p) = (hd(i+1,p) - hd(i,p))/dtav;
    }
  }
  // takes care of periodicity
  // this has to be changed when other boundary condition is used
  for (int p = 0; p < ndim_; p++)
  {
    hd(nint_,p) = hd(0,p);
  }
  // eqf contains the integral which has to be equidistributed
  P_ERROR_X1( eqf.Size() == nint_+1, "EQF has wrong size.");
  // when the derivatives are too small;
  if (small_deri) for (int i = 0; i < nint_+1; ++i) eqf(i) = i;
  else eqf(0) = 0.0;
  // computing eqf
  const double pwr=1.0/(ndeg_+1.0);
  for (int j=0; j < nint_; ++j)
  {
    double EP=0;
    if (j == 0)
    {
      for (int i = 0; i < ndim_; ++i)
      {
        EP+=pow(fabs(hd(nint_,i)),pwr);
      }
    } else
    {
      for (int i = 0; i < ndim_; ++i)
      {
        EP+=pow(fabs(hd(j-1,i)),pwr);
      }
    }
    double E=0;
    for (int i = 0; i < ndim_; ++i)
    {
      E+=pow(fabs(hd(j,i)),pwr);
    }
    eqf(j+1)=eqf(j)+0.5*(mesh(j+1)-mesh(j))*(E+EP);
  }
}

void PerSolColloc::meshAdapt_internal( Vector& newmesh, const Vector& profile )
{
  Vector eqf(NINT+1);
  meshAssess(eqf, mesh, profile, lgr);
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
      int idx = PerSolColloc::meshlookup( mesh, t );
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

void PerSolColloc::meshAdapt(Vector& newprofile, const Vector& profile, Vector& newtangent, const Vector& tangent)
{
  // saving the solution into solNu
#ifdef DEBUG
  Vector profile_tmp(profile);
#endif

  Vector newmesh(mesh);
  meshAdapt_internal(newmesh, profile);
  profileConvert(newprofile, newmesh, profile, mesh, lgr, NDIM);
  profileConvert(newtangent, newmesh, tangent, mesh, lgr, NDIM);
#ifdef DEBUG
  // printing the adapted profile
  std::ofstream file1("prof");
  file1 << std::scientific;
  file1.precision(12);
  std::ofstream file2("newprof");
  file2 << std::scientific;
  file2.precision(12);

  std::ofstream file3("gradient");
  file3 << std::scientific;
  file3.precision(12);
  for (int i=0; i<NINT; ++i)
  {
    for (int j=0; j<NDEG+1; ++j)
    {
      const double t1 = mesh(i) + j*(mesh(i+1)-mesh(i))/NDEG;
      const double t2 = newmesh(i) + j*(newmesh(i+1)-newmesh(i))/NDEG;
      file1<<t1<<"\\t";
      file2<<t2<<"\\t";
      file3<<t1<<"\\t"<<t2<<"\\t";
      for (int p=0; p<NDIM; ++p)
      {
	file1<<profile_tmp(p+NDIM*(j+NDEG*i))<<"\\t";
	file2<<profile(p+NDIM*(j+NDEG*i))<<"\\t";
      }
      file1<<"\\n";
      file2<<"\\n";
      file3<<"\\n";
    }
  }
#endif //DEBUG
  mesh = newmesh;
}

//----------------------------------------------------------------------------
//
// here are the integration routines
//
//----------------------------------------------------------------------------

void PerSolColloc::getMetric(Matrix& mt, const Vector& t)
{
  poly_int(mt, t);
}

void PerSolColloc::getDiffMetric(Matrix& mt, const Vector& t)
{
  poly_int(mt, t);
}

void PerSolColloc::star(Vector& out, const Vector& in, const Matrix& mt, const Vector& msh, int dim)
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

double PerSolColloc::integrate(const Vector& v1, const Vector& v2, const Matrix& mt, const Vector& msh, int dim)
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

void PerSolColloc::Star(Vector& V1, const Vector& V2)
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

  // Now we add (M udot)^* M\\phi + <udot, \\phi > ...
  for (int j = 0; j < NDIM; j++)
  {
    V1(NINT*NDEG*NDIM + j) += V2(NINT * NDEG * NDIM + j);
  }

#endif
}

double PerSolColloc::Integrate(const Vector& V1, const Vector& V2)
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

double PerSolColloc::IntegrateCont(const Vector& V1, const Vector& V2, const Vector& V3)
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

void PerSolColloc::PhaseStar(Vector& V1, const Vector& V2)
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

void PerSolColloc::PhaseRotStar(Vector& V1, const Vector& V2, const Array1D<int>& Re, const Array1D<int>& Im)
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

void PerSolColloc::pdMeshConvert(Vector& newprofile, Vector& newtangent, const Vector& oldprofile, const Vector& oldtangent)
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

void PerSolColloc::Import(Vector& newprofile, const Vector& oldprofile, const Vector& oldmesh, int old_ndeg)
{
  Vector old_meshINT(old_ndeg+1);
  repr_mesh(old_meshINT);

  Array1D< Array1D<double> > old_lgr(old_ndeg+1, old_ndeg+1);
  for (int i = 0; i < old_ndeg+1; i++)
  {
    poly_coeff_lgr(old_lgr(i), old_meshINT, i);
  }

  Vector eqf(oldmesh.Size());
  meshAssess(eqf, oldmesh, oldprofile, old_lgr);
  meshConstruct(mesh, oldmesh, eqf);
  profileConvert(newprofile, mesh, oldprofile, oldmesh, old_lgr, NDIM);
}

// it exports for CollocTR and PointTR, so no last value is necessary
void PerSolColloc::Export(Vector& outs, const Vector& mshint, const Vector& mshdeg, const Vector& in)
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
      // std::cout<<"int "<<i<<" "<<k<<"\\n";
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
