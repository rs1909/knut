// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef SYSTEM_H
#define SYSTEM_H

#include "matrix.h"

#ifndef _WIN32
extern "C"
{
#include <dlfcn.h>
}
#else
#include <windows.h>
#endif

class KNSystem
{
  public:

    KNSystem(const std::string& shobj, int usederi);
    ~KNSystem();
    static void compileSystem(const std::string& cxxfile, const std::string& shobj, const std::string& executableDir);
    static void generateSystem(const std::string& vffile, const std::string& executableDir);
    static void makeSystem(const std::string& shobj, const std::string& executableDir);

    int    ndim() const 
    {
      return (*v_ndim)();
    }
    int    npar() const 
    {
      return (*v_npar)();
    }
    int    ntau() const
    {
      return (*v_ntau)();
    }
    // Vectorized versions
    void   p_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par )
    {
      if (v_p_tau != 0) (*v_p_tau)(out, time, par);
      else
      {
        for (int i=0; i<time.size(); ++i)
        {
          KNVector tout(out, i);
          (*v_tau)(tout, time(i), par);
        }
      }
    }
    void   p_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, int vp )
    {
      if (v_p_dtau != 0) (*v_p_dtau)(out, time, par, vp);
      else
      {
        for (int i=0; i<time.size(); ++i)
        {
          KNVector tout(out, i);
          (*v_dtau)(tout, time(i), par, vp);
        }
      }
    }
    void   p_rhs( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, int sel )
    {
      if (v_p_rhs != 0) (*v_p_rhs)(out, time, x, par, sel);
      else
      {
        for (int i=0; i<time.size(); ++i)
        {
          KNVector vout(out, i);
          KNMatrix xxin(x, i);
          (*v_rhs)(vout, time(i), xxin, par);
        }
      }
    }
    void   p_deri( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par,
                   int sel, int nx, const int* vx, int np, const int* vp, const KNArray3D<double>& vv )
    {
      if ((nderi == 2) || ((nderi == 1) && ((nx == 1 && np == 0) || (nx == 0 && np == 1))))
      {
        if (v_p_deri != 0)
        {
          (*v_p_deri)(out, time, x, par, sel, nx, vx, np, vp, vv);
        } else
        {
          for (int i=0; i<time.size(); ++i)
          {
            KNMatrix mout(out, i);
            KNMatrix xxin(x, i);
            KNMatrix vvin(vv, i);
            (*v_deri)(mout, time(i), xxin, par, nx, vx, np, vp, vvin);
          }
        }
      } else
      {
        p_discrderi(out, time, x, par, sel, nx, vx, np, vp, vv);
      }
    }
    // Setting the starting point
    void   stpar(KNVector& par) const
    {
      (*v_stpar)(par);
    }
    void   stsol(KNVector& out, double t) const
    {
      (*v_stsol)(out, t);
    }
    void   parnames(const char *names[]) const
    {
      if (v_parnames != 0) (*v_parnames)(names);
    }

  private:
    void   p_discrderi( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& p_xx, const KNVector& par, int sel,
                        int nx, const int* vx, int np, const int* vp, const KNArray3D<double>& p_vv );

    void   p_resize(int sz)
    {
      if ( sz > p_size )
      {
        p_size = sz;
        p_fx.init(ndim(), p_size);
        p_fx_eps.init(ndim(), p_size);
        p_dfx.init(ndim(), ndim(), p_size);
        p_dfx_eps.init(ndim(), ndim(), p_size);
        p_xx_eps.init(ndim(), 2*ntau()+1, p_size);
        p2_dfx.init(ndim(), ndim(), p_size);
        p2_dfx_eps.init(ndim(), ndim(), p_size);
        p2_xx_eps.init(ndim(), 2*ntau()+1, p_size);
      }
    }

    typedef int(*tp_sys_ndim)();
    typedef int(*tp_sys_npar)();
    typedef int(*tp_sys_ntau)();
    typedef int(*tp_sys_nderi)();
    typedef int(*tp_sys_nevent)();
    typedef void(*tp_sys_tau)(KNVector& out, double t, const KNVector& par);
    typedef void(*tp_sys_dtau)(KNVector& out, double t, const KNVector& par, int vp);
    typedef void(*tp_sys_rhs)(KNVector& out, double t, const KNMatrix& x, const KNVector& par);
    typedef void(*tp_sys_deri)(KNMatrix& out, double t, const KNMatrix& x, const KNVector& par, int nx, const int* vx, int np, const int* vp, const KNMatrix& v);
    typedef void(*tp_sys_p_tau)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par );
    typedef void(*tp_sys_p_dtau)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, int vp );
    typedef void(*tp_sys_p_rhs)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, int sel );
    typedef void(*tp_sys_p_deri)( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, int sel, int nx, const int* vx, int np, const int* vp, const KNArray3D<double>& vv );
    typedef void(*tp_sys_p_event)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par );
    typedef void(*tp_sys_stpar)(KNVector& par);
    typedef void(*tp_sys_stsol)(KNVector& out, double t);
    typedef void(*tp_sys_parnames)(const char *names[]);

    typedef void(*FPTR)();
    union punned
    {
      void *obj;
      FPTR fun;
    };
    FPTR fptr(void * ptr);

#ifndef _WIN32
    typedef void*   tdlhand;
#else
    typedef HMODULE tdlhand;
#endif

    tdlhand     tdlopen(const char* fname);
    void*       tdlsym(tdlhand h, const char* sname);
    int         tdlclose(tdlhand h);
    const char* tdlerror();

    tdlhand handle;
    const char* error;

    int     nderi;

    KNVector  f, f_eps;
    KNVector  f2, f_eps2;
    KNMatrix  xx_eps;
    KNMatrix  xx_eps2;
    KNVector  par_eps;
    KNMatrix  dxx2, dxx_eps2;
    KNVector  vt;

    int             p_size;
    KNArray2D<double> p_fx,  p_fx_eps;
    KNArray3D<double> p_dfx, p_dfx_eps;
    KNArray3D<double> p_xx_eps;
    KNArray3D<double> p2_dfx, p2_dfx_eps;
    KNArray3D<double> p2_xx_eps;


    tp_sys_ndim     v_ndim;
    tp_sys_npar     v_npar;
    tp_sys_ntau     v_ntau;
    tp_sys_nderi    v_nderi;
    tp_sys_nevent   v_nevent;
    tp_sys_tau      v_tau;
    tp_sys_dtau     v_dtau;
    tp_sys_rhs      v_rhs;
    tp_sys_deri     v_deri;
    tp_sys_p_tau    v_p_tau;
    tp_sys_p_dtau   v_p_dtau;
    tp_sys_p_rhs    v_p_rhs;
    tp_sys_p_deri   v_p_deri;
    tp_sys_p_event  v_p_event;
    tp_sys_stpar    v_stpar;
    tp_sys_stsol    v_stsol;
    tp_sys_parnames v_parnames;

};

inline KNSystem::FPTR KNSystem::fptr(void * ptr)
{
  punned tmp;
  tmp.obj = ptr;
  return tmp.fun;
}

inline KNSystem::tdlhand KNSystem::tdlopen(const char* fname)
{
#ifndef _WIN32
  return dlopen(fname, RTLD_LAZY);
#else
  return LoadLibrary(fname);
#endif
}

inline void* KNSystem::tdlsym(tdlhand h, const char* sname)
{
#ifndef _WIN32
  return dlsym(h, sname);
#else
  return (void*) GetProcAddress(h, sname);
#endif
}

inline int KNSystem::tdlclose(KNSystem::tdlhand h)
{
#ifndef _WIN32
  return dlclose(h);
#else
  return (int) FreeLibrary(h);
#endif
}

inline const char* KNSystem::tdlerror()
{
#ifndef _WIN32
  return dlerror();
#else
  DWORD errcode = GetLastError();
  SetLastError(0);
  if (errcode != 0) return "Windows KNSystem Error\n";
  else return 0;
#endif
}

#endif
