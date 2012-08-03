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
#include "config.h"
#ifdef GINAC_FOUND
#include "vf.h"
#endif

#ifndef _WIN32
extern "C"
{
#include <dlfcn.h>
}
#else
#include <windows.h>
#endif

extern "C" {
  typedef size_t(*tp_sys_ndim)();
  typedef size_t(*tp_sys_npar)();
  typedef size_t(*tp_sys_ntau)();
  typedef size_t(*tp_sys_nderi)();
  typedef size_t(*tp_sys_nevent)();
  typedef void(*tp_sys_tau)(KNVector& out, double t, const KNVector& par);
  typedef void(*tp_sys_dtau)(KNVector& out, double t, const KNVector& par, size_t vp);
  typedef void(*tp_sys_mass)(KNArray1D<double>& out);
  typedef void(*tp_sys_rhs)(KNVector& out, double t, const KNMatrix& x, const KNVector& par);
  typedef void(*tp_sys_deri)(KNMatrix& out, double t, const KNMatrix& x, const KNVector& par, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNMatrix& v);
  typedef void(*tp_sys_p_tau)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par );
  typedef void(*tp_sys_p_dtau)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray1D<double>& par, size_t vp );
  typedef void(*tp_sys_p_rhs)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel );
  typedef void(*tp_sys_p_deri)( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par, size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv );
  typedef void(*tp_sys_p_event)( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNArray1D<double>& par );
  typedef void(*tp_sys_stpar)(KNVector& par);
  typedef void(*tp_sys_stsol)(KNVector& out, double t);
  typedef void(*tp_sys_p_stsol)(KNArray2D<double>& out, const KNArray1D<double>& time);
  typedef void(*tp_sys_parnames)(const char *names[]);
}

class KNSystem
{
  public:

    KNSystem(const std::string& sysName, const std::string& sysType, size_t usederi);
    ~KNSystem();
    static void compileSystem(const std::string& cxxfile, const std::string& shobj, const std::string& executableDir);
    static void generateSystem(const std::string& vffile, const std::string& shobj, const std::string& executableDir);
    static bool makeSystem(std::string& soname, const std::string& sysName, const std::string& sysType, const std::string& executableDir);

    size_t  ndim() const ;
    size_t  npar() const ;
    size_t  ntau() const ;
    // Vectorized versions
    void   p_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par );
    void   p_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp );
    void   mass(KNArray1D<double>& out);
    void   p_rhs( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel );
    void   p_deri( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par,
                   size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv );
    // Setting the starting point
    void   stpar(KNVector& par) const;
    void   stsol(KNArray2D<double>& out, const KNArray1D<double>& time) const;
    void   parnames(const char *names[]) const;

  private:
    void   p_discrderi( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& p_xx, const KNVector& par, size_t sel,
                        size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& p_vv );
    void   p_resize(size_t sz)
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

    size_t     nderi;

    KNVector  f, f_eps;
    KNVector  f2, f_eps2;
    KNMatrix  xx_eps;
    KNMatrix  xx_eps2;
    KNVector  par_eps;
    KNMatrix  dxx2, dxx_eps2;
    KNVector  vt;

    size_t            p_size;
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
    tp_sys_mass     v_mass;
    tp_sys_rhs      v_rhs;
    tp_sys_deri     v_deri;
    tp_sys_p_tau    v_p_tau;
    tp_sys_p_dtau   v_p_dtau;
    tp_sys_p_rhs    v_p_rhs;
    tp_sys_p_deri   v_p_deri;
    tp_sys_p_event  v_p_event;
    tp_sys_stpar    v_stpar;
    tp_sys_stsol    v_stsol;
    tp_sys_p_stsol  v_p_stsol;
    tp_sys_parnames v_parnames;

#ifdef GINAC_FOUND
    void makeSymbolic(const std::string& vffile);
    bool useVectorField;  // use the vector field file directly, no compilation
    size_t ex_ndim;
    size_t ex_npar;
    size_t ex_ntau;
    size_t ex_nevent;
    std::vector<GiNaC::ex> ex_tau;
    std::vector<GiNaC::ex> ex_tau_p;
    std::vector<GiNaC::ex> ex_rhs;
    std::vector<double> ex_mass;
    std::vector<GiNaC::ex> ex_rhs_p;
    std::vector<GiNaC::ex> ex_rhs_x;
    std::vector<GiNaC::ex> ex_rhs_xp;
    std::vector<GiNaC::ex> ex_rhs_xx;
    std::vector<double> ex_stpar;
    std::vector<GiNaC::ex> ex_stsol;
    std::vector<std::string> ex_parnames;
    
    size_t sym_ndim();
    size_t sym_npar();
    size_t sym_ntau();
    void sym_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par );
    void sym_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp );
    void sym_rhs( KNArray2D<double>& out_p, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel );
    void sym_mass(KNArray1D<double>& out) const;
    void sym_deri( KNArray3D<double>& out_p, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par,
                size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv );
    void sym_stpar(KNVector& par) const;
    void sym_stsol(KNArray2D<double>& out, const KNArray1D<double>& time) const;
    void sym_parnames(const char *names[]) const;
#endif
    static bool workingCompiler;
};

inline KNSystem::tdlhand KNSystem::tdlopen(const char* fname)
{
#ifndef _WIN32
  return dlopen(fname, RTLD_NOW);
#else
  return LoadLibrary(fname);
#endif
}

inline void* KNSystem::tdlsym(tdlhand h, const char* sname)
{
#ifndef _WIN32
  void * res = dlsym(h, sname);
  return res;
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
