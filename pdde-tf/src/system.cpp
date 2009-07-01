// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sys/stat.h>
#include <sstream>

#include "system.h"
#include "matrix.h"
#include "pointtype.h"
#include "config.h"

#ifdef __APPLE__
#include <CoreFoundation/CFBase.h>
#include <CoreFoundation/CFString.h>
#include <CoreFoundation/CFBundle.h>
#endif

System::System(const std::string& shobj)
{
  std::string objname(shobj);
  if (shobj.find('/') == std::string::npos) objname = "./" + shobj;
#ifdef WIN32
  std::string::size_type index = 0;
  while ((index = objname.find('/', index)) != std::string::npos) objname.at(index) = '\\';
#endif

  // Finding the executable directory
  std::string executableFile;
#ifdef __APPLE__
  CFURLRef bundleURL = CFBundleCopyExecutableURL(CFBundleGetMainBundle());
  if(bundleURL)
  {
    CFStringRef cfPath = CFURLCopyFileSystemPath(bundleURL, kCFURLPOSIXPathStyle);
    size_t ssize = static_cast<size_t>(CFStringGetLength(cfPath)+1);
    char *buf = new char[ssize];
    if(cfPath) CFStringGetCString(cfPath, buf, 512, kCFStringEncodingASCII);
    executableFile = buf;
    delete[] buf;
  }
#elif __linux__
  std::string workingDir;
  std::ostringstream procf;
  procf << "/proc/" << getpid() << "/exe";
  char *buf = new char[512];
  const ssize_t bsfn = readlink(procf.str().c_str(), buf, 511);
  if ( bsfn != -1) { buf[bsfn] = '\0'; executableFile = buf; }
  delete[] buf;
#endif

  std::string executableDir(executableFile);
  std::string::size_type sidx = executableDir.find_last_of('/');
  if (sidx != std::string::npos) executableDir.erase(sidx,std::string::npos);
  
  makeSystem(shobj, executableDir);
  handle = tdlopen(objname.c_str());
  P_ERROR_X5(handle != 0, "Cannot open system definition file. Error code", tdlerror(), ". The offending file was '", objname, "'.");

  tdlerror();    /* Clear any existing error */
  v_ndim = (tp_sys_ndim) fptr(tdlsym(handle, "sys_ndim"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_ndim(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_npar = (tp_sys_npar) fptr(tdlsym(handle, "sys_npar"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_npar(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_ntau = (tp_sys_ntau) fptr(tdlsym(handle, "sys_ntau"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_ntau(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_nderi = (tp_sys_nderi) fptr(tdlsym(handle, "sys_nderi"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_nderi(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_tau = (tp_sys_tau) fptr(tdlsym(handle, "sys_tau"));
  found_tau = ((error = tdlerror()) == 0);

  tdlerror();    /* Clear any existing error */
  v_dtau = (tp_sys_dtau) fptr(tdlsym(handle, "sys_dtau"));
  found_dtau = ((error = tdlerror()) == 0);

  tdlerror();    /* Clear any existing error */
  v_rhs = (tp_sys_rhs) fptr(tdlsym(handle, "sys_rhs"));
  found_rhs = ((error = tdlerror()) == 0);

  tdlerror();    /* Clear any existing error */
  v_deri = (tp_sys_deri) fptr(tdlsym(handle, "sys_deri"));
  found_deri = ((error = tdlerror()) == 0);

  tdlerror();    /* Clear any existing error */
  v_stpar = (tp_sys_stpar) fptr(tdlsym(handle, "sys_stpar"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_stpar(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_stsol = (tp_sys_stsol) fptr(tdlsym(handle, "sys_stsol"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_stsol(): ", error, ".");

  /* Vectorized versions */
  tdlerror();    /* Clear any existing error */
  v_p_tau = (tp_sys_p_tau) fptr(tdlsym(handle, "sys_p_tau"));
  found_p_tau = ((error = tdlerror()) == 0);
  if (!found_tau && !found_p_tau) P_MESSAGE3("Cannot find either sys_tau() or sys_p_tau(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_p_dtau = (tp_sys_p_dtau) fptr(tdlsym(handle, "sys_p_dtau"));
  found_p_dtau = ((error = tdlerror()) == 0);
  if (!found_dtau && !found_p_dtau) P_MESSAGE3("Cannot find either sys_dtau() or sys_p_dtau(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_p_rhs = (tp_sys_p_rhs) fptr(tdlsym(handle, "sys_p_rhs"));
  found_p_rhs = ((error = tdlerror()) == 0);
  if (!found_rhs && !found_p_rhs) P_MESSAGE3("Cannot find either sys_rhs() or sys_p_rhs(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_p_deri = (tp_sys_p_deri) fptr(tdlsym(handle, "sys_p_deri"));
  found_p_deri = ((error = tdlerror()) == 0);
  if (!found_deri && !found_p_deri) P_MESSAGE3("Cannot find either sys_deri() or sys_p_deri(): ", error, ".");

  nderi = (*v_nderi)();
//  std::cout<<"The order of supplied derivatives is "<<nderi<<".\n";

  f.Init(ndim()), f_eps.Init(ndim());
  f2.Init(ndim()), f_eps2.Init(ndim());
  xx_eps.Init(ndim(), 2 * ntau() + 1);
  xx_eps2.Init(ndim(), 2 * ntau() + 1);
  par_eps.Init(npar() + ParEnd);
  dxx2.Init(ndim(), ndim()), dxx_eps2.Init(ndim(), ndim());
  vt.Init(ndim());

  p_size = 0;
  p_resize(1);
}

System::~System()
{
  if (handle != 0) tdlclose(handle);
}

static inline void AX(Vector & res, const Matrix& M, const Vector& v)
{
  for (int i = 0; i < M.Row(); i++)
  {
    for (int j = 0; j < M.Col(); j++)
    {
      res(i) = 0.0;
      res(i) += M(i, j) * v(j);
    }
  }
}

void System::compileSystem(const std::string& cxxfile, const std::string& shobj, const std::string& executableDir)
{
#ifndef WIN32
  std::string compiler(CMAKE_CXX_COMPILER);
#else
  std::string compiler("g++");
#endif
  std::string::size_type slashpos = compiler.find_last_of('/');
  // truncate the string
  if (slashpos != std::string::npos) compiler = compiler.substr(slashpos+1);
  std::string cmdline(compiler);
  cmdline += " " CMAKE_CXX_FLAGS " " CMAKE_SHARED_LIBRARY_C_FLAGS " " 
              CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS " -I\"" + executableDir;
  cmdline += "/../" KNUT_INCLUDE_DIR "\"";
  cmdline += " \"" + cxxfile + "\" -o \"" + shobj + "\" 2>&1";
  
//   std::cout << "The command line: " << cmdline << "\n";
  // want to see the results if possible...
  FILE* fd = popen(cmdline.c_str(),"r");
  P_ERROR_X3(fd != 0, "Pipe cannot be opened for '", cmdline, "'.");
  std::string result;
  char* buf = new char[1024];
  char* str = 0;
  do{
    str = std::fgets(buf, 1024, fd);
    if (str != 0) result.append(str);
  } while (str != 0);
  delete[] buf;
  int cres = pclose(fd);
  if (cres != 0)  P_MESSAGE4("The output of the compile command '", cmdline, "' is ", result);
}

// Static member. Compile if necessary
void System::makeSystem(const std::string& shobj, const std::string& executableDir)
{
  // convert the name first from .so to .cpp
  std::string cxxfile(shobj);
  if (cxxfile.substr(cxxfile.size()-3,cxxfile.size()) == ".so")
  {
    cxxfile.erase(cxxfile.size()-3); cxxfile.append(".cpp");
//     std::cout << cxxfile << "\n";
  } else {
    P_MESSAGE3("The file name '", shobj, "' does not have the '.so' extension.");
  }
  // It is not portable to Windows!!!
  struct stat *sbuf_so  = new struct stat;
  struct stat *sbuf_cxx = new struct stat;
  int res_so = stat(shobj.c_str(), sbuf_so);
  int res_cxx = stat(cxxfile.c_str(), sbuf_cxx);
  // if there's no .so, but there's a .cpp
  bool compile = (res_so != 0)&&(res_cxx == 0);
  // if both .so and .cpp exist, the date decides
#ifdef __APPLE__
  if ((res_cxx == 0)&&(res_cxx == 0))
    compile |= (sbuf_so->st_mtimespec.tv_sec <= sbuf_cxx->st_mtimespec.tv_sec)&&
            (sbuf_so->st_mtimespec.tv_nsec <= sbuf_cxx->st_mtimespec.tv_nsec);
#else
  // for Linux and possibly for windows
  if ((res_cxx == 0)&&(res_cxx == 0))
    compile |= (sbuf_so->st_mtime <= sbuf_cxx->st_mtime);
#endif
  delete sbuf_so;
  delete sbuf_cxx;
  
  if (compile)
  {
    compileSystem(cxxfile, shobj, executableDir);
  }
  else if (res_so != 0)
  {
    P_MESSAGE5("The file '", shobj ,
      "' doesn't exist and it cannot be compiled from '", cxxfile, "'.");
  }
}

void System::discrderi(Matrix &out, double t, const Matrix& xx, const Vector& par,
                       int nx, const int* vx, int np, const int* vp, const Matrix& vv)
{
  const double abs_eps_x1 = 1e-6;
  const double rel_eps_x1 = 1e-6;
  const double abs_eps_p1 = 1e-6;
  const double rel_eps_p1 = 1e-6;
  const double abs_eps_x2 = 1e-6;
  const double rel_eps_x2 = 1e-6;

  const int n = ndim();
  const int m = 2 * ntau() + 1;

  // f, f_eps, xx_eps
  // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
  if ((nx == 1) && (np == 0))
  {
    rhs(f, t, xx, par);
    for (int j = 0; j < n; j++)
    {
      // xx_eps = xx;
      for (int q = 0; q < m; q++) for (int p = 0; p < n; p++) xx_eps(p, q) = xx(p, q);
      const double eps = abs_eps_x1 + rel_eps_x1 * fabs(xx(j, vx[0]));
      xx_eps(j, vx[0]) = xx(j, vx[0]) + eps;
      rhs(f_eps, t, xx_eps, par);
      for (int p = 0; p < n; p++) out(p, j) = (f_eps(p) - f(p)) / eps;
    }
  }
  // f, f_eps, par_eps
  // derivatives w.r.t. the parameters, purely, so this results a vector
  if ((nx == 0) && (np == 1))
  {
    rhs(f, t, xx, par);
    par_eps = par;
    const double eps = abs_eps_p1 + rel_eps_p1 * fabs(par(vp[0]));
    par_eps(vp[0]) = par(vp[0]) + eps;
    rhs(f_eps, t, xx, par_eps);
    for (int p = 0; p < n; p++) out(p) = (f_eps(p) - f(p)) / eps;
  }
  // f2, f_eps2, dxx2, dxx_eps2, xx_eps2, vt
  // second derivatives w.r.t. x
  if ((nx == 2) && (np == 0))
  {
    for (int j = 0; j < n; j++)
    {
      deri(dxx2, t, xx, par, 1, &vx[0], 0, vp, vv);
      xx_eps2 = xx;
      const double eps2 = abs_eps_x2 + rel_eps_x2 * fabs(xx(j, vx[1]));
      xx_eps2(j, vx[1]) += eps2;
      deri(dxx_eps2, t, xx_eps2, par, 1, &vx[0], 0, vp, vv);
      for (int p = 0; p < n; p++)
      {
        out(p, j) = 0.0;
        for (int q = 0; q < n; q++)
        {
          out(p, j) += (dxx_eps2(p, q) - dxx2(p, q)) * vv(q, vx[0]);
        }
        out(p, j) /= eps2;
      }
    }
  }
  // mixed derivative w.r.t. x and the parameters
  if ((nx == 1) && (np == 1))
  {
    deri(dxx2, t, xx, par, 1, vx, 0, vp, vv);
    par_eps = par;
    const double eps = abs_eps_p1 + rel_eps_p1 * fabs(par(vp[0]));
    par_eps(vp[0]) = par(vp[0]) + eps;
    deri(dxx_eps2, t, xx, par_eps, 1, vx, 0, vp, vv);
    for (int p = 0; p < n; p++)
    {
      for (int q = 0; q < n; q++)
      {
        out(p, q) = (dxx_eps2(p, q) - dxx2(p, q)) / eps;
      }
    }
  }
}

void System::p_discrderi( Array3D<double>& out, const Array1D<double>& time, const Array3D<double>& p_xx, const Vector& par, 
                          int nx, const int* vx, int np, const int* vp, const Array3D<double>& p_vv )
{
  const double abs_eps_x1 = 1e-6;
  const double rel_eps_x1 = 1e-6;
  const double abs_eps_p1 = 1e-6;
  const double rel_eps_p1 = 1e-6;
  const double abs_eps_x2 = 1e-6;
  const double rel_eps_x2 = 1e-6;

  const int n = ndim();
  const int m = 2 * ntau() + 1;

  p_resize(time.Size());
  // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
  if ((nx == 1) && (np == 0))
  {
//     std::cout<<"@";
    p_rhs(p_fx, time, p_xx, par);
    for (int j = 0; j < n; j++)
    {
      for (int p = 0; p < n; p++) 
        for (int q = 0; q < m; q++) 
          for (int idx=0; idx < time.Size(); ++idx) p_xx_eps(p,q,idx) = p_xx(p,q,idx);
      for (int idx=0; idx < time.Size(); ++idx)
      {
        const double eps = abs_eps_x1 + rel_eps_x1 * fabs(p_xx(j, vx[0],idx));
        p_xx_eps(j, vx[0],idx) = p_xx(j, vx[0],idx) + eps;
      }
      p_rhs(p_fx_eps, time, p_xx_eps, par);
      for (int p = 0; p < n; p++)
      {
        for (int idx=0; idx < time.Size(); ++idx)
          out(p, j, idx) = (p_fx_eps(p, idx) - p_fx(p, idx)) / (p_xx_eps(j, vx[0],idx) - p_xx(j, vx[0],idx));
      }
    }
  }
  // derivatives w.r.t. the parameters, purely, so this results a vector
  if ((nx == 0) && (np == 1))
  {
//     std::cout<<"&";
    par_eps = par;
    p_rhs(p_fx, time, p_xx, par);
    const double eps = abs_eps_p1 + rel_eps_p1 * fabs(par(vp[0]));
    par_eps(vp[0]) = par(vp[0]) + eps;
    p_rhs(p_fx_eps, time, p_xx, par_eps);
    for (int p = 0; p < n; p++)
    {
      for (int idx=0; idx < time.Size(); ++idx) out(p, 0, idx) = (p_fx_eps(p, idx) - p_fx(p, idx)) / eps;
    }
  }
  // second derivatives w.r.t. x
  if ((nx == 2) && (np == 0))
  {
//     std::cout<<"<?";
    for (int j = 0; j < n; j++)
    {
      p_deri(p2_dfx, time, p_xx, par, 1, &vx[0], 0, vp, p_vv);
      for (int p = 0; p < n; p++) 
        for (int q = 0; q < m; q++) 
          for (int idx=0; idx < time.Size(); ++idx) p2_xx_eps(p,q,idx) = p_xx(p,q,idx);
      for (int idx=0; idx < time.Size(); ++idx)
      {
        const double eps = abs_eps_x2 + rel_eps_x2 * fabs(p_xx(j, vx[1], idx));
        p2_xx_eps(j, vx[1], idx) += eps;
      }
      p_deri(p2_dfx_eps, time, p2_xx_eps, par, 1, &vx[0], 0, vp, p_vv);
      for (int p = 0; p < n; p++)
      {
        for (int idx=0; idx < time.Size(); ++idx) out(p, j, idx) = 0.0;
        for (int q = 0; q < n; q++)
        {
          for (int idx=0; idx < time.Size(); ++idx) out(p, j, idx) += (p2_dfx_eps(p, q, idx) - p2_dfx(p, q, idx)) * p_vv(q, vx[0], idx);
        }
        for (int idx=0; idx < time.Size(); ++idx) out(p, j, idx) /= (p2_xx_eps(j, vx[1], idx) - p_xx(j, vx[1], idx));
      }
    }
//     std::cout<<">";
  }
  // mixed derivative w.r.t. x and the parameters
  if ((nx == 1) && (np == 1))
  {
//     std::cout<<"<!";
    par_eps = par;
    p_deri(p2_dfx, time, p_xx, par, 1, vx, 0, vp, p_vv);
    const double eps = abs_eps_p1 + rel_eps_p1 * fabs(par(vp[0]));
    par_eps(vp[0]) = par(vp[0]) + eps;
    p_deri(p2_dfx_eps, time, p_xx, par_eps, 1, vx, 0, vp, p_vv);
    for (int p = 0; p < n; p++)
    {
      for (int q = 0; q < n; q++)
      {
        for (int idx=0; idx < time.Size(); ++idx) out(p, q, idx) = (p2_dfx_eps(p, q, idx) - p2_dfx(p, q, idx)) / eps;
      }
    }
//     std::cout<<">";
  }
}
