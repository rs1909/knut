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
#include <cstring>
#include <list>
#include <map>
#include <sstream>
#include <fstream>
extern "C"{
#include <sys/stat.h>
#include <sys/wait.h>
#include <errno.h>
}

#include "config.h"
#ifdef GINAC_FOUND
#include "vf.h"
#endif
#include "system.h"
#include "matrix.h"
#include "pointtype.h"

#ifdef __APPLE__
#include <CoreFoundation/CFBase.h>
#include <CoreFoundation/CFString.h>
#include <CoreFoundation/CFBundle.h>
#endif

int pipeOpen(std::list<std::string>& arglist,
                           int* input, int* output, int* error)
{
  if ( (input==0)&&(output==0)&&(error==0) ) return -1;
  char *argv[arglist.size()+1];
  int i = 0;
  for (std::list<std::string>::const_iterator it=arglist.begin(); it != arglist.end(); ++it)
  {
    argv[i] = new char[it->size()+1];
    strcpy(argv[i],it->c_str());
    ++i;
  }
  argv[i] = 0;

  int fds_output[2], fds_input[2], fds_error[2];
  int fdc_output[2], fdc_input[2], fdc_error[2];

  /* Create a pipe.  File descriptors for the two ends of the pipe are
     placed in fds.  */
  if (output) if (pipe (fds_output) == -1) P_MESSAGE2("pipe: ", strerror(errno));
  if (input)  if (pipe (fds_input) == -1) P_MESSAGE2("pipe: ", strerror(errno));
  if (error)  if (pipe (fds_error) == -1) P_MESSAGE2("pipe: ", strerror(errno));
  for (int i=0; i<2; ++i)
  {
    fdc_output[i] = fds_output[i]; 
    fdc_input[i] = fds_input[i];
    fdc_error[i] = fds_error[i];
  }
  /* Fork a child process.  */
  pid_t pid = fork ();
  if (pid == (pid_t) 0)
  {
    /* This is the child process.  Close our copy of the read (write) end of
       the file descriptor.  */
    if (input)  if (close (fds_input[1]) == -1) P_MESSAGE2("close: ", strerror(errno));
    if (output) if (close (fds_output[0]) == -1) P_MESSAGE2("close: ", strerror(errno));
    if (error)  if (close (fds_error[0]) == -1) P_MESSAGE2("close ", strerror(errno));
    /* Connect the read(write) end of the pipe to standard input.  */
    if (input)  if (dup2 (fds_input[0], STDIN_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
    if (output) if (dup2 (fds_output[1], STDOUT_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
    if (error)  if (dup2 (fds_error[1], STDERR_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
    /* Replace the child process with the "cat" program.  */
    int st = execvp (argv[0], argv);
    // This wouldn't return unless there's a problem
    if (st == -1) P_MESSAGE4("Error executing command ", argv[0], ": ", strerror(errno));
  }
  else
  {
    /* Close our copy of the write (read) end of the file descriptor.  */
    if (input)  if (close (fdc_input[0]) == -1) P_MESSAGE2("close ", strerror(errno));
    if (output) if (close (fdc_output[1]) == -1) P_MESSAGE2("close ", strerror(errno));
    if (error)  if (close (fdc_error[1]) == -1) P_MESSAGE2("close ", strerror(errno));
    if (input)  *input = fdc_input[1];
    if (output) *output = fdc_output[0];
    if (error)  *error = fdc_error[0];
  }
  for (int k=0; k<1; ++k) delete[] argv[k];
  return pid;
}

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
#elif WIN32
  char *buf = new char[MAX_PATH]; //always use MAX_PATH for filepaths
  GetModuleFileName(NULL, buf, MAX_PATH*sizeof(char));
  executableFile = buf;
  delete[] buf;
#endif

  std::string executableDir(executableFile);
#ifndef WIN32
  std::string::size_type sidx = executableDir.find_last_of('/');
#else
  std::string::size_type sidx = executableDir.find_last_of('\\');
#endif
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
  v_nevent = (tp_sys_nevent) fptr(tdlsym(handle, "sys_nevent"));
  if ((error = tdlerror()) != 0) v_nevent = 0;

  tdlerror();    /* Clear any existing error */
  v_tau = (tp_sys_tau) fptr(tdlsym(handle, "sys_tau"));
  if ((error = tdlerror()) != 0) v_tau = 0;

  tdlerror();    /* Clear any existing error */
  v_dtau = (tp_sys_dtau) fptr(tdlsym(handle, "sys_dtau"));
  if ((error = tdlerror()) != 0) v_dtau = 0;
  
  tdlerror();    /* Clear any existing error */
  v_rhs = (tp_sys_rhs) fptr(tdlsym(handle, "sys_rhs"));
  if ((error = tdlerror()) != 0) v_rhs = 0;

  tdlerror();    /* Clear any existing error */
  v_deri = (tp_sys_deri) fptr(tdlsym(handle, "sys_deri"));
  if ((error = tdlerror()) != 0) v_deri = 0;

  tdlerror();    /* Clear any existing error */
  v_stpar = (tp_sys_stpar) fptr(tdlsym(handle, "sys_stpar"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_stpar(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_stsol = (tp_sys_stsol) fptr(tdlsym(handle, "sys_stsol"));
  P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_stsol(): ", error, ".");
  
  tdlerror();    /* Clear any existing error */
  v_parnames = (tp_sys_parnames) fptr(tdlsym(handle, "sys_parnames"));
  if ((error = tdlerror()) != 0) v_parnames = 0;

  /* Vectorized versions */
  tdlerror();    /* Clear any existing error */
  v_p_tau = (tp_sys_p_tau) fptr(tdlsym(handle, "sys_p_tau"));
  if ((error = tdlerror()) != 0) v_p_tau = 0;
  if ((v_tau == 0) && (v_p_tau == 0)) P_MESSAGE3("Cannot find either sys_tau() or sys_p_tau(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_p_dtau = (tp_sys_p_dtau) fptr(tdlsym(handle, "sys_p_dtau"));
  if ((error = tdlerror()) != 0) v_p_dtau = 0;
  if ((v_dtau == 0) && (v_p_dtau == 0)) P_MESSAGE3("Cannot find either sys_dtau() or sys_p_dtau(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_p_rhs = (tp_sys_p_rhs) fptr(tdlsym(handle, "sys_p_rhs"));
  if ((error = tdlerror()) != 0) v_p_rhs = 0;
  if ((v_rhs == 0) && (v_p_rhs == 0)) P_MESSAGE3("Cannot find either sys_rhs() or sys_p_rhs(): ", error, ".");

  tdlerror();    /* Clear any existing error */
  v_p_deri = (tp_sys_p_deri) fptr(tdlsym(handle, "sys_p_deri"));
  if ((error = tdlerror()) != 0) v_p_deri = 0;
  if ((v_deri == 0) && (v_p_deri == 0)) P_MESSAGE3("Cannot find either sys_deri() or sys_p_deri(): ", error, ".");

  nderi = (*v_nderi)();

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

static inline void addArgList(std::list<std::string>& argv, const std::string& str)
{
  size_t start = 0;
  size_t i=0;
  while (i < str.size())
  {    
    while (((str[i] == ' ')||(str[i] == '\t'))&&(i < str.size())) ++i;
    start = i;
    while (((str[i] != ' ')&&(str[i] != '\t'))&&(i < str.size())) ++i;
    argv.push_back(str.substr(start, i-start));
  }
}

static inline void toCommandLine(std::string& cmdline, const std::list<std::string>& arglist)
{
  cmdline.erase();
  for(std::list<std::string>::const_iterator it=arglist.begin(); it != arglist.end(); ++it)
  {
    if (it != arglist.begin()) cmdline += ' ';
    cmdline += *it;
  }
}

#if WIN32
#define DIRSEP '\\'
#else
#define DIRSEP '/'
#endif

static void runCompiler(const std::string& cxxstring, const std::string& shobj, const std::string& executableDir)
{  
  std::string cxxcomp(CMAKE_CXX_COMPILER);
  // constructing the command line
  std::list<std::string> arglist;
  arglist.push_back(cxxcomp.substr(cxxcomp.find_last_of(DIRSEP)+1,std::string::npos));
  addArgList(arglist, std::string( 
    CMAKE_CXX_FLAGS " "
    CMAKE_SHARED_LIBRARY_C_FLAGS " " 
    CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS));
  std::string includeArg("-I");
  includeArg += executableDir.substr(0,executableDir.find_last_of(DIRSEP));
  includeArg += DIRSEP;
  includeArg += KNUT_INCLUDE_DIR;
  arglist.push_back(includeArg);
  arglist.push_back("-x");
  arglist.push_back("c++");
  arglist.push_back("-o");
  arglist.push_back(shobj);
  arglist.push_back("-");

  std::string cmdline;
  toCommandLine(cmdline, arglist);
//   std::cout << cmdline << std::endl;
  // running the command
  int input, output, error;
  int pid = pipeOpen(arglist, &input, &output, &error);
  write (input, cxxstring.c_str(), cxxstring.size());
  close (input);
  int status = 0;
  pid_t retval = waitpid(-1, &status, 0);
  if (retval == -1) P_MESSAGE2("Error while waiting for the compiler: ", strerror(errno));
  if (retval == 0) P_MESSAGE1("A mysterious error. (No child wanted to report status.)");
  // Output
  const size_t bufsize = 2048;
  char *out_buf = new char[bufsize];
  char *err_buf = new char[bufsize];
  size_t bytes = read(output, out_buf, bufsize);
  out_buf[bytes] = '\0';
  // Standard Error
  bytes = read(error, err_buf, bufsize);
  err_buf[bytes] = '\0';
  // Check the exist status
  close (output);
  close (error);
  if (WIFEXITED(status))
  {
    if (WEXITSTATUS(status) != 0) P_MESSAGE6("The error output of the compile command '",
      cmdline, "' is ", err_buf, " and the standard output is ", out_buf);
  } else
  {
    P_MESSAGE1("The compiler was unexpectedly terminated.");
  }
}

void System::compileSystem(const std::string& cxxfile, const std::string& shobj, const std::string& executableDir)
{
  std::ifstream file(cxxfile.c_str());
  std::string cxxcode, buf;
  while (!std::getline(file,buf).eof())
  {
    cxxcode += buf + '\n';
  }
  runCompiler(cxxcode, shobj, executableDir);
}

void System::generateSystem(const std::string& vffile, const std::string& executableDir)
{
#ifdef GINAC_FOUND
  // parsing the vector field file
  std::ostringstream cxxcode;
  std::map<std::string, std::string> options;
  VectorField vf;
  vf.ReadXML(vffile);
  int pserr = vf.ProcessSymbols();
  if (pserr == -1) P_MESSAGE1("Could not parse the vector filed definition.");
  if (vf.testHasNonconstantDelay()) P_MESSAGE1("Nonconstant delays are not suported yet.");
  vf.PrintKnut(cxxcode, options);
  runCompiler(cxxcode.str(), vffile + ".so", executableDir);
#endif
} 

// Static member. Compile if necessary
void System::makeSystem(const std::string& shobj, const std::string& executableDir)
{
  // convert the name first from .so to .cpp
  std::string cxxfile(shobj);
  std::string vffile;
  if (cxxfile.substr(cxxfile.size()-3,cxxfile.size()) == ".so")
  {
    cxxfile.erase(cxxfile.size()-3);
    vffile = cxxfile;
    cxxfile.append(".cpp");
  } else {
    P_MESSAGE3("The file name '", shobj, "' does not have the '.so' extension.");
  }
  // It is not portable to Windows!!!
  struct stat *sbuf_so  = new struct stat;
  struct stat *sbuf_cxx = new struct stat;
  struct stat *sbuf_vf = new struct stat;
  int res_so = stat(shobj.c_str(), sbuf_so);
  int res_cxx = stat(cxxfile.c_str(), sbuf_cxx);
  int res_vf = stat(vffile.c_str(), sbuf_vf);
  // if there's no .so, but there's a .cpp
  bool compile = (res_so != 0)&&(res_cxx == 0);
  bool generate = (res_so != 0)&&(res_vf == 0);
  // if both .so and .cpp exist, the date decides
#ifdef __APPLE__
  if (res_cxx == 0)
    compile |= (sbuf_so->st_mtimespec.tv_sec <= sbuf_cxx->st_mtimespec.tv_sec)&&
            (sbuf_so->st_mtimespec.tv_nsec <= sbuf_cxx->st_mtimespec.tv_nsec);
  if ((res_vf == 0) && (vffile.substr(vffile.size()-3,vffile.size()) == ".vf"))
    generate |= (sbuf_so->st_mtimespec.tv_sec <= sbuf_vf->st_mtimespec.tv_sec)&&
            (sbuf_so->st_mtimespec.tv_nsec <= sbuf_vf->st_mtimespec.tv_nsec);
#else
  // for Linux and possibly for windows
  if (res_cxx == 0)
    compile |= (sbuf_so->st_mtime <= sbuf_cxx->st_mtime);
  if ((res_vf == 0) && (vffile.substr(vffile.size()-3,vffile.size()) == ".vf"))
    generate |= (sbuf_so->st_mtime <= sbuf_cxx->st_mtime);
#endif
  delete sbuf_so;
  delete sbuf_cxx;
  
  if (generate)
  {
    generateSystem(vffile, executableDir);
  } else
  {
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
}

void System::p_discrderi( Array3D<double>& out, const Array1D<double>& time, const Array3D<double>& p_xx, const Vector& par, int sel, 
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
    p_rhs(p_fx, time, p_xx, par, sel);
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
      p_rhs(p_fx_eps, time, p_xx_eps, par, sel);
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
    p_rhs(p_fx, time, p_xx, par, sel);
    const double eps = abs_eps_p1 + rel_eps_p1 * fabs(par(vp[0]));
    par_eps(vp[0]) = par(vp[0]) + eps;
    p_rhs(p_fx_eps, time, p_xx, par_eps, sel);
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
      p_deri(p2_dfx, time, p_xx, par, sel, 1, &vx[0], 0, vp, p_vv);
      for (int p = 0; p < n; p++) 
        for (int q = 0; q < m; q++) 
          for (int idx=0; idx < time.Size(); ++idx) p2_xx_eps(p,q,idx) = p_xx(p,q,idx);
      for (int idx=0; idx < time.Size(); ++idx)
      {
        const double eps = abs_eps_x2 + rel_eps_x2 * fabs(p_xx(j, vx[1], idx));
        p2_xx_eps(j, vx[1], idx) += eps;
      }
      p_deri(p2_dfx_eps, time, p2_xx_eps, par, sel, 1, &vx[0], 0, vp, p_vv);
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
    p_deri(p2_dfx, time, p_xx, par, sel, 1, vx, 0, vp, p_vv);
    const double eps = abs_eps_p1 + rel_eps_p1 * fabs(par(vp[0]));
    par_eps(vp[0]) = par(vp[0]) + eps;
    p_deri(p2_dfx_eps, time, p_xx, par_eps, sel, 1, vx, 0, vp, p_vv);
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
