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
#ifndef WIN32
#include <sys/wait.h>
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#else
#include <windows.h>
#endif 
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

// Separates the string into arguments and add to the list of arguments.
// It does not work if the argument has a space inside.
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

#ifdef WIN32
#define DIRSEP '\\'
#else
#define DIRSEP '/'
#endif

// Constructs the command line and the argument list of the GNU C++ compiler (g++).
static inline void mkArgListCommandLine(std::list<std::string>& arglist, std::string& cmdline, const std::string& shobj, const std::string& executableDir)
{
#ifdef WIN32
  std::string cxxcomp("g++");
  arglist.push_back(cxxcomp);
#else
  std::string cxxcomp(CMAKE_CXX_COMPILER);
  arglist.push_back(cxxcomp.substr(cxxcomp.find_last_of(DIRSEP)+1,std::string::npos));
#endif
  // constructing the command line
  addArgList(arglist, std::string( 
    CMAKE_CXX_FLAGS " "
    CMAKE_SHARED_LIBRARY_C_FLAGS " " 
    CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS " "
    KNUT_NO_EXCEPTIONS));
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

  cmdline.erase();
  for(std::list<std::string>::const_iterator it=arglist.begin(); it != arglist.end(); ++it)
  {
    if (it != arglist.begin()) cmdline += ' ';
#ifdef WIN32
    if (it != arglist.begin()) cmdline += "\"";
    cmdline += *it;
    if (it != arglist.begin()) cmdline += "\"";
#else
    cmdline += *it;
#endif
  }
}

#ifdef WIN32

// Run the compiler with cxxstring as the program to compile.
// Catch the standard input, output and error streams to interact.
// This is the Windows version.
static void runCompiler(const std::string& cxxstring, const std::string& shobj, const std::string& executableDir)
{
  // setting up the command line
  std::list<std::string> arglist;
  std::string cmdline;
  mkArgListCommandLine(arglist, cmdline, shobj, executableDir);
  
  HANDLE g_hChildStd_IN_Rd = NULL;
  HANDLE g_hChildStd_IN_Wr = NULL;
  HANDLE g_hChildStd_OUT_Rd = NULL;
  HANDLE g_hChildStd_OUT_Wr = NULL;
  HANDLE g_hChildStd_ERR_Rd = NULL;
  HANDLE g_hChildStd_ERR_Wr = NULL;

  HANDLE g_hInputFile = NULL;

  SECURITY_ATTRIBUTES saAttr;
 
  // Set the bInheritHandle flag so pipe handles are inherited. 
  saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
  saAttr.bInheritHandle = TRUE;
  saAttr.lpSecurityDescriptor = NULL;

  // Create a pipe for the child process's STDOUT. 
  if (!CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &saAttr, 0)) 
    P_MESSAGE2("StdoutRd CreatePipe", (int)GetLastError());
  // Ensure the read handle to the pipe for STDOUT is not inherited.
  if (!SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0))
    P_MESSAGE2("Stdout SetHandleInformation", (int)GetLastError());
  // Create a pipe for the child process's STDIN. 
  if (!CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &saAttr, 0))
    P_MESSAGE2("Stdin CreatePipe", (int)GetLastError());
  // Ensure the write handle to the pipe for STDIN is not inherited.
  if (!SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0))
    P_MESSAGE2("Stdin SetHandleInformation", (int)GetLastError());
  // Create a pipe for the child process's STDERR.
  if (!CreatePipe(&g_hChildStd_ERR_Rd, &g_hChildStd_ERR_Wr, &saAttr, 0)) 
    P_MESSAGE2("Stdin CreatePipe", (int)GetLastError());
  // Ensure the read handle to the pipe for STDERR is not inherited.
  if (!SetHandleInformation(g_hChildStd_ERR_Rd, HANDLE_FLAG_INHERIT, 0))
    P_MESSAGE2("Stdin SetHandleInformation", (int)GetLastError());
 
   // Create the child process. 
  CHAR *szCmdline = new CHAR[cmdline.size()+1];
  strcpy(szCmdline, cmdline.c_str());
  
  PROCESS_INFORMATION piProcInfo; 
  STARTUPINFO siStartInfo;
  BOOL bSuccess = FALSE;
 
  // Set up members of the PROCESS_INFORMATION structure. 
  ZeroMemory( &piProcInfo, sizeof(PROCESS_INFORMATION) );
 
  // Set up members of the STARTUPINFO structure. 
  // This structure specifies the STDIN and STDOUT handles for redirection.
  ZeroMemory( &siStartInfo, sizeof(STARTUPINFO) );
  siStartInfo.cb = sizeof(STARTUPINFO); 
  siStartInfo.hStdError = g_hChildStd_ERR_Wr;
  siStartInfo.hStdOutput = g_hChildStd_OUT_Wr;
  siStartInfo.hStdInput = g_hChildStd_IN_Rd;
  siStartInfo.dwFlags |= STARTF_USESTDHANDLES;
 
  // Create the child process. 
  bSuccess = CreateProcess(NULL, 
    szCmdline,     // command line 
    NULL,          // process security attributes 
    NULL,          // primary thread security attributes 
    TRUE,          // handles are inherited 
    0,             // creation flags 
    NULL,          // use parent's environment 
    NULL,          // use parent's current directory 
    &siStartInfo,  // STARTUPINFO pointer 
    &piProcInfo);  // receives PROCESS_INFORMATION 
   
  // If an error occurs, exit the application. 
  if ( ! bSuccess )
    P_MESSAGE2("CreateProcess", (int)GetLastError());
 
  // Write to the pipe that is the standard input for a child process. 
  // Data is written to the pipe's buffers, so it is not necessary to wait
  // until the child process is running before writing data.
  DWORD dwRead = cxxstring.size(), dwWritten;
  bSuccess = WriteFile(g_hChildStd_IN_Wr, cxxstring.c_str(), dwRead, &dwWritten, NULL);
  if ( ! bSuccess ) 
    P_MESSAGE2("WriteFile", (int)GetLastError());

  if (!CloseHandle(g_hChildStd_IN_Wr)) 
    P_MESSAGE2("StdInWr CloseHandle", (int)GetLastError());
  if (!CloseHandle(g_hChildStd_IN_Rd)) 
    P_MESSAGE2("StdInRd CloseHandle", (int)GetLastError());
 
  // Read from pipe that is the standard output for child process.
  const size_t bufsize = 2048;
  CHAR *out_buf = new char[bufsize];
  CHAR *err_buf = new char[bufsize];
  dwRead = 0;
  dwWritten = 0;
  bSuccess = FALSE;
  HANDLE hParentStdOut = GetStdHandle(STD_OUTPUT_HANDLE);

  // Close the write end of the pipe before reading from the 
  // read end of the pipe, to control child process execution.
  // The pipe is assumed to have enough buffer space to hold the
  // data the child process has already written to it.
  if (!CloseHandle(g_hChildStd_OUT_Wr)) 
    P_MESSAGE2("StdOutWr CloseHandle", (int)GetLastError());
   
  std::string outstr;
  std::string errstr;
  for (;;) 
  {
    bSuccess = ReadFile( g_hChildStd_OUT_Rd, out_buf, bufsize-1, &dwRead, NULL);
    err_buf[dwRead] = '\0';
    if( ! bSuccess || dwRead == 0 ) break;
    outstr.append(out_buf);
  }
  
  if (!CloseHandle(g_hChildStd_OUT_Rd)) 
    P_MESSAGE2("StdOutRd CloseHandle", (int)GetLastError());
  
  if (!CloseHandle(g_hChildStd_ERR_Wr)) 
    P_MESSAGE2("StdErrWr CloseHandle", (int)GetLastError());
   
  for (;;) 
  {
    bSuccess = ReadFile( g_hChildStd_ERR_Rd, err_buf, bufsize-1, &dwRead, NULL);
    err_buf[dwRead] = '\0';
    if( ! bSuccess || dwRead == 0 ) break;
    errstr.append(err_buf);
  }   
  
  if (!CloseHandle(g_hChildStd_ERR_Rd)) 
    P_MESSAGE2("StdErrRd CloseHandle", (int)GetLastError());
  
  WaitForSingleObject(piProcInfo.hProcess, INFINITE);
  DWORD exitStatus;
  GetExitCodeProcess(piProcInfo.hProcess, &exitStatus); 
  CloseHandle(piProcInfo.hProcess);
  CloseHandle(piProcInfo.hThread);

  if (exitStatus != 0) P_MESSAGE6("The error output of the compile command '",
      cmdline, "' is ", errstr.c_str(), " and the standard output is ", outstr.c_str());
  
  delete[] out_buf;
  delete[] err_buf;
}

#else

// Opens the three pipes for interacting with a program.
// arglist is the argument list starting with the program name.
// the int* - s are the file numbers for the three pipes.
int pipeOpen(std::list<std::string>& arglist, int* input, int* output, int* error)
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
    if (input)  
    if (output) 
    if (error)  
    if (input)
    {
   	  if (close (fdc_input[0]) == -1) P_MESSAGE2("close ", strerror(errno));
   	  int iflags = fcntl(fdc_input[1], F_GETFL, 0);
      iflags |= O_NONBLOCK;
      fcntl(fdc_input[1], F_SETFL, iflags);
   	  *input = fdc_input[1];
    }   	  
    if (output)
    {
   	  if (close (fdc_output[1]) == -1) P_MESSAGE2("close ", strerror(errno));
      int oflags = fcntl(fdc_output[0], F_GETFL, 0);
      oflags |= O_NONBLOCK;
      fcntl(fdc_output[0], F_SETFL, oflags);
      *output = fdc_output[0];
    }
    if (error)
    {
      if (close (fdc_error[1]) == -1) P_MESSAGE2("close ", strerror(errno));
      int eflags = fcntl(fdc_error[0], F_GETFL, 0);
      eflags |= O_NONBLOCK;
      fcntl(fdc_error[0], F_SETFL, eflags);
      *error = fdc_error[0];
    }
  }
  for (size_t k=0; k<arglist.size(); ++k) delete[] argv[k];
  return pid;
}

// Run the compiler with cxxstring as the program to compile.
// Catch the standard input, output and error streams to interact.
// This is the Linux/Mac OS X version.
static void runCompiler(const std::string& cxxstring, const std::string& shobj, const std::string& executableDir)
{
  std::list<std::string> arglist;
  std::string cmdline;
  mkArgListCommandLine(arglist, cmdline, shobj, executableDir);
  
  //   std::cout << cmdline << std::endl;
  // running the command
  int input, output, error;
  int pid = pipeOpen(arglist, &input, &output, &error);
  // we need to write and read until all of them are finished
  bool outfin = false, infin = false, errfin = false;
  // input
  const char* cxxbuf = cxxstring.c_str();
  ssize_t cxxlen = cxxstring.size(), wbytes = 0;
  // output and err
  pollfd fds[3] = {{input, POLLOUT, 0},{output, POLLIN, 0},{error, POLLIN, 0}};
  const size_t bufsize = 2048;
  char *out_buf = new char[bufsize];
  std::string out_str, err_str;
  int rct = 0, ect = 0;
  do {
    int pact = poll(fds, 3, 1000);
//    std::cerr << "Pipes ready: " << pact 
//    	<< " EV0 " <<  fds[0].revents 
//    	<< " EV1 " <<  fds[1].revents
//    	<< " EV2 " <<  fds[2].revents << " " << POLLOUT << " " << POLLIN <<"\n";
    P_ERROR_X2(pact != -1, "Error polling pipe status: ", strerror(errno));
  	if (!infin && (fds[0].revents != 0))
  	{
      ssize_t obytes = write (input, cxxbuf+wbytes, cxxlen-wbytes);
//      std::cerr << "runCompiler: Written bytes " << obytes << " out of " << cxxlen << "\n";
      if (obytes == -1) P_MESSAGE2("Error feeding the compiler: ", strerror(errno));
      else wbytes += obytes;
	  if (wbytes == cxxlen) 
	  { 
      infin = true;
      fds[0].events = 0;
      if (close (input) == -1) P_MESSAGE2("Error closing compiler input pipe: ", strerror(errno));
	  }
  	}
  	if (!outfin && (fds[1].revents != 0))
  	{
      ssize_t bytes = read(output, out_buf, bufsize-1);
      ++rct;
      if (bytes > 0)
      {
//        std::cerr << "runCompiler: stdout bytes " << bytes << "\n";
      	out_buf[bytes] = '\0'; 
      	out_str.append(out_buf);
      } 
      else if (bytes == 0)
      {
        outfin = true; // finished
        fds[1].events = 0;
      	if (close (output) == -1) P_MESSAGE2("Error closing output pipe: ", strerror(errno));
      }
      else if ((bytes == -1) && (errno != EAGAIN)) P_MESSAGE2("Error reading standard output: ", strerror(errno));
    }
  	if (!errfin && (fds[2].revents != 0))
  	{
      ssize_t bytes = read(error, out_buf, bufsize-1);
      ++ect;
      if (bytes > 0)
      {
//        std::cerr << "runCompiler: stderr bytes " << bytes << "\n";
      	out_buf[bytes] = '\0'; 
      	err_str.append(out_buf);
      } 
      else if (bytes == 0) 
      {
      	errfin = true; // finished
	  	fds[2].events = 0;
      	if (close (error) == -1) P_MESSAGE2("Error closing error pipe: ", strerror(errno));
      }
      else if ((bytes == -1) && (errno != EAGAIN)) P_MESSAGE2("Error reading standard error: ", strerror(errno));
  	}
  } while (!outfin || !infin || !errfin);
//  std::cerr << "Waited for read " << rct << " and error " << ect << " cycles."; 
  // checking status of compiler
  int status = 0;
  pid_t retval = waitpid(-1, &status, 0);
  if (retval == -1) P_MESSAGE2("Error while waiting for the compiler: ", strerror(errno));
  if (retval == 0) P_MESSAGE1("A mysterious error. (No child wanted to report status.)");
  // Check the exist status
  if (WIFEXITED(status))
  {
    if (WEXITSTATUS(status) != 0) P_MESSAGE8("Error code is ", WEXITSTATUS(status), " and the error output of the compile command '",
      cmdline, "' is\n", err_str.c_str(), "\nand the standard output is\n", out_str.c_str());
//    std::ofstream cxxfile("error.cxx");
//    cxxfile << cxxbuf;
  } else
  {
    P_MESSAGE1("The compiler was unexpectedly terminated.");
  }
}

#endif

KNSystem::KNSystem(const std::string& shobj, int usederi)
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
  std::string::size_type sidx = executableDir.find_last_of(DIRSEP);
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
  v_mass = (tp_sys_mass) fptr(tdlsym(handle, "sys_mass"));
  if ((error = tdlerror()) != 0) v_mass = 0;
  
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

  nderi = std::min ( (*v_nderi)(), usederi );
  
  f.init(ndim()), f_eps.init(ndim());
  f2.init(ndim()), f_eps2.init(ndim());
  xx_eps.init(ndim(), 2 * ntau() + 1);
  xx_eps2.init(ndim(), 2 * ntau() + 1);
  par_eps.init(VarToIndex(VarEnd,npar()));
  dxx2.init(ndim(), ndim()), dxx_eps2.init(ndim(), ndim());
  vt.init(ndim());

  p_size = 0;
  p_resize(1);
}

KNSystem::~KNSystem()
{
  if (handle != 0) tdlclose(handle);
}

static inline void timesX(KNVector & res, const KNMatrix& M, const KNVector& v)
{
  for (int i = 0; i < M.row(); i++)
  {
    for (int j = 0; j < M.col(); j++)
    {
      res(i) = 0.0;
      res(i) += M(i, j) * v(j);
    }
  }
}

void KNSystem::compileSystem(const std::string& cxxfile, const std::string& shobj, const std::string& executableDir)
{
  std::ifstream file(cxxfile.c_str());
  std::string cxxcode, buf;
  while (!std::getline(file,buf).eof())
  {
    cxxcode += buf + '\n';
  }
  runCompiler(cxxcode, shobj, executableDir);
}

void KNSystem::generateSystem(const std::string& vffile, const std::string& executableDir)
{
#ifdef GINAC_FOUND
  // parsing the vector field file
  std::ostringstream cxxcode;
  std::map<std::string, std::string> options;
  VectorField vf;
  vf.ReadXML(vffile);
  int pserr = vf.ProcessSymbols();
  if (pserr == -1) P_MESSAGE1("Could not parse the vector filed definition.");
  if (vf.isStateDependentDelay()) P_MESSAGE1("State dependent delays are not suported yet.");
  vf.PrintKnut(cxxcode, options);
  runCompiler(cxxcode.str(), vffile + ".so", executableDir);
#endif
}

// Static member. Compile if necessary
void KNSystem::makeSystem(const std::string& shobj, const std::string& executableDir)
{
  // convert the name first from .so to .cpp
  std::string cxxfile(shobj);
  std::string vffile;
  if (cxxfile.substr(std::max((int)cxxfile.size()-3,0),cxxfile.size()) == ".so")
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

void KNSystem::p_discrderi( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& p_xx, const KNVector& par, int sel, 
                          int nx, const int* vx, int np, const int* vp, const KNArray3D<double>& p_vv )
{
  const double abs_eps_x1 = 1e-6;
  const double rel_eps_x1 = 1e-6;
  const double abs_eps_p1 = 1e-6;
  const double rel_eps_p1 = 1e-6;
  const double abs_eps_x2 = 1e-6;
  const double rel_eps_x2 = 1e-6;

  const int n = ndim();
  const int m = 2 * ntau() + 1;

  p_resize(time.size());
  // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
  if ((nx == 1) && (np == 0))
  {
//     std::cout<<"@";
    p_rhs(p_fx, time, p_xx, par, sel);
    for (int j = 0; j < n; j++)
    {
      for (int p = 0; p < n; p++) 
          for (int idx=0; idx < time.size(); ++idx) p_xx_eps(p,vx[0],idx) = p_xx(p,vx[0],idx);
      for (int idx=0; idx < time.size(); ++idx)
      {
        const double eps = abs_eps_x1 + rel_eps_x1 * fabs(p_xx(j, vx[0],idx));
        p_xx_eps(j, vx[0],idx) = p_xx(j, vx[0],idx) + eps;
      }
      p_rhs(p_fx_eps, time, p_xx_eps, par, sel);
      for (int p = 0; p < n; p++)
      {
        for (int idx=0; idx < time.size(); ++idx)
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
      for (int idx=0; idx < time.size(); ++idx) out(p, 0, idx) = (p_fx_eps(p, idx) - p_fx(p, idx)) / eps;
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
          for (int idx=0; idx < time.size(); ++idx) p2_xx_eps(p,q,idx) = p_xx(p,q,idx);
      for (int idx=0; idx < time.size(); ++idx)
      {
        const double eps = abs_eps_x2 + rel_eps_x2 * fabs(p_xx(j, vx[1], idx));
        p2_xx_eps(j, vx[1], idx) += eps;
      }
      p_deri(p2_dfx_eps, time, p2_xx_eps, par, sel, 1, &vx[0], 0, vp, p_vv);
      for (int p = 0; p < n; p++)
      {
        for (int idx=0; idx < time.size(); ++idx) out(p, j, idx) = 0.0;
        for (int q = 0; q < n; q++)
        {
          for (int idx=0; idx < time.size(); ++idx) out(p, j, idx) += (p2_dfx_eps(p, q, idx) - p2_dfx(p, q, idx)) * p_vv(q, vx[0], idx);
        }
        for (int idx=0; idx < time.size(); ++idx) out(p, j, idx) /= (p2_xx_eps(j, vx[1], idx) - p_xx(j, vx[1], idx));
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
        for (int idx=0; idx < time.size(); ++idx) out(p, q, idx) = (p2_dfx_eps(p, q, idx) - p2_dfx(p, q, idx)) / eps;
      }
    }
//     std::cout<<">";
  }
}
