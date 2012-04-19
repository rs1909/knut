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
#include <algorithm>
#include <limits>
extern "C"{
#include <sys/stat.h>
#ifndef WIN32
#include <sys/wait.h>
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <unistd.h>
#else
#include <windows.h>
#endif 
}

#include "system.h"
#include "matrix.h"
#include "pointtype.h"
// #include "config.h" // already included from system.h
#ifdef GINAC_FOUND
#include "vf.h"
#endif

#ifdef __APPLE__
#include <CoreFoundation/CFBase.h>
#include <CoreFoundation/CFString.h>
#include <CoreFoundation/CFBundle.h>
#endif

bool KNSystem::workingCompiler = true;

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
  arglist.push_back("-pipe");
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
  
  delete[] out_buf; out_buf = 0;
  delete[] err_buf; err_buf = 0;
  delete[] szCmdline; szCmdline = 0;
}

#else

// Opens the three pipes for interacting with a program.
// arglist is the argument list starting with the program name.
// the int* - s are the file numbers for the three pipes.
static pid_t pipeOpen(std::list<std::string>& arglist, int* input, int* output, int* error)
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
//    std::cout << "Closing ends\n";
    /* This is the child process.  Close our copy of the read (write) end of
       the file descriptor.  */
    if (input)  if (close (fds_input[1]) == -1) P_MESSAGE2("close: ", strerror(errno));
    if (output) if (close (fds_output[0]) == -1) P_MESSAGE2("close: ", strerror(errno));
    if (error)  if (close (fds_error[0]) == -1) P_MESSAGE2("close ", strerror(errno));
//    std::cout << "Closed ends\n";
    /* Connect the read(write) end of the pipe to standard input.  */
    if (input)  if (dup2 (fds_input[0], STDIN_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
//    std::cout << "Dup input\n";
    if (output) if (dup2 (fds_output[1], STDOUT_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
//    std::cout << "Dup output\n";
    if (error)  if (dup2 (fds_error[1], STDERR_FILENO) == -1) P_MESSAGE2("dup2: ", strerror(errno));
    /* Replace the child process with the "cat" program.  */
//    std::cout << "Starting process\n"; std::cout.flush();
    int st = execvp (argv[0], argv);
    // This wouldn't return unless there's a problem
    if (st == -1) P_MESSAGE4("Error executing command ", argv[0], ": ", strerror(errno));
//    std::cout << "Somehow returned\n"; std::cout.flush();
  }
  else
  {
    if (pid == -1) P_MESSAGE2("Failed to fork process: ", strerror(errno));
    /* Close our copy of the write (read) end of the file descriptor.  */
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
//    std::cout << "Somehow returned - other!!\n"; std::cout.flush();
  }
  for (size_t k=0; k<arglist.size(); ++k) { delete[] argv[k]; argv[k] = 0; }
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
  pid_t pid = pipeOpen(arglist, &input, &output, &error);
  // we need to write and read until all of them are finished
  bool outfin = false, infin = false, errfin = false;
  // input
  const char* cxxbuf = cxxstring.c_str();
  size_t cxxlen = cxxstring.size(), wbytes = 0;
  // output and err
  pollfd fds[3] = {{input, POLLOUT|POLLHUP, 0},{output, POLLIN|POLLHUP, 0},{error, POLLIN|POLLHUP, 0}};
  const size_t bufsize = 2048;
  char *out_buf = new char[bufsize];
  std::string out_str, err_str;
  int rct = 0, ect = 0;
  size_t iin = 0, iout = 1, ierr = 2, nfds = 3;
  do {
    int pact = poll(fds, nfds, 1000);
//     std::cerr << "Pipes ready: " << pact 
//      << " EV0 " <<  fds[iin].revents 
//      << " EV1 " <<  fds[iout].revents
//      << " EV2 " <<  fds[ierr].revents 
//      << " POLLOUT=" << POLLOUT << " POLLIN=" << POLLIN << " POLLPRI=" << POLLPRI 
//      << " POLLHUP=" << POLLHUP << " POLLERR=" << POLLERR 
//      << " POLLWRNORM=" << POLLWRNORM << " POLLWRBAND=" << POLLWRBAND 
//      << " POLLRDNORM=" << POLLRDNORM << " POLLRDBAND=" << POLLRDBAND << " POLLNVAL=" << POLLNVAL << "\n";
    P_ERROR_X2(pact != -1, "Error polling pipe status: ", strerror(errno));
    if (!infin)
    {
      P_ERROR_X2((fds[iin].revents & (POLLNVAL|POLLERR)) == 0, "Input pipe error -> ", err_str);
      if ((fds[iin].revents & POLLOUT) != 0)
      {
        ssize_t obytes = write (input, cxxbuf+wbytes, cxxlen-wbytes);
//        std::cerr << "runCompiler: Written bytes " << obytes << " out of " << cxxlen << "\n";
        if (obytes == -1) P_MESSAGE2("Error feeding the compiler: ", strerror(errno));
        else wbytes += static_cast<size_t>(obytes);
        if (wbytes == cxxlen) infin = true;
      }
      if ((fds[iin].revents & POLLHUP) != 0) infin = true;
      if (infin)
      {
        if (!outfin) { fds[iin] = fds[iout]; iout = iin; }
        if (!errfin) { fds[iout] = fds[ierr]; ierr = iout; }
        nfds -= 1;
        if (close (input) == -1) P_MESSAGE2("Error closing compiler input pipe: ", strerror(errno));
      }
    }
    if (!outfin)
    {
      P_ERROR_X2((fds[iout].revents & (POLLNVAL|POLLERR)) == 0, "Output pipe error -> ", err_str);
      if ((fds[iout].revents & POLLIN) != 0)
      {
        ssize_t bytes = read(output, out_buf, bufsize-1);
        ++rct;
        if (bytes > 0)
        {
//          std::cerr << "runCompiler: stdout bytes " << bytes << "\n";
          out_buf[bytes] = '\0'; 
          out_str.append(out_buf);
        } 
        else if ((bytes == -1) && (errno != EAGAIN)) P_MESSAGE2("Error reading standard output: ", strerror(errno));
      }
      if ((fds[iout].revents & POLLHUP) != 0)
      {
        outfin = true;
        if (!errfin) { fds[iout] = fds[ierr]; ierr = iout; }
        nfds -= 1;
        if (close (output) == -1) P_MESSAGE2("Error closing compiler output pipe: ", strerror(errno));
      }
    }
    if (!errfin)
    {
      P_ERROR_X2((fds[ierr].revents & (POLLNVAL|POLLERR)) == 0, "Error pipe error -> ", err_str);
      if ((fds[ierr].revents & POLLIN) != 0)
      {
        ssize_t bytes = read(error, out_buf, bufsize-1);
        ++ect;
        if (bytes > 0)
        {
//          std::cerr << "runCompiler: stderr bytes " << bytes << "\n";
          out_buf[bytes] = '\0'; 
          err_str.append(out_buf);
        } 
        else if ((bytes == -1) && (errno != EAGAIN)) P_MESSAGE2("Error reading standard error: ", strerror(errno));
      }
      if ((fds[ierr].revents & POLLHUP) != 0)
      {
        errfin = true;
        nfds -= 1;
        if (close (error) == -1) P_MESSAGE2("Error closing compiler error pipe: ", strerror(errno));
      }
    }
  } while (!outfin || !infin || !errfin);
  delete[] out_buf; out_buf = 0;
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

KNSystem::KNSystem(const std::string& sysName, const std::string& sysType, size_t usederi)
#ifdef GINAC_FOUND
 : useVectorField(false)
#endif
{
  std::string sname(sysName);
  if (sysName.find('/') == std::string::npos) sname = "./" + sysName;
#ifdef WIN32
  // replace the unix '/' with windows '\'
  std::string::size_type index = 0;
  while ((index = sname.find('/', index)) != std::string::npos) sname.at(index) = '\\';
#endif

  // Finding the executable directory
  std::string executableFile;
#ifdef __APPLE__
  CFURLRef bundleURL = CFBundleCopyExecutableURL(CFBundleGetMainBundle());
  if(bundleURL)
  {
    CFStringRef cfPath = CFURLCopyFileSystemPath(bundleURL, kCFURLPOSIXPathStyle);
    if(cfPath) 
    {
      CFIndex ssize = CFStringGetLength(cfPath)+1;
      char *buf = new char[ssize];
      CFStringGetCString(cfPath, buf, ssize, kCFStringEncodingMacRoman);
      executableFile = buf;
      delete[] buf; buf = 0;
      CFRelease(cfPath);
    }
    CFRelease(bundleURL);
  }
#elif __linux__
  std::string workingDir;
  std::ostringstream procf;
  procf << "/proc/" << getpid() << "/exe";
  char *buf = new char[512];
  const ssize_t bsfn = readlink(procf.str().c_str(), buf, 511);
  if ( bsfn != -1) { buf[bsfn] = '\0'; executableFile = buf; }
  delete[] buf; buf = 0;
#elif WIN32
  char *buf = new char[MAX_PATH]; //always use MAX_PATH for filepaths
  GetModuleFileName(NULL, buf, MAX_PATH*sizeof(char));
  executableFile = buf;
  delete[] buf; buf = 0;
#endif
  std::string executableDir(executableFile);
  std::string::size_type sidx = executableDir.find_last_of(DIRSEP);
  if (sidx != std::string::npos) executableDir.erase(sidx,std::string::npos);
  
  // making the shared object or loading the vector field file
  std::string objname;
#ifdef GINAC_FOUND
  useVectorField = makeSystem(objname, sname, sysType, executableDir);
  if (useVectorField)
  {
    try {
      makeSymbolic(sname);
    }
    catch (KNException& ex)
    {
      useVectorField = false;
      throw (ex);
    }
  }
#else
  makeSystem(objname, sname, sysType, executableDir);
#endif
  // filling up the function pointers
#ifdef GINAC_FOUND
  if (!useVectorField)
  {
#endif
    char * c_objname = new char[objname.size()+1];
    strncpy (c_objname, objname.c_str(), objname.size()+1);
    handle = tdlopen(c_objname);
    delete[] c_objname; c_objname = 0;
    P_ERROR_X5(handle != 0, "Cannot open system definition file. Error code", tdlerror(), ". The offending file was '", objname, "'.");
    
    tdlerror();    /* Clear any existing error */
    void* nd_res = tdlsym(handle, "sys_ndim");
    v_ndim = reinterpret_cast<tp_sys_ndim>( nd_res );
    P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_ndim(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_npar = reinterpret_cast<tp_sys_npar>(tdlsym(handle, "sys_npar"));
    P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_npar(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_ntau = reinterpret_cast<tp_sys_ntau>(tdlsym(handle, "sys_ntau"));
    P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_ntau(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_nderi = reinterpret_cast<tp_sys_nderi>(tdlsym(handle, "sys_nderi"));
    P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_nderi(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_nevent = reinterpret_cast<tp_sys_nevent>(tdlsym(handle, "sys_nevent"));
    if ((error = tdlerror()) != 0) v_nevent = 0;
    
    tdlerror();    /* Clear any existing error */
    v_tau = reinterpret_cast<tp_sys_tau>(tdlsym(handle, "sys_tau"));
    if ((error = tdlerror()) != 0) v_tau = 0;
    
    tdlerror();    /* Clear any existing error */
    v_dtau = reinterpret_cast<tp_sys_dtau>(tdlsym(handle, "sys_dtau"));
    if ((error = tdlerror()) != 0) v_dtau = 0;
    
    tdlerror();    /* Clear any existing error */
    v_mass = reinterpret_cast<tp_sys_mass>(tdlsym(handle, "sys_mass"));
    if ((error = tdlerror()) != 0) v_mass = 0;
    
    tdlerror();    /* Clear any existing error */
    v_rhs = reinterpret_cast<tp_sys_rhs>(tdlsym(handle, "sys_rhs"));
    if ((error = tdlerror()) != 0) v_rhs = 0;
    
    tdlerror();    /* Clear any existing error */
    v_deri = reinterpret_cast<tp_sys_deri>(tdlsym(handle, "sys_deri"));
    if ((error = tdlerror()) != 0) v_deri = 0;
    
    tdlerror();    /* Clear any existing error */
    v_stpar = reinterpret_cast<tp_sys_stpar>(tdlsym(handle, "sys_stpar"));
    P_ERROR_X3((error = tdlerror()) == 0, "Cannot find sys_stpar(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_stsol = reinterpret_cast<tp_sys_stsol>(tdlsym(handle, "sys_stsol"));
    if ((error = tdlerror()) != 0) v_stsol = 0;
    
    tdlerror();    /* Clear any existing error */
    v_parnames = reinterpret_cast<tp_sys_parnames>(tdlsym(handle, "sys_parnames"));
    if ((error = tdlerror()) != 0) v_parnames = 0;
    
    /* Vectorized versions */
    tdlerror();    /* Clear any existing error */
    v_p_tau = reinterpret_cast<tp_sys_p_tau>(tdlsym(handle, "sys_p_tau"));
    if ((error = tdlerror()) != 0) v_p_tau = 0;
    if ((v_tau == 0) && (v_p_tau == 0)) P_MESSAGE3("Cannot find either sys_tau() or sys_p_tau(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_p_dtau = reinterpret_cast<tp_sys_p_dtau>(tdlsym(handle, "sys_p_dtau"));
    if ((error = tdlerror()) != 0) v_p_dtau = 0;
    if ((v_dtau == 0) && (v_p_dtau == 0)) P_MESSAGE3("Cannot find either sys_dtau() or sys_p_dtau(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_p_rhs = reinterpret_cast<tp_sys_p_rhs>(tdlsym(handle, "sys_p_rhs"));
    if ((error = tdlerror()) != 0) v_p_rhs = 0;
    if ((v_rhs == 0) && (v_p_rhs == 0)) P_MESSAGE3("Cannot find either sys_rhs() or sys_p_rhs(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_p_deri = reinterpret_cast<tp_sys_p_deri>(tdlsym(handle, "sys_p_deri"));
    if ((error = tdlerror()) != 0) v_p_deri = 0;
    if ((v_deri == 0) && (v_p_deri == 0)) P_MESSAGE3("Cannot find either sys_deri() or sys_p_deri(): ", error, ".");
    
    tdlerror();    /* Clear any existing error */
    v_p_stsol = reinterpret_cast<tp_sys_p_stsol>(tdlsym(handle, "sys_p_stsol"));
    if ((error = tdlerror()) == 0) v_p_stsol = 0;
    if ((v_stsol == 0) && (v_p_stsol == 0)) P_MESSAGE3("Cannot find either sys_stsol() or sys_p_stsol(): ", error, ".");
    nderi = std::min ( (*v_nderi)(), usederi );  
#ifdef GINAC_FOUND
  } else
  {
    handle = 0;
    v_ndim = 0;
    v_npar = 0;
    v_ntau = 0;
    v_nderi = 0;
    v_nevent = 0;
    v_tau = 0; 
    v_dtau = 0; 
    v_mass = 0; 
    v_rhs = 0;
    v_deri = 0;
    v_p_tau = 0; 
    v_p_dtau = 0; 
    v_p_rhs = 0; 
    v_p_deri = 0; 
    v_p_event = 0;
    v_stpar = 0; 
    v_stsol = 0; 
    v_parnames = 0;
    nderi = 2;
  }
#endif
  f.init(this->ndim()), f_eps.init(this->ndim());
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
  for (size_t i = 0; i < M.row(); i++)
  {
    for (size_t j = 0; j < M.col(); j++)
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

void KNSystem::generateSystem(const std::string& vffile, const std::string& shobj, const std::string& executableDir)
{
#ifdef GINAC_FOUND
  // parsing the vector field file
  std::ostringstream cxxcode;
  std::map<std::string, std::string> options;
  VectorField vf;
  vf.ReadXML(vffile);
  int pserr = vf.ProcessSymbols();
  if (pserr == -1) P_MESSAGE1("Could not parse the vector field definition.");
  if (vf.isStateDependentDelay()) P_MESSAGE1("State dependent delays are not suported yet.");
  vf.PrintKnut(cxxcode, options);
  try {
//  	std::cout << "Trying compiler ?? is ?? working??\n";
    runCompiler(cxxcode.str(), shobj, executableDir);
  }
  catch (KNException& ex)
  {
//  	std::cout << "Compiler is not working! " << ex.str() << "\n";
  	std::cout.flush();
    workingCompiler = false;
  }
#endif
}

static inline bool to_compile(const struct stat *sbuf_so, const struct stat *sbuf_src)
{
#ifdef __APPLE__
  if (sbuf_so->st_mtimespec.tv_sec < sbuf_src->st_mtimespec.tv_sec) return true;
  else if (sbuf_so->st_mtimespec.tv_sec == sbuf_src->st_mtimespec.tv_sec) return (sbuf_so->st_mtimespec.tv_nsec <= sbuf_src->st_mtimespec.tv_nsec);
  else return false;
#else
  return (sbuf_so->st_mtime <= sbuf_src->st_mtime);
#endif
}

// Static member. Compile if necessary
// input : sysName : system definition file
//         sysType : type of system definition, .vf, .cpp, .so
//         executableDir : the directory of the software
// output : soname : the object file to be loaded
// return value : true if vector field file is to be used directly
//                false if object file is to be loaded
bool KNSystem::makeSystem(std::string& soname, const std::string& sysName, const std::string& sysType, const std::string& executableDir)
{
  std::string type(sysType);
  std::string cxxcode;
  std::string cxxfile;
  std::string sofile;
  std::string vffile;
  struct stat sbuf_so;
  struct stat sbuf_cxx;
  struct stat sbuf_vf;
  bool use_vf = false;
  bool compile_vf = false, compile_cxx = false;

  std::transform(type.begin(), type.end(), type.begin(), ::tolower);
  if (type == "")
  {
    std::string ext = sysName.substr(sysName.find_last_of('.')+1,std::string::npos);
    if (ext == "so") type = SYS_TYPE_SO;
    if (ext == "cpp") type = SYS_TYPE_CXX;
    if (ext == "cc") type = SYS_TYPE_CXX;
    if (ext == "C") type = SYS_TYPE_CXX;
    if (ext == "vf") type = SYS_TYPE_VFC;
  }
  if ((type == SYS_TYPE_VFC)&&(!workingCompiler)) type = SYS_TYPE_VF0;
  if (type == SYS_TYPE_SO)
  {
    int res_so = stat(sysName.c_str(), &sbuf_so);
    if (res_so != 0)
    P_ERROR_X3 (res_so == 0, "Shared Object system definition '", sysName, "' is missing.");
    soname = sysName;
    return use_vf;
  } else if (type == SYS_TYPE_CXX)
  {
    cxxfile = sysName;
    sofile = sysName + ".so";
    int res_so = stat(sofile.c_str(), &sbuf_so);
    int res_cxx = stat(cxxfile.c_str(), &sbuf_cxx);
    P_ERROR_X3 (res_cxx == 0, "C++ system definition '", sysName, "' is missing.");
    if (res_so != 0) compile_cxx = true;
    else compile_cxx = to_compile(&sbuf_so, &sbuf_cxx);
    // try to compile even if compiler is know not to work
    if (compile_cxx) compileSystem(cxxfile, sofile, executableDir);
    soname = sofile;
    return use_vf;
  } else if (type == SYS_TYPE_VFC)
  {
    vffile = sysName;
    sofile = sysName + ".so";
    int res_so = stat(sofile.c_str(), &sbuf_so);
    int res_vf = stat(vffile.c_str(), &sbuf_vf);
    P_ERROR_X3 (res_vf == 0, "Vector field system definition '", sysName, "' is missing.");
    if (res_so == 0) compile_vf = to_compile(&sbuf_so, &sbuf_vf);
    else compile_vf = true;
    if (compile_vf)
    {
      generateSystem(vffile, sofile, executableDir);
      soname = sofile;
      if (workingCompiler) return use_vf;
      type = SYS_TYPE_VF0;
    } else
    {
      soname = sofile;
      return use_vf;
    }
  } 
  if (type == SYS_TYPE_VF0)
  {
    vffile = sysName;
    int res_vf = stat(vffile.c_str(), &sbuf_cxx);
    P_ERROR_X3 (res_vf == 0, "Vector field definition '", sysName, "' is missing.");
    soname = "";
    use_vf = true;
    return use_vf;
  }
  P_MESSAGE4("Invalid system definition type '", type, "' of ", sysName);
//   delete sbuf_vf; sbuf_vf = 0;
//   delete sbuf_cxx; sbuf_cxx = 0; 
//   delete sbuf_so; sbuf_so = 0;
  return use_vf;
}

void KNSystem::p_discrderi( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& p_xx, const KNVector& par, size_t sel, 
                          size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& p_vv )
{
  const double abs_eps_x1 = 1e-6;
  const double rel_eps_x1 = 1e-6;
  const double abs_eps_p1 = 1e-6;
  const double rel_eps_p1 = 1e-6;
  const double abs_eps_x2 = 1e-6;
  const double rel_eps_x2 = 1e-6;

  const size_t n = ndim();
  const size_t m = 2 * ntau() + 1;

  p_resize(time.size());
  // derivatives w.r.t. the dependent variables: x(t), x(t-tau1), etc.
  if ((nx == 1) && (np == 0))
  {
//    std::cout<<"@ " << vx[0] << "\n";
    p_xx_eps = p_xx;    
    p_rhs(p_fx, time, p_xx, par, sel);
    for (size_t j = 0; j < n; j++)
    {
      for (size_t idx=0; idx < time.size(); ++idx)
      {
        const double eps = abs_eps_x1 + rel_eps_x1 * fabs(p_xx(j, vx[0],idx));
        p_xx_eps(j, vx[0],idx) = p_xx(j, vx[0],idx) + eps;
      }
      p_rhs(p_fx_eps, time, p_xx_eps, par, sel);
      for (size_t p = 0; p < n; p++)
      {
        for (size_t idx=0; idx < time.size(); ++idx)
          out(p, j, idx) = (p_fx_eps(p, idx) - p_fx(p, idx)) / (p_xx_eps(j, vx[0],idx) - p_xx(j, vx[0],idx));
      }
      // redo the alteration
      for (size_t idx=0; idx < time.size(); ++idx) p_xx_eps(j, vx[0],idx) = p_xx(j, vx[0],idx);
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
    for (size_t p = 0; p < n; p++)
    {
      for (size_t idx=0; idx < time.size(); ++idx) out(p, 0, idx) = (p_fx_eps(p, idx) - p_fx(p, idx)) / eps;
    }
  }
  // second derivatives w.r.t. x
  if ((nx == 2) && (np == 0))
  {
//     std::cout<<"<?";
    p2_xx_eps = p_xx;
    for (size_t j = 0; j < n; j++)
    {
      p_deri(p2_dfx, time, p_xx, par, sel, 1, &vx[0], 0, vp, p_vv);
      // perturbing solution
      for (size_t idx=0; idx < time.size(); ++idx) p2_xx_eps(j, vx[1], idx) += abs_eps_x2 + rel_eps_x2 * fabs(p_xx(j, vx[1], idx));
      // perturbed rhs deri
      p_deri(p2_dfx_eps, time, p2_xx_eps, par, sel, 1, &vx[0], 0, vp, p_vv);
      for (size_t p = 0; p < n; p++)
      {
        for (size_t idx=0; idx < time.size(); ++idx) out(p, j, idx) = 0.0;
        for (size_t q = 0; q < n; q++)
        {
          for (size_t idx=0; idx < time.size(); ++idx) out(p, j, idx) += (p2_dfx_eps(p, q, idx) - p2_dfx(p, q, idx)) * p_vv(q, vx[0], idx);
        }
        for (size_t idx=0; idx < time.size(); ++idx) out(p, j, idx) /= (p2_xx_eps(j, vx[1], idx) - p_xx(j, vx[1], idx));
      }
      // restoring eps
      for (size_t idx=0; idx < time.size(); ++idx) p2_xx_eps(j, vx[1], idx) = p_xx(j, vx[1], idx);
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
    for (size_t p = 0; p < n; p++)
    {
      for (size_t q = 0; q < n; q++)
      {
        for (size_t idx=0; idx < time.size(); ++idx) out(p, q, idx) = (p2_dfx_eps(p, q, idx) - p2_dfx(p, q, idx)) / eps;
      }
    }
//     std::cout<<">";
  }
}

#ifdef GINAC_FOUND
void KNSystem::makeSymbolic(const std::string& vffile)
{
  // load the vector field file and generate expressions
//   std::cout << "This is from a vector field file.\n";
  VectorField vf;
  vf.ReadXML(vffile.c_str());
  int pserr = vf.ProcessSymbols();
  if (pserr == -1) P_MESSAGE1("Could not parse the vector field definition.");
  if (vf.isStateDependentDelay()) P_MESSAGE1("State dependent delays are not suported yet.");
  ex_ndim = vf.Knut_ndim();
  ex_npar = vf.Knut_npar();
  ex_ntau = vf.Knut_ntau();
  ex_nevent = vf.Knut_nevent();
  vf.Knut_tau(ex_tau);
  vf.Knut_tau_p(ex_tau_p);
  vf.Knut_RHS(ex_rhs);
  vf.Knut_RHS_p(ex_rhs_p);
  vf.Knut_RHS_x(ex_rhs_x);
  vf.Knut_RHS_xp(ex_rhs_xp);
  vf.Knut_RHS_xx(ex_rhs_xx);
  vf.Knut_stpar(ex_stpar);
  vf.Knut_stsol(ex_stsol);
  vf.Knut_parnames(ex_parnames);
}
#endif

//******** INTERFACE WITH THE SYSTEM DEFINITION ********//
#define ERR_THRESHOLD 8e6

size_t    KNSystem::ndim() const 
{
#ifdef GINAC_FOUND
  if (useVectorField) return ex_ndim;
#endif
//  if (ex_ndim != (*v_ndim)()) std::cout << "ndim: " << ex_ndim << "," << (*v_ndim)()<< "\n";
  return (*v_ndim)();
}

size_t    KNSystem::npar() const 
{
#ifdef GINAC_FOUND
  if (useVectorField) return ex_npar;
#endif
//  if (ex_npar != (*v_npar)()) std::cout << "npar: " << ex_npar << "," << (*v_npar)() << "\n";
  return (*v_npar)();
}

size_t    KNSystem::ntau() const
{
#ifdef GINAC_FOUND
  if (useVectorField) return ex_ntau;
#endif
//  if (ex_ntau != (*v_ntau)()) std::cout << "ntau: " << ex_ntau << "," << (*v_ntau)()<< "\n";
  return (*v_ntau)();
}

// Vectorized versions
void   KNSystem::p_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par )
{
#ifdef GINAC_FOUND
  if (useVectorField) sym_tau(out, time, par);
  else
  {
#endif
    if (v_p_tau != 0) (*v_p_tau)(out, time, par);
    else
    {
      for (size_t i=0; i<time.size(); ++i)
      {
        KNVector tout(out, i);
        (*v_tau)(tout, time(i), par);
      }
    }
//    KNArray2D<double> out_p(out);
//    sym_tau(out_p, time, par);
//    double max = 0;
//    size_t ip=0, jp=0;
//    for (size_t i=0; i < out.dim1(); ++i)
//    for (size_t j=0; j < out.dim2(); ++j) if (fabs(out_p(i,j)-out(i,j)) > max) { max = (fabs(out_p(i,j)-out(i,j))); ip=i; jp=j; }
//    if ((max) > ERR_THRESHOLD*std::numeric_limits<double>::epsilon( ) ) std::cout << "tau: " << max << " at (" << ip << "," << jp << ") sym="<< out_p(ip,jp) << " comp=" << out(ip,jp) << "\n";
#ifdef GINAC_FOUND
  }
#endif
}

void   KNSystem::p_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp )
{
#ifdef GINAC_FOUND
  if (useVectorField) sym_dtau(out, time, par, vp);
  else
  {
#endif
    if (v_p_dtau != 0) (*v_p_dtau)(out, time, par, vp);
    else
    {
      for (size_t i=0; i<time.size(); ++i)
      {
        KNVector tout(out, i);
        (*v_dtau)(tout, time(i), par, vp);
      }
    }
//    KNArray2D<double> out_p(out);
//    sym_dtau(out_p, time, par, vp);
//    double max = 0;
//    for (size_t i=0; i < out.dim1(); ++i)
//    for (size_t j=0; j < out.dim2(); ++j) if (fabs(out_p(i,j)-out(i,j)) > max) max = fabs((out_p(i,j)-out(i,j)));
//    if ((max) > ERR_THRESHOLD*std::numeric_limits<double>::epsilon( ) ) std::cout << "dtau: " << max << "\n";
#ifdef GINAC_FOUND
  }
#endif
}

void   KNSystem::mass(KNArray1D<double>& out)
{
  if (v_mass != 0)(*v_mass)(out);
  else
  {
    for(size_t i=0; i<out.size(); ++i) out(i) = 1.0;
  }
}

void   KNSystem::p_rhs( KNArray2D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel )
{
#ifdef GINAC_FOUND
  if (useVectorField) sym_rhs(out, time, x, par, sel);
  else
  {
#endif
    if (v_p_rhs != 0) (*v_p_rhs)(out, time, x, par, sel);
    else
    {
      for (size_t i=0; i<time.size(); ++i)
      {
        KNVector vout(out, i);
        KNMatrix xxin(x, i);
        (*v_rhs)(vout, time(i), xxin, par);
      }
    }
//    KNArray2D<double> out_p(out);
//    sym_rhs(out_p, time, x, par, sel);
//    double max = 0;
//    for (size_t i=0; i < out.dim1(); ++i)
//    for (size_t j=0; j < out.dim2(); ++j) if (fabs(out_p(i,j)-out(i,j)) > max) max = (fabs(out_p(i,j)-out(i,j)));
//    if ((max) > ERR_THRESHOLD*std::numeric_limits<double>::epsilon( ) ) std::cout << "p_rhs: " << max << "\n";
//    std::cout << "+";
#ifdef GINAC_FOUND
  }
#endif
}

void   KNSystem::p_deri( KNArray3D<double>& out, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par,
                size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv )
{
#ifdef GINAC_FOUND
  if (useVectorField) sym_deri(out, time, x, par, sel, nx, vx, np, vp, vv);
  else
  {
#endif
    if ((nderi >= 2) || ((nderi == 1) && ((nx == 1 && np == 0) || (nx == 0 && np == 1))))
    {
      if (v_p_deri != 0)
      {
        (*v_p_deri)(out, time, x, par, sel, nx, vx, np, vp, vv);
      } else
      {
        for (size_t i=0; i<time.size(); ++i)
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
//    KNArray3D<double> out_p(out);
//    sym_deri(out_p, time, x, par, sel, nx, vx, np, vp, vv);
//    double max = 0;
//    size_t ip=0, jp=0, kp=0;
//    for (size_t i=0; i < out.dim1(); ++i)
//    for (size_t j=0; j < out.dim2(); ++j) 
//    for (size_t k=0; k < out.dim3(); ++k) if (fabs(out_p(i,j,k)-out(i,j,k)) > max) { max = (fabs(out_p(i,j,k)-out(i,j,k))); ip=i; jp=j; kp=k; }
//    if ((max) > pow(ERR_THRESHOLD,nx+np)*std::numeric_limits<double>::epsilon( ) ) 
//    { 
//      std::cout << "p_deri("<< nx << "," << np << ") " << max << " at (" << ip << "," << jp << "," << kp << ") ";
//      if (np) std::cout << "vp=[";
//      for (size_t i=0; i<np; i++) std::cout << vp[i] << ",";
//      if (np) std::cout << "]";
//      if (nx) std::cout << "vx=[";
//      for (size_t i=0; i<nx; i++) std::cout << vx[i] << ",";
//      if (nx) std::cout << "] ";
//      std::cout << "sym="<< out_p(ip,jp,kp) << " comp=" << out(ip,jp,kp);
//      if (nx) std::cout << " f=" << ex_rhs_x[ip + (vx[0] + jp*ex_ntau)*ex_ndim];
//      std::cout << "\n";
//    }
//    std::cout << "+";
#ifdef GINAC_FOUND
  }
#endif
}

// Setting the starting point
void   KNSystem::stpar(KNVector& par) const
{
#ifdef GINAC_FOUND
  if (useVectorField) sym_stpar(par);
  else
  {
#endif
    (*v_stpar)(par);
//    KNVector par_p(par);
//    sym_stpar(par_p);
//    double max = 0;
//    for (size_t i=0; i < par_p.size(); ++i) if (fabs(par_p(i)-par(i)) > max) max = (fabs(par_p(i)-par(i)));
//    if ((max) > ERR_THRESHOLD*std::numeric_limits<double>::epsilon( ) ) std::cout << "stpar: " << max << "\n";
#ifdef GINAC_FOUND
  }
#endif
}

void   KNSystem::stsol(KNArray2D<double>& out, const KNArray1D<double>& time) const
{
#ifdef GINAC_FOUND
  if (useVectorField) sym_stsol(out, time);
  else
  {
#endif
    if (v_p_stsol != 0) (*v_p_stsol)(out, time);
    else if (v_stsol != 0)
    {
      for (size_t i=0; i<time.size(); ++i)
      {
        KNVector vf(out, i);
        (*v_stsol)(vf, time(i));
      }
    }
//    KNArray2D<double>& out_p(out);
//    sym_stsol(out_p, time);
//    double max = 0;
//    for (size_t i=0; i < out.dim1(); ++i)
//    for (size_t j=0; j < out.dim2(); ++j) if (fabs(out_p(i,j)-out(i,j)) > max) max = (fabs(out_p(i,j)-out(i,j)));
//    if ((max) > ERR_THRESHOLD*std::numeric_limits<double>::epsilon( ) ) std::cout << "stsol: " << max << "\n";
#ifdef GINAC_FOUND
  }
#endif
}

void   KNSystem::parnames(const char *names[]) const
{
#ifdef GINAC_FOUND
  if (useVectorField) sym_parnames(names);
  else
  {
#endif
    if (v_parnames != 0)
    {
      (*v_parnames)(names);
    }
#ifdef GINAC_FOUND
  }
#endif
}

#ifdef GINAC_FOUND

size_t KNSystem::sym_ndim() 
{
  return ex_ndim;
}

size_t KNSystem::sym_npar() 
{
  return ex_npar;
}

size_t KNSystem::sym_ntau() 
{
  return ex_ntau;
}

void KNSystem::sym_tau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par ) 
{
  VectorField::Knut_tau_eval( ex_tau, out, time, par, ex_ndim, ex_npar, ex_ntau );
}

void KNSystem::sym_dtau( KNArray2D<double>& out, const KNArray1D<double>& time, const KNVector& par, size_t vp )
{
  VectorField::Knut_tau_p_eval( ex_tau_p, out, time, par, vp, ex_ndim, ex_npar, ex_ntau );
}

void KNSystem::sym_rhs( KNArray2D<double>& out_p, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par, size_t sel )
{
#ifdef GINAC_FOUND
  VectorField::Knut_RHS_eval(ex_rhs, out_p, time, x, par, sel );
#endif
}

void KNSystem::sym_deri( KNArray3D<double>& out_p, const KNArray1D<double>& time, const KNArray3D<double>& x, const KNVector& par,
                size_t sel, size_t nx, const size_t* vx, size_t np, const size_t* vp, const KNArray3D<double>& vv )
{
#ifdef GINAC_FOUND
  if ((np == 1)&&(nx == 0))
  {
    VectorField::Knut_RHS_p_eval( ex_rhs_p, out_p, time, x, par, sel, vp[0], ex_ndim, ex_npar, ex_ntau );
  }
  if ((np == 0)&&(nx == 1))
  {
    VectorField::Knut_RHS_x_eval( ex_rhs_x, out_p, time, x, par, sel, vx[0], ex_ndim, ex_npar, ex_ntau );
  }
  if ((np == 1)&&(nx == 1))
  {
    VectorField::Knut_RHS_xp_eval( ex_rhs_xp, out_p, time, x, par, sel, vx[0], vp[0], ex_ndim, ex_npar, ex_ntau );
  }
  if ((np == 0)&&(nx == 2))
  {
    VectorField::Knut_RHS_xx_eval( ex_rhs_xx, out_p, time, x, vv, par, sel, vx[0], vx[1], ex_ndim, ex_npar, ex_ntau );
  }
#endif                
}

void KNSystem::sym_stpar(KNVector& par) const
{
  for (size_t i = 0; i < ex_npar; ++i)
  {
    par(i) = ex_stpar[i];
  }
}

void KNSystem::sym_stsol(KNArray2D<double>& out_p, const KNArray1D<double>& time) const
{
  VectorField::Knut_RHS_stsol_eval( ex_stsol, out_p, time );
}

void KNSystem::sym_parnames(const char *names[]) const
{
  for (size_t i = 0; i < ex_npar; ++i)
  {
    names[i] = ex_parnames[i].c_str();
  }
}
#endif
