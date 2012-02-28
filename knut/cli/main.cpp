// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "basecomp.h"
#include "knerror.h"
#include <unistd.h>
// #include "constants.h" included in basecomp
// #include <iostream> included in basecomp

int main(int argc, const char** argv)
{
  // parameters
  KNConstants*  params = 0;
  const char*  branchFile = 0;

  bool save = false;
  try
  {
    // argument parsing
    for (int acnt = 1; acnt < argc;  acnt++)
    {
      if (argv[acnt][0] == '-')
      {
        switch (argv[acnt][1])
        {
          case 's':
            save = true;
          case 'c':
            params = new KNConstants;
            {
              std::string constFile(argv[++acnt]);
              std::string cfdir(constFile);
              char* cwd_ptr = new char[512];
              getcwd(cwd_ptr, 511);
              P_ERROR_X1(cwd_ptr != 0, "Cannot obtain CWD.");
              std::string  cwd(cwd_ptr);
              delete[] cwd_ptr;
              cwd += '/';
              
              size_t found = cfdir.rfind('/');
              if (found != std::string::npos) cfdir.erase(found);
              else cfdir = cwd;
              int err = chdir (cfdir.c_str());
//              std::cout << "Changed directory " << cfdir.c_str() << "\n";
//              std::cout << "Previous directory " << cwd.c_str() << "\n";
              P_ERROR_X1(err == 0, "Error changing directory.");
              if (constFile[0] != '/') constFile.insert(0, cwd);
              params->loadXmlFile(constFile);
            }
            if (save) { params->printXmlFile(std::cout); exit(0); }
            break;
          case 'b':
            branchFile = argv[++acnt];
            break;
          case 'v':
#ifdef HAVE_CONFIG_H
            std::cout << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << " (" << PACKAGE_REVISION << ")\n";
#endif
            exit(0);
            break;
          default:
            P_MESSAGE1("Unexpected command line argument.");
            break;
        }
      }
      else
      {
        P_MESSAGE1("Unexpected command line argument.");
      }
    }
  
    if (params == 0)
    {
      P_MESSAGE1("Missing constants file.");
      exit(-1);
    }
  
    KNCliContinuation comp(*params);
    if (branchFile == 0) comp.run();
    else comp.run(branchFile);
  }
  catch (KNException ex)
  {
    KNCliContinuation::printException(ex);
#ifdef DEBUG
    throw(-1);
#endif
    exit(-1);
  }

  delete params;
  return 0;
}
