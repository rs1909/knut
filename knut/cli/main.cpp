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
#include <new>
// #include "constants.h" included in basecomp
// #include <iostream> included in basecomp

int main(int argc, const char** argv)
{
  // parameters
  KNConstants* params = nullptr;
  KNExprSystem* sys = nullptr;

  int save {-1};
  int run {-1};
  try
  {
    // argument parsing
    for (int acnt = 1; acnt < argc;  acnt++)
    {
      if ( (argv[acnt][0] == '-') && (run != acnt - 1) )
      {
        switch (argv[acnt][1])
        {
          case 's':
            save = acnt;
            P_ERROR_X1(argv[acnt][2] == '\0', "Unexpected command line argument.");
            break;
          case 'c':
            run = acnt;
            P_ERROR_X1(argv[acnt][2] == '\0', "Unexpected command line argument.");
            break;
          case 'v':
            P_ERROR_X1(argv[acnt][2] == '\0', "Unexpected command line argument.");
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
        // the previous argument was "-c"
        // so this must be a file name
        if (run == acnt - 1)
        {
          params = new KNConstants;
          {
            std::string constFile{argv[acnt]};
            std::string cfdir{argv[acnt]};
            auto* cwd_ptr = new char[512];
            getcwd(cwd_ptr, 511);
            P_ERROR_X1(cwd_ptr != nullptr, "Cannot obtain CWD.");
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
            if (!constFile.empty ()) params->loadXmlFileV5(constFile);
            std::cout << "Need to compile " << params->getCompile () << "\n";
            if (params->getSystem ().empty ())
            {
              sys = new KNSystem (params->getSysName (), params->getCompile ());
//                 std::string tmp;
//                 sys -> toString (tmp);
//                 params->setSystemText (tmp);
//                 params->setSysName ("");
            } else
            {
        //      std::cout << params->getSystem () << "\n";
              sys = new KNExprSystem (params->getSystem (), params->getCompile ());
            }
          }
          if (save > 0) { params->saveXmlFileV5(""); exit(0); }
        }
      }
    } // end for

    if (params == nullptr)
    {
      P_MESSAGE1("Missing constants file.");
      exit(-1);
    }

    KNCliContinuation comp;
    KNDataFile *inputData = nullptr;
    if (params->getLabel() != 0)
    {
      inputData = new KNDataFile (params->getInputFile());
    }
    comp.run(sys, params, inputData);
    delete inputData;
    delete params;
    delete sys;
  }
  catch (KNException ex)
  {
    KNCliContinuation::printException(ex);
    delete params;
    delete sys;
#ifdef DEBUG
    throw(-1);
#endif
    exit(-1);
  }
  catch (std::bad_alloc& ba)
  {
    std::cerr << "bad_alloc caught: " << ba.what() << '\n';
  }
  return 0;
}

