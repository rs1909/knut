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
// #include "constants.h" included in basecomp
// #include <iostream> included in basecomp

inline void parNamePrint(Vector& /*par*/, int npar, Array1D<Var>& var)
{
  for (int j = 1; j < var.Size(); j++) std::cout << "\t" << parType(npar, var(j) - VarPAR0) << parNum(npar, var(j) - VarPAR0) << "\t";
}

inline void parValuePrint(Vector& par, int /*npar*/, Array1D<Var>& var)
{
  for (int j = 1; j < var.Size(); j++) std::cout << "\t" << par(var(j) - VarPAR0);
}


int main(int argc, const char** argv)
{
  // parameters
  NConstants*  params = 0;
  const char*  constFile = 0;
  const char*  branchFile = 0;

  try
  {
    // argument parsing
    for (int acnt = 1; acnt < argc;  acnt++)
    {
      if (argv[acnt][0] == '-')
      {
        switch (argv[acnt][1])
        {
          case 'c':
            constFile = argv[++acnt];
            params = new NConstants;
            params->loadXmlFile(constFile);
//            params->printXmlFile(std::cout);
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
      exit(1);
    }
  
    CLIComp comp(*params);
    if (branchFile == 0) comp.run();
    else comp.run(branchFile);
  }
  catch (knutException ex)
  {
    CLIComp::printException(ex);
#ifdef DEBUG
    throw(-1);
#endif
    exit(-1);
  }

  delete params;
  return 0;
}
