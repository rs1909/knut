// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the root package's directory
//
// ------------------------------------------------------------------------- //

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "basecomp.h"
#include "pderror.h"
// #include "constants.h" included in basecomp
// #include <iostream> included in basecomp

/*

sys-milltwofull.so   SYSNAME
0       LABEL
0 1 0      TYPE, CP, NPARX, PARX ....
200 4 5     NINT, NDEG, NMUL, STAB, NMAT
12 12 4 4    NINT1, NINT2, NDEG1, NDEG2 (for torus computations only)
126 -100.0 100.0  STEPS, CPMIN, CPMAX
0.1 0.01 0.2 0.1  DS, DSMIN, DSMAX, DSSTART
1e-4 1e-4 1e-4      EPSC, EPSR, EPSS
30 30 30     NITC, NITR, NITS

*/

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
  const char*  outFile = 0;
  const char*  inFile = 0;
  const char*  branchFile = 0;

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
          params->loadFile(constFile);
          params->printFile(std::cout);
          break;
        case 'i':
          inFile = argv[++acnt];
          break;
        case 'o':
          outFile = argv[++acnt];
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
          P_MESSAGE("Unexpected command line argument.\n");
          break;
      }
    }
    else
    {
      P_MESSAGE("Unexpected command line argument.\n");
    }
  }

  if (params == 0)
  {
    std::cout << "Error: missing constants file.\n";
    exit(1);
  }
  if ((inFile == 0) && (params->getLabel() != 0))
  {
    std::cout << "Error: missing input file.\n";
    exit(1);
  }
  if (outFile == 0)
  {
    outFile = "pdde.mat";
    std::cout << "Warning: missing output file, using \"" << outFile << "\" instead.\n";
  }
  if (branchFile == 0)
  {
    branchFile = "pdde.br";
    std::cout << "Warning: missing branch file, using \"" << branchFile << "\" instead.\n";
  }
  if (inFile == 0) params->setInputFileText("");
  else params->setInputFileText(inFile);
  params->setOutputFileText(outFile);
  
  CLIComp comp(*params);
  comp.run(branchFile);

  delete params;
}
