// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef NPOINT_H
#define NPOINT_H

#include "basepoint.h"
#include "matrix.h"
#include "hypermatrix.h"
#include "testfunct.h"

#include "ncolloc.h"

#include "plot.h"

#include <fstream>
#include <iomanip>
#include <list>
#include <cmath>

#define MULTIPLIERS 8
#define BMATRICES 1

#include "pointtype.h"
#include "mat4data.h"

class Point : public PerSolPoint
{
  public:

    // constructor
    Point(System& sys, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg, int nmul = MULTIPLIERS, int nmat = BMATRICES);
    virtual ~Point();

    void    Stability();

    // algorithms
    void    SwitchTFLP(BranchSW type, double ds);   // switches branch with testFunct
    void    SwitchTFPD(double ds);   // switches branch with testFunct
    void    SwitchTFHB(double ds, std::ostream& out);   // switches branch with testFunct
    void    SwitchTRSol(Vector& Sol, const Vector& mshint, const Vector& mshdeg)
    {
      colloc->Export(Sol, mshint, mshdeg, sol);
    } // starting data for tori: solution
    void    SwitchTFTRTan(Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg);   // starting data for tori: tangent WITH testFunct
 
    void    Plot(GnuPlot& pl);
    inline void    Print(char* file)
    {
      std::ofstream ff(file);
      for (int i = 0; i < colloc->Nint()*colloc->Ndeg() + 1; i++)
      {
        for (int j = 0; j < colloc->Ndim(); j++) ff << sol(colloc->Ndim()*i + j) << "\t";
        ff << "\n";
      }
    }
    void ReadNull(std::ifstream& file);
    void Read(std::ifstream& file);
    void Write(std::ofstream& file);

  private:

    void    Construct();
    void    Destruct();

    // internal member-functions
//   inline double SolNorm( Vector& sol, Vector& par );

/// MAKE IT VIRTUAL
    void Jacobian
    (
      HyperMatrix& AA, HyperVector& RHS,                      // output
      Vector& parPrev, Vector& par,                           // parameters
      Vector& solPrev, Vector& sol,                           // solution
      Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
      double ds, bool cont                                    // ds stepsize, cont: true if continuation
    );


    // multipliers
//    Vector       mRe;
//    Vector       mIm;
    // number of trivial multipliers
//    int   nTrivMulLP, nTrivMulPD, nTrivMulNS;

    // for the rotation phase conditions
//    Array1D<int> rotRe;
//    Array1D<int> rotIm;

    // stability matrix
    StabMatrix   jacStab;

    // test functionals
    Array1D<baseTestFunct*> testFunct;
    // for collocation
    NColloc*     colloc;
};

//-------------------------------------------------------//
// 5R   2 13R         2 13R         2 13R         1 2R 1 2R
//
// LABEL  NORM          CP            P0             U   IT
//     0  0.000000e+00  0.000000e+00  0.000000e+00   0    0
//-------------------------------------------------------//

inline void parNamePrint(std::ostream& out, int npar, Var cp, const Array1D<Var>& var, const std::vector<std::string>& parNames)
{
//   out.fill(' ');
  out.unsetf(std::ios::adjustfield);
//   out.setf(std::ios::left);

  out << std::left << std::setfill(' ') << std::setw(5) <<  "LABEL" << " "; // 2 at the end
  out << std::left << std::setfill(' ') << std::setw(2) <<  "TP" << ' '; // 3 at the end
  out << ' ' << std::left << std::setfill(' ') << std::setw(12) << "NORM" << "  "; // 2 at the end
  out << ' ' << std::left << std::setfill(' ') << std::setw(12) << parNames.at(cp - VarPAR0) << "  "; // 2 at the end
  for (int j = 1; j < var.Size(); j++)
  {
    out << ' ' << std::left << std::setfill(' ') << std::setw(12) << parNames.at(var(j) - VarPAR0) << "  "; // 2 at the end
  }
  out << std::left << std::setfill(' ') << std::setw(2) << "U" << " ";
  out << std::left << std::setfill(' ') << std::setw(2) << "IT";
}

inline void parValuePrint(std::ostream& out, const Vector& par, Var cp, const Array1D<Var>& var, int lb, BifType tp, double norm, int ustab, int it)
{
  out.fill(' ');
  out.unsetf(std::ios::adjustfield);
  out.precision(6);
  out.setf(std::ios::uppercase | std::ios::scientific | std::ios::right);

  out << std::right << std::setfill(' ') << std::setw(5) << lb << " "; // 2 at the end
  if (tp == BifLP) out << std::left << std::setfill(' ') << std::setw(2) << "LP" << ' '; // 2 at the end
  else if (tp == BifPD) out << std::left << std::setfill(' ') << std::setw(2) << "PD" << ' '; // 2 at the end
  else if (tp == BifNS) out << std::left << std::setfill(' ') << std::setw(2) << "NS" << ' '; // 2 at the end
  else out << std::left << std::setfill(' ') << std::setw(2) << "  " << ' '; // 2 at the end
  out << std::right << std::setfill(' ') << std::setw(13) << norm << "  "; // 2 at the end
  out << std::right << std::setfill(' ') << std::setw(13) << par(cp - VarPAR0) << "  "; // 2 at the end
  for (int j = 1; j < var.Size(); j++)
  {
    out << std::right << std::setfill(' ') << std::setw(13) << par(var(j) - VarPAR0) << "  "; // 2 at the end
  }
  out << std::left << std::setfill(' ') << std::setw(2) << ustab << " ";
  out << std::left << std::setfill(' ') << std::setw(2) << it;
}

#endif
