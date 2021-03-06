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

#include "pointtype.h"
#include "mat4data.h"

class KNDdePeriodicSolution : public KNAbstractPeriodicSolution
{
  public:

    // constructor
    KNDdePeriodicSolution(KNAbstractContinuation* cnt, KNExprSystem& sys, 
      KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_, size_t nint, size_t ndeg, size_t nmul = MULTIPLIERS);
    ~KNDdePeriodicSolution() override;

    void    Stability(bool init) override;

    // algorithms
    void    SwitchTFLP(BranchSW type, double ds) override;   // switches branch with testFunct
    void    SwitchTFPD(double ds) override;   // switches branch with testFunct
    void    SwitchTFHB(double ds) override;   // switches branch with testFunct
    void    SwitchTRSol(KNVector& Sol, const KNVector& mshint, const KNVector& mshdeg)
    {
      colloc->exportProfile(Sol, mshint, mshdeg, sol);
    } // starting data for tori: solution
    void    SwitchTFTRTan(KNVector& Re, KNVector& Im, double& alpha, const KNVector& mshint, const KNVector& mshdeg);   // starting data for tori: tangent WITH testFunct
 
    void    Plot(GnuPlot& pl);
    inline void    print(char* file)
    {
      std::ofstream ff(file);
      for (size_t i = 0; i < colloc->nInt()*colloc->nDeg() + 1; i++)
      {
        for (size_t j = 0; j < colloc->nDim(); j++) ff << sol(colloc->nDim()*i + j) << "\t";
        ff << "\n";
      }
    }
    void ReadNull(std::ifstream& file);
    void Read(std::ifstream& file);
    void Write(std::ofstream& file);

  private:
    // These are shadowing the construct/destruct of base classes
    void    construct() override;
    void    destruct() override;

    // internal member-functions
//   inline double SolNorm( KNVector& sol, KNVector& par );

/// MAKE IT VIRTUAL
    void jacobian(
      KNSparseBlockMatrix& AA, KNBlockVector& RHS,                      // output
      KNVector& parPrev, KNVector& par,                           // parameters
      KNVector& solPrev, KNVector& sol,                           // solution
      KNArray1D<size_t>& varMap,                                // contains the variables. If cont => contains the P0 too.
      double ds, bool cont                                    // ds stepsize, cont: true if continuation
    ) override;
    void postProcess() override;

    // multipliers
//    KNVector       mRe;
//    KNVector       mIm;
    // number of trivial multipliers
//    int   nTrivMulLP, nTrivMulPD, nTrivMulNS;

    // for the rotation phase conditions
//    KNArray1D<int> rotRe;
//    KNArray1D<int> rotIm;

    // stability matrix
    KNSparseMatrixPolynomial*  jacStab;

    // test functionals
    KNArray1D<KNAbstractTestFunctional*> testFunct;
    // for collocation
    KNDdeBvpCollocation*     colloc;
};

//-------------------------------------------------------//
// 5R   2 13R         2 13R         2 13R         1 2R 1 2R
//
// LABEL  NORM          CP            P0             U   IT
//     0  0.000000e+00  0.000000e+00  0.000000e+00   0    0
//-------------------------------------------------------//

inline void parNamePrint(std::ostream& out, size_t npar, Var cp, const KNArray1D<Var>& var, const std::vector<std::string>& parNames)
{
//   out.fill(' ');
  out.unsetf(std::ios::adjustfield);
//   out.setf(std::ios::left);

  out << std::left << std::setfill(' ') << std::setw(5) <<  "LABEL" << " "; // 2 at the end
  out << std::left << std::setfill(' ') << std::setw(2) <<  "TP" << ' '; // 3 at the end
  out << ' ' << std::left << std::setfill(' ') << std::setw(12) << "NORM" << "  "; // 2 at the end
  out << ' ' << std::left << std::setfill(' ') << std::setw(12) << parNames.at(VarToIndex(cp,npar)) << "  "; // 2 at the end
  for (size_t j = 1; j < var.size(); j++)
  {
    out << ' ' << std::left << std::setfill(' ') << std::setw(12) << parNames.at(VarToIndex(var(j),npar)) << "  "; // 2 at the end
  }
  out << std::left << std::setfill(' ') << std::setw(2) << "U" << " ";
  out << std::left << std::setfill(' ') << std::setw(2) << "IT";
}

inline void parValuePrint(std::ostream& out, const KNVector& par, Var cp, const KNArray1D<Var>& var, size_t lb, BifType tp, double norm, size_t ustab, size_t it)
{
  const size_t npar = par.size() - (VarEnd - VarInternal);
  out.fill(' ');
  out.unsetf(std::ios::adjustfield);
  out.precision(6);
  out.setf(std::ios::uppercase | std::ios::scientific | std::ios::right);

  out << std::right << std::setfill(' ') << std::setw(5) << lb+1 << " "; // 2 at the end
  if (tp == BifType::BifLP) out << std::left << std::setfill(' ') << std::setw(2) << "LP" << ' '; // 2 at the end
  else if (tp == BifType::BifPD) out << std::left << std::setfill(' ') << std::setw(2) << "PD" << ' '; // 2 at the end
  else if (tp == BifType::BifNS) out << std::left << std::setfill(' ') << std::setw(2) << "NS" << ' '; // 2 at the end
  else if (tp == BifType::BifUN) out << std::left << std::setfill(' ') << std::setw(2) << "UN" << ' '; // 2 at the end
  else if (tp == BifType::BifMax) out << std::left << std::setfill(' ') << std::setw(2) << "EP" << ' '; // 2 at the end
  else if (tp == BifType::BifEndPoint) out << std::left << std::setfill(' ') << std::setw(2) << "EP" << ' '; // 2 at the end
  else if (tp == BifType::BifNoConvergence) out << std::left << std::setfill(' ') << std::setw(2) << "MX" << ' '; // 2 at the end
  else out << std::left << std::setfill(' ') << std::setw(2) << "  " << ' '; // 2 at the end
  out << std::right << std::setfill(' ') << std::setw(13) << norm << "  "; // 2 at the end
  out << std::right << std::setfill(' ') << std::setw(13) << par(VarToIndex(cp,npar)) << "  "; // 2 at the end
  for (size_t j = 1; j < var.size(); j++)
  {
    out << std::right << std::setfill(' ') << std::setw(13) << par(VarToIndex(var(j),npar)) << "  "; // 2 at the end
  }
  out << std::left << std::setfill(' ') << std::setw(2) << ustab << " ";
  out << std::left << std::setfill(' ') << std::setw(2) << it;
}

#endif
