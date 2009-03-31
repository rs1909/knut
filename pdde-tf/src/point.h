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
#include <list>
#include <cmath>

#define MULTIPLIERS 8
#define BMATRICES 1

#include "pointtype.h"
#include "mat4data.h"
#include "multipliers.h"


class Point : public BasePoint
{
  public:

    // constructor
    Point(System& sys, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg, int nmul = MULTIPLIERS, int nmat = BMATRICES);
    ~Point();

    void    Stability();

    // algorithms
    void    SwitchTFLP(BranchSW type, double ds);   // switches branch with testFunct
    void    SwitchTFPD(double ds);   // switches branch with testFunct
    void    SwitchTFHB(double ds);   // switches branch with testFunct
    void    SwitchTRSol(Vector& Sol, const Vector& mshint, const Vector& mshdeg)
    {
      colloc->Export(Sol, mshint, mshdeg, sol);
    } // starting data for tori: solution
    void    SwitchTFTRTan(Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg);   // starting data for tori: tangent WITH testFunct
    int     StartTF(Eqn FN, std::ostream& out = std::cout);
    inline void    setSym(int n, int* sRe, int* sIm)
    {
      rotRe.Init(n);
      rotIm.Init(n);
      for (int i = 0; i < n; i++)
      {
        rotRe(i) = sRe[i];
        rotIm(i) = sIm[i];
      }
    }
    inline void    setSym(Array1D<int>& sRe, Array1D<int>& sIm)
    {
      P_ASSERT(sRe.Size() == sIm.Size());
      rotRe.Init(sRe.Size());
      rotIm.Init(sRe.Size());
      for (int i = 0; i < sRe.Size(); i++)
      {
        rotRe(i) = sRe(i);
        rotIm(i) = sIm(i);
      }
    }


    inline double  NormMX()
    {
      double max = 0.0, min = 1.0e32;
      for (int i = 0; i < colloc->Nint()*colloc->Ndeg() + 1; i++)
      {
        double e = 0.0;
        for (int j = 0; j < colloc->Ndim(); j++) e += sqrt(sol(colloc->Ndim() * i + j) * sol(colloc->Ndim() * i + j));
        if (e > max) max = e;
        if (e < min) min = e;
      }
      return max -min;
    }
    int     UStab() { return unstableMultipliers(mRe, mIm, nTrivMulLP, nTrivMulPD, nTrivMulNS); }
    PtType  testBif() { return bifurcationType(mRe, mIm, nTrivMulLP, nTrivMulPD, nTrivMulNS); }
    void    clearStability() { mRe.Clear(); mIm.Clear(); }

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
    void BinaryRead(mat4Data& data, int n);
    void BinaryWrite(mat4Data& data, int n);

  private:

    void    Construct();
    void    Destruct();
    void    FillSol(System& sys_);

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
    Vector       mRe;
    Vector       mIm;
    // number of trivial multipliers
    int   nTrivMulLP, nTrivMulPD, nTrivMulNS;

    // for the rotation phase conditions
    Array1D<int> rotRe;
    Array1D<int> rotIm;

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

inline void parNamePrint(std::ostream& out, int npar, Var cp, const Array1D<Var>& var)
{
  out.fill(' ');
  out.unsetf(std::ios::adjustfield);
  out.setf(std::ios::left);

  out << "LABEL   NORM          ";
  out << ' ';
  out << parType(npar, cp - VarPAR0);
  out.width(11);
  out << parNum(npar, cp - VarPAR0) << "  ";
  for (int j = 1; j < var.Size(); j++)
  {
    out << ' ';
    out << parType(npar, var(j) - VarPAR0);
    out.width(11);
    out << parNum(npar, var(j) - VarPAR0) << "  ";
  }
  out << " U IT";
}

inline void parValuePrint(std::ostream& out, const Vector& par, Var cp, const Array1D<Var>& var, int lb, double norm, int ustab, int it)
{
  out.fill(' ');
  out.unsetf(std::ios::adjustfield);
  out.precision(6);
  out.setf(std::ios::uppercase | std::ios::scientific | std::ios::right);

  out.width(5);
  out.setf(std::ios::right);
  out << lb << "  ";
  out.width(13);
  out.setf(std::ios::right);
  out << norm << "  ";
  out.width(13);
  out.setf(std::ios::right);
  out << par(cp - VarPAR0) << "  ";
  for (int j = 1; j < var.Size(); j++)
  {
    out.width(13);
    out << par(var(j) - VarPAR0) << "  ";
  }
  out.width(2);
  out << ustab << " ";
  out.width(2);
  out << it;
}

#endif
