// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef NPOINT_H
#define NPOINT_H

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

class Point
{
  public:

    // constructor
    Point(System& sys, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg, int nmul = MULTIPLIERS, int nmat = BMATRICES);
    ~Point();

    // algorithms
    void    Reset(Array1D<Eqn>& eqn_, Array1D<Var>& var_);
    void    SwitchTFLP(double ds);   // switches branch with testFunct
    void    SwitchTFPD(double ds);   // switches branch with testFunct
    void    SwitchTFHB(double ds);   // switches branch with testFunct
    void    SwitchTRSol(Vector& Sol, const Vector& mshint, const Vector& mshdeg)
    {
      colloc.Export(Sol, mshint, mshdeg, sol);
    } // starting data for tori: solution
    void    SwitchTFTRTan(Vector& Re, Vector& Im, double& alpha, const Vector& mshint, const Vector& mshdeg);   // starting data for tori: tangent WITH testFunct
    int     StartTF(Eqn FN, std::ostream& out = std::cout);
    int     Refine(std::ostream& str = std::cout, bool adapt = false);
    int     Tangent(bool adapt = false);
    int     Continue(double ds, bool jacstep);
    void    Stability();

    // supplementary
    inline void    setSol(Vector& s)
    {
      sol = s;
    }
    inline Vector& getSol()
    {
      return sol;
    }

    inline void    setPar(Vector& p)
    {
      for (int i = 0; (i < p.Size()) && (i < par.Size()); i++) par(i) = p(i);
    }
    inline Vector& getPar()
    {
      return par;
    }

    inline void    setCont(int p)
    {
      varMapCont(varMap.Size()) = p;
    }
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

    inline void    setRefIter(int i)
    {
      RefIter = i;
    }
    inline void    setContIter(int i)
    {
      ContIter = i;
    }
    inline void    setKernIter(int i)
    {
      KernIter = i;
    }
    inline void    setRefEps(double d)
    {
      RefEps = d;
    }
    inline void    setContEps(double d)
    {
      ContEps = d;
    }
    inline void    setKernEps(double d)
    {
      KernEps = d;
    }

    inline double  Norm()
    {
      // double e=0;
      // for( int j=0; j<colloc.Ndim(); j++ ) e += sol(j)*sol(j); // does not matter that headpoint or tailpoint, its periodic...
      return sqrt(colloc.Integrate(sol, sol));
    }
    inline double  NormMX()
    {
      double max = 0.0, min = 1.0e32;
      for (int i = 0; i < colloc.Nint()*colloc.Ndeg() + 1; i++)
      {
        double e = 0.0;
        for (int j = 0; j < colloc.Ndim(); j++) e += sqrt(sol(colloc.Ndim() * i + j) * sol(colloc.Ndim() * i + j));
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
      for (int i = 0; i < colloc.Nint()*colloc.Ndeg() + 1; i++)
      {
        for (int j = 0; j < colloc.Ndim(); j++) ff << sol(colloc.Ndim()*i + j) << "\t";
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
    void    Dispose();
    void    FillSol(System& sys_);

    // internal member-functions
//   inline double SolNorm( Vector& sol, Vector& par );

    void Jacobian
    (
      HyperMatrix& AA, HyperVector& RHS,                      // output
      Vector& parPrev, Vector& par,                           // parameters
      Vector& solPrev, Vector& sol, Array3D<double>& solData, // solution
      Array1D<int>&    varMap,                                // contains the variables. If cont => contains the P0 too.
      double ds, bool cont                                    // ds stepsize, cont: true if continuation
    );

    inline void   Update(HyperVector& X);                            // sol,   par              += X
    inline void   ContUpdate(HyperVector& X);                        // solNu, parNu, parNu(cp) += X
    inline void   AdaptUpdate(HyperVector& X);                        // sol,   par,   par(cp)   += X

    // convergence parameters
    double       RefEps;
    double       ContEps;
    double       KernEps;
    int          RefIter;
    int          ContIter;
    int          KernIter;

    // variables and equations
    Array1D<Var> var;
    Array1D<Eqn> eqn;
    Array1D<int> varMap;
    Array1D<int> varMapCont;

    int          dim1;
    int          dim3;

    // solutions
    Vector       sol;
    Vector       par;

    // solution in the continuation
    Vector       solNu;
    Vector       parNu;

    // tangents
    HyperVector* xxDot;  // -> solDot, qqDot, parDot

    // multipliers
    Vector       mRe;
    Vector       mIm;
    // number of trivial multipliers
    int   nTrivMulLP, nTrivMulPD, nTrivMulNS;

    // for the rotation phase conditions
    Array1D<int> rotRe;
    Array1D<int> rotIm;

    // int the Newton iterations
    HyperVector* xx;     // -> sol,   --- this will be the solution in the iterations i.e., Dxx
    HyperVector* rhs;

    // store the Jacobian
    HyperMatrix* jac;

    // for collocation
    NColloc      colloc;

    // the solutions interpolated values
    Array3D<double> solData;

    // stability matrix
    StabMatrix   jacStab;

    // test functionals
    Array1D<baseTestFunct*> testFunct;
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
