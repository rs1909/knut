// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002 - 2008 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef BASEPOINT_H
#define BASEPOINT_H

#include "matrix.h"
#include "hypermatrix.h"
#include "pointtype.h"
#include "basecolloc.h"
#include "multipliers.h"
#include <cfloat>

class System;

class BasePoint
{
  public:
    BasePoint(System& sys, const Array1D<Eqn>& eqn_, const Array1D<Var>& var_, 
          const int solsize, const int nz_jac_);
    virtual ~BasePoint();
    void    Reset(const Array1D<Eqn>& eqn_, const Array1D<Var>& var_);
    int     Refine(std::ostream& str = std::cout, bool adapt = false);
    int     Tangent(bool adapt = false);
    int     Continue(double ds, bool jacstep);

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
      for (int i = 0; (i < p.size()) && (i < par.size()); i++) par(i) = p(i);
    }
    inline Vector& getPar()
    {
      return par;
    }

    inline void    setCont(int p)
    {
      varMapCont(varMap.size()) = p;
    }
    inline int     getCont()
    {
      return varMapCont(varMap.size());
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
      // for( int j=0; j<colloc->Ndim(); j++ ) e += sol(j)*sol(j); // does not matter that headpoint or tailpoint, its periodic...
      return sqrt(basecolloc->Integrate(sol, sol));
    }
/// END BASE CLASS
  protected:
    virtual void Construct();
    virtual void Destruct();

    // should be moved to BaseColloc? Does it use any internal data?
    // this is a wrapper to collocation routines, so it can be here...
    virtual void Jacobian
    (
      HyperMatrix& AA, HyperVector& RHS, // output
      Vector& parPrev, Vector& par,      // parameters
      Vector& solPrev, Vector& sol,      // solution
      Array1D<int>&    varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    ) = 0;
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
    int          nz_jac;

    // solutions
    Vector       sol;
    Vector       par;
    HyperVector* xxDot;

    // solution in the continuation
    Vector       solNu;
    Vector       parNu;
    HyperVector* xxDotNu;

    // int the Newton iterations
    HyperVector* xx;     // -> sol,   --- this will be the solution in the iterations i.e., Dxx
    HyperVector* rhs;

    // store the Jacobian
    HyperMatrix* jac;

    // for collocation
    BaseColloc*  basecolloc;
};

inline void BasePoint::Update(HyperVector& X)
{
  sol += X.getV1();
  for (int i = 1; i < varMap.size(); i++) par(varMap(i)) += X.getV3()(i - 1);
}

inline void BasePoint::ContUpdate(HyperVector& X)
{
  solNu += X.getV1();
  for (int i = 1; i < varMapCont.size(); i++) parNu(varMapCont(i)) += X.getV3()(i - 1);
}

inline void BasePoint::AdaptUpdate(HyperVector& X)
{
  sol += X.getV1();
  for (int i = 1; i < varMapCont.size(); i++) par(varMapCont(i)) += X.getV3()(i - 1);
}

class mat4Data;
class PerSolColloc;

class PerSolPoint : public BasePoint
{
  public:
    PerSolPoint(System& sys, const Array1D<Eqn>& eqn_, const Array1D<Var>& var_, const int solsize, const int nz_jac_, const int nmul) 
      : BasePoint(sys, eqn_, var_, solsize, nz_jac_), mRe(nmul), mIm(nmul), nTrivMulLP(0), nTrivMulPD(0), nTrivMulNS(0) { }
    virtual ~PerSolPoint() {}
    
    virtual void Stability() = 0;
    virtual void SwitchTFLP(BranchSW type, double ds) = 0;   // switches branch with testFunct
    virtual void SwitchTFPD(double ds) = 0;   // switches branch with testFunct
    virtual void SwitchTFHB(double ds, std::ostream& out) = 0;   // switches branch with testFunct
    
    int StartTF(bool findangle, std::ostream& out);
    inline void    setSym(int n, int* sRe, int* sIm)
    {
      rotRe.init(n);
      rotIm.init(n);
      for (int i = 0; i < n; i++)
      {
        rotRe(i) = sRe[i];
        rotIm(i) = sIm[i];
      }
    }

    inline void    setSym(Array1D<int>& sRe, Array1D<int>& sIm)
    {
      P_ASSERT(sRe.size() == sIm.size());
      rotRe.init(sRe.size());
      rotIm.init(sRe.size());
      for (int i = 0; i < sRe.size(); i++)
      {
        rotRe(i) = sRe(i);
        rotIm(i) = sIm(i);
      }
    }

    int     UStab() { return unstableMultipliers(mRe, mIm, nTrivMulLP, nTrivMulPD, nTrivMulNS); }
    BifType  testBif() { return bifurcationType(mRe, mIm, nTrivMulLP, nTrivMulPD, nTrivMulNS); }
    void    clearStability() { mRe.clear(); mIm.clear(); }
    static inline double  Amplitude(const Vector& sol, int ndim, int ndeg, int nint)
    {
      double nrm = 0.0;
      for (int p = 0; p < ndim; ++p)
      {
        double min = DBL_MAX;
        double max = -DBL_MAX;
        for (int j = 0; j < nint; ++j)
        {
          for (int k = 0; k < ndeg; ++k)
          {
            if (min > sol(p + ndim*(k + ndeg*j))) min = sol(p + ndim*(k + ndeg*j));
            if (max < sol(p + ndim*(k + ndeg*j))) max = sol(p + ndim*(k + ndeg*j));
          }
        }
        nrm += (max - min) * (max - min);
      }
      return nrm;
    }
    inline double NormMX()
    {
      const int ndeg = persolcolloc->Nint();
      const int nint = persolcolloc->Nint();
      const int ndim = persolcolloc->Ndim();
      return Amplitude(sol, ndim, ndeg, nint);
    }
    
    void BinaryRead(mat4Data& data, int n);
    void BinaryWrite(mat4Data& data, int n);

  protected:
//    virtual void Construct();
//    virtual void Destruct();
    void         FillSol(System& sys_);

    // multipliers
    Vector       mRe;
    Vector       mIm;
    // number of trivial multipliers
    int   nTrivMulLP, nTrivMulPD, nTrivMulNS;

    // for the rotation phase conditions
    Array1D<int> rotRe;
    Array1D<int> rotIm;
    
    PerSolColloc* persolcolloc;
};

#endif /*BASEPOINT_H*/
