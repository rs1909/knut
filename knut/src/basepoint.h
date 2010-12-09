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

class KNSystem;
class KNAbstractContinuation;

class KNAbstractPoint
{
  public:
    // members
    KNAbstractPoint(); // not defined
    KNAbstractPoint(KNAbstractContinuation* cnt, KNSystem& sys, 
          const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_, 
          const int solsize, const int nz_jac_);
    virtual ~KNAbstractPoint();
    void    reset(const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_);
    int     refine(bool adapt = false);
    int     tangent(bool adapt = false);
    int     nextStep(double ds, double& angle, bool jacstep);

    // supplementary
    inline void    setSol(KNVector& s)
    {
      sol = s;
    }
    inline KNVector& getSol()
    {
      return sol;
    }

    inline void    setPar(KNVector& p)
    {
      for (int i = 0; (i < p.size()) && (i < par.size()); i++) par(i) = p(i);
    }
    inline KNVector& getPar()
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

    inline double  norm()
    {
      // double e=0;
      // for( int j=0; j<colloc->nDim(); j++ ) e += sol(j)*sol(j); // does not matter that headpoint or tailpoint, its periodic...
      return sqrt(basecolloc->integrate(sol, sol));
    }
    std::ostream& outStream();
    void printStream();
    void printClearLastLine();
    void setNotifyObject(KNAbstractContinuation* cnt);
/// END BASE CLASS
  protected:
    virtual void construct();
    virtual void destruct();

    // should be moved to KNAbstractCollocation? Does it use any internal data?
    // this is a wrapper to collocation routines, so it can be here...
    virtual void jacobian(
      KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
      KNVector& parPrev, KNVector& par,      // parameters
      KNVector& solPrev, KNVector& sol,      // solution
      KNArray1D<int>&    varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    ) = 0;
    inline void   update(KNBlockVector& X);                            // sol,   par              += X
    inline void   updateWithCp(KNBlockVector& X);                        // solNu, parNu, parNu(cp) += X
    inline void   updateWithAdaptation(KNBlockVector& X);                        // sol,   par,   par(cp)   += X

    // convergence parameters
    double       RefEps;
    double       ContEps;
    double       KernEps;
    int          RefIter;
    int          ContIter;
    int          KernIter;

    // variables and equations
    KNArray1D<Var> var;
    KNArray1D<Eqn> eqn;
    KNArray1D<int> varMap;
    KNArray1D<int> varMapCont;

    int          dim1;
    int          dim3;
    int          nz_jac;

    // solutions
    KNVector       sol;
    KNVector       par;
    KNBlockVector* xxDot;

    // solution in the continuation
    KNVector       solNu;
    KNVector       parNu;
    KNBlockVector* xxDotNu;

    // int the Newton iterations
    KNBlockVector* xx;     // -> sol,   --- this will be the solution in the iterations i.e., Dxx
    KNBlockVector* rhs;

    // store the Jacobian
    KNSparseBlockMatrix* jac;

    // for collocation
    KNAbstractCollocation*  basecolloc;
    
    KNAbstractContinuation* notifyObject;
};

inline void KNAbstractPoint::update(KNBlockVector& X)
{
  sol += X.getV1();
  for (int i = 1; i < varMap.size(); i++) par(varMap(i)) += X.getV3()(i - 1);
}

inline void KNAbstractPoint::updateWithCp(KNBlockVector& X)
{
  solNu += X.getV1();
  for (int i = 1; i < varMapCont.size(); i++) parNu(varMapCont(i)) += X.getV3()(i - 1);
}

inline void KNAbstractPoint::updateWithAdaptation(KNBlockVector& X)
{
  sol += X.getV1();
  for (int i = 1; i < varMapCont.size(); i++) par(varMapCont(i)) += X.getV3()(i - 1);
}

class KNDataFile;
class KNAbstractBvpCollocation;

class KNAbstractPeriodicSolution : public KNAbstractPoint
{
  public:
    KNAbstractPeriodicSolution(KNAbstractContinuation* cnt, KNSystem& sys, 
      const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_, 
      const int solsize, const int nz_jac_, const int nmul) 
      : KNAbstractPoint(cnt, sys, eqn_, var_, solsize, nz_jac_), mRe(nmul), mIm(nmul), 
        nTrivMulLP(0), nTrivMulPD(0), nTrivMulNS(0) { }
    virtual ~KNAbstractPeriodicSolution() {}
    
    virtual void Stability(bool init) = 0;
    virtual void SwitchTFLP(BranchSW type, double ds) = 0;   // switches branch with testFunct
    virtual void SwitchTFPD(double ds) = 0;   // switches branch with testFunct
    virtual void SwitchTFHB(double ds) = 0;   // switches branch with testFunct
    
    virtual void findAngle();
    inline  void setSym(int n, int* sRe, int* sIm)
    {
      rotRe.init(n);
      rotIm.init(n);
      for (int i = 0; i < n; i++)
      {
        rotRe(i) = sRe[i];
        rotIm(i) = sIm[i];
      }
    }

    inline void    setSym(KNArray1D<int>& sRe, KNArray1D<int>& sIm)
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
    static inline double  Amplitude(const KNVector& sol, int ndim, int ndeg, int nint)
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
      const int ndeg = persolcolloc->nInt();
      const int nint = persolcolloc->nInt();
      const int ndim = persolcolloc->nDim();
      return Amplitude(sol, ndim, ndeg, nint);
    }
    
    void BinaryRead(KNDataFile& data, int n);
    void BinaryWrite(KNDataFile& data, int n);

  protected:
//    virtual void construct();
//    virtual void destruct();
    void         FillSol(KNSystem& sys_);

    // multipliers
    KNVector       mRe;
    KNVector       mIm;
    // number of trivial multipliers
    int   nTrivMulLP, nTrivMulPD, nTrivMulNS;

    // for the rotation phase conditions
    KNArray1D<int> rotRe;
    KNArray1D<int> rotIm;
    
    KNAbstractBvpCollocation* persolcolloc;
};

#endif /*BASEPOINT_H*/
