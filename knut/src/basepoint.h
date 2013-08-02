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


class KNDataFile;
class KNSystem;
class KNAbstractContinuation;

enum class IterateTangent { no, yes };

class KNAbstractPoint
{
  public:
    // members
    KNAbstractPoint(); // not defined
    KNAbstractPoint(KNAbstractContinuation* cnt, KNSystem& sys, 
          const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_, 
          const size_t solsize, const size_t nz_jac_);
    virtual ~KNAbstractPoint();
    void    reset(const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_);
    void    adapt()
    {
      basecolloc->meshAdapt(solNu, sol, xxDotNu->getV1(), xxDot->getV1());
      sol = solNu;
      xxDot->getV1() = xxDotNu->getV1();
    }
    size_t  refine(bool adapt = false);
    size_t  tangent(bool adapt = false);
    size_t  nextStep(double ds, double& angle, const IterateTangent jacstep);

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
      for (size_t i = 0; (i < p.size()) && (i < par.size()); i++) par(i) = p(i);
    }
    inline KNVector& getPar()
    {
      return par;
    }

    inline void    setCont(Var p)
    {
      varMapCont(varMap.size()) = VarToIndex(p, npar);
    }
    inline size_t  getCont()
    {
      return varMapCont(varMap.size());
    }
    inline void    setRefIter(size_t i)
    {
      RefIter = i;
    }
    inline void    setContIter(size_t i)
    {
      ContIter = i;
    }
    inline void    setKernIter(size_t i)
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
    inline void    setContCurvature(double d)
    {
      ContCurvature = d;
    }

    inline double  norm()
    {
      // double e=0;
      // for( int j=0; j<colloc->nDim(); j++ ) e += sol(j)*sol(j); // does not matter that headpoint or tailpoint, its periodic...
      return sqrt(basecolloc->integrate(sol, sol));
    }
    virtual inline double  normMX() { return norm(); }
    
    std::ostream& outStream();
    void printStream();
    void setNotifyObject(KNAbstractContinuation* cnt);
    
    /// I/O to REIMPLEMENT
    virtual void BinaryRead(KNDataFile& data, size_t n) = 0;
    virtual void BinaryWrite(KNDataFile& data, BifType bif, size_t n) = 0;

/// END BASE CLASS
  protected:
    virtual void construct();
    virtual void destruct();

    // should be moved to KNAbstractCollocation? Does it use any internal data?
    // this is a wrapper to collocation routines, so it can be here...
    // BEGIN --- TO REIMPLEMENT
    virtual void jacobian(
      KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
      KNVector& parPrev, KNVector& par,      // parameters
      KNVector& solPrev, KNVector& sol,      // solution
      KNArray1D<size_t>& varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    ) = 0;
    virtual void  postProcess() = 0;
    // END --- TO REIMPLEMENT
    inline void   update(KNBlockVector& X);                            // sol,   par              += X
    inline void   updateWithCp(KNBlockVector& X);                        // solNu, parNu, parNu(cp) += X
    inline void   updateWithAdaptation(KNBlockVector& X);                        // sol,   par,   par(cp)   += X

    // convergence parameters
    double       RefEps;
    double       ContEps;
    double       KernEps;
    size_t       RefIter;
    size_t       ContIter;
    size_t       KernIter;
    double       ContCurvature;

    // variables and equations
    KNArray1D<Var> var;
    KNArray1D<Eqn> eqn;
    KNArray1D<size_t> varMap;
    KNArray1D<size_t> varMapCont;

    size_t       dim1;
    size_t       dim3;
    size_t       nz_jac;
    const size_t npar;

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
  for (size_t i = 1; i < varMap.size(); i++) par(varMap(i)) += X.getV3()(i - 1);
}

inline void KNAbstractPoint::updateWithCp(KNBlockVector& X)
{
  solNu += X.getV1();
  for (size_t i = 1; i < varMapCont.size(); i++) parNu(varMapCont(i)) += X.getV3()(i - 1);
}

inline void KNAbstractPoint::updateWithAdaptation(KNBlockVector& X)
{
  sol += X.getV1();
  for (size_t i = 1; i < varMapCont.size(); i++) par(varMapCont(i)) += X.getV3()(i - 1);
}

class KNDataFile;
class KNAbstractBvpCollocation;

class KNAbstractPeriodicSolution : public KNAbstractPoint
{
  public:
    KNAbstractPeriodicSolution(KNAbstractContinuation* cnt, KNSystem& sys, 
      const KNArray1D<Eqn>& eqn_, const KNArray1D<Var>& var_, 
      const size_t solsize, const size_t nz_jac_, const size_t nmul) 
      : KNAbstractPoint(cnt, sys, eqn_, var_, solsize, nz_jac_), mRe(nmul), mIm(nmul), mRePrev(nmul), mImPrev(nmul),
        nTrivMulLP(0), nTrivMulPD(0), nTrivMulNS(0) { }
    virtual ~KNAbstractPeriodicSolution() {}
    
    virtual void Stability(bool init) = 0;
    virtual void SwitchTFLP(BranchSW type, double ds) = 0;   // switches branch with testFunct
    virtual void SwitchTFPD(double ds) = 0;   // switches branch with testFunct
    virtual void SwitchTFHB(double ds) = 0;   // switches branch with testFunct
    
    virtual void findAngle();
    inline  void setSym(size_t n, size_t* sRe, size_t* sIm)
    {
      rotRe.init(n);
      rotIm.init(n);
      for (size_t i = 0; i < n; i++)
      {
        rotRe(i) = sRe[i];
        rotIm(i) = sIm[i];
      }
    }

    inline void    setSym(KNArray1D<size_t>& sRe, KNArray1D<size_t>& sIm)
    {
      P_ASSERT(sRe.size() == sIm.size());
      rotRe.init(sRe.size());
      rotIm.init(sRe.size());
      for (size_t i = 0; i < sRe.size(); i++)
      {
        rotRe(i) = sRe(i);
        rotIm(i) = sIm(i);
      }
    }

    size_t  UStab(size_t pt = 0) { return unstableMultipliers(mRe, mIm, nTrivMulLP, nTrivMulPD, nTrivMulNS, pt); }
    BifType testBif(size_t pt = 0)
    {
      return bifurcationType(mRePrev, mImPrev, mRe, mIm, nTrivMulLP, nTrivMulPD, nTrivMulNS, pt-1, pt);
    }
    void    storeMultiplier() { mRePrev = mRe; mImPrev = mIm; }
    void    clearStability() { mRe.clear(); mIm.clear(); }
    static inline double  Amplitude(const KNVector& sol, size_t ndim, size_t ndeg, size_t nint)
    {
      double nrm = 0.0;
      for (size_t p = 0; p < ndim; ++p)
      {
        double min = DBL_MAX;
        double max = -DBL_MAX;
        for (size_t j = 0; j < nint; ++j)
        {
          for (size_t k = 0; k < ndeg; ++k)
          {
            if (min > sol(p + ndim*(k + ndeg*j))) min = sol(p + ndim*(k + ndeg*j));
            if (max < sol(p + ndim*(k + ndeg*j))) max = sol(p + ndim*(k + ndeg*j));
          }
        }
        nrm += (max - min) * (max - min);
      }
      return nrm;
    }
    // redefinition of the abstract version
    inline double normMX()
    {
      const size_t ndeg = persolcolloc->nInt();
      const size_t nint = persolcolloc->nInt();
      const size_t ndim = persolcolloc->nDim();
      return Amplitude(sol, ndim, ndeg, nint);
    }
    
    void BinaryRead(KNDataFile& data, size_t n);
    void BinaryWrite(KNDataFile& data, BifType bif, size_t n);

  protected:
//    virtual void construct();
//    virtual void destruct();
    void         FillSol(KNSystem& sys_);

    // multipliers
    KNVector       mRe;
    KNVector       mIm;
    // previous multipliers
    KNVector       mRePrev;
    KNVector       mImPrev;
    // number of trivial multipliers
    size_t   nTrivMulLP, nTrivMulPD, nTrivMulNS;

    // for the rotation phase conditions
    KNArray1D<size_t> rotRe;
    KNArray1D<size_t> rotIm;
    
    KNAbstractBvpCollocation* persolcolloc;
};

#endif /*BASEPOINT_H*/
