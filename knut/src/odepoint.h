// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef ODEPOINT_H
#define ODEPOINT_H

#include "pointtype.h"
#include "basepoint.h"

class KNSystem;
class KNVector;
class KNBlockVector;
class KNOdeBvpCollocation;
template< class T > class KNArray1D;

class KNOdePeriodicSolution : public KNAbstractPeriodicSolution
{
  public:

    KNOdePeriodicSolution(KNAbstractContinuation* cnt, KNSystem& sys, 
      KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_, size_t nint, size_t ndeg);
    ~KNOdePeriodicSolution();
    
    virtual void Stability(bool init);
    virtual void SwitchTFLP(BranchSW , double ) { P_MESSAGE1("Not implemented"); }     // switches branch with testFunct
    virtual void SwitchTFPD(double ) { P_MESSAGE1("Not implemented"); }                    // switches branch with testFunct
    virtual void SwitchTFHB(double ) { P_MESSAGE1("Not implemented"); } // switches branch with testFunct
  private:
    virtual void jacobian(
      KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
      KNVector& parPrev, KNVector& par,      // parameters
      KNVector& solPrev, KNVector& sol,      // solution
      KNArray1D<size_t>& varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    );
    virtual void postProcess() {}
    
    KNOdeBvpCollocation* colloc;
    KNLuSparseMatrix     jacStab;
    KNMatrix     matrixInitialCondition;
    KNMatrix     matrixSolution;
    KNMatrix     monodromyMatrix;
};

#endif // ODEPOINT_H
